"""Add allele-specific FACETS annotations from a vcf2maf MAF to VCF INFO."""
import argparse
import csv
import math
import os
from pathlib import Path
import re
import sys
import tempfile

import vcfpy

from vatools.utils import open_maybe_gz, write_record
from vatools.vcf_info_annotator import column_mappings_type, coerce_value


def facets_mappings():
    columns = [
        ('tcn', 'Integer', 'Total copy number'),
        ('lcn', 'Integer', 'Minor allele copy number'),
        ('cf', 'Float', 'Cellular fraction of the copy number segment'),
        ('purity', 'Float', 'Tumor purity'),
        ('t_var_freq', 'Float', 'Observed tumor variant allele frequency'),
        ('expected_alt_copies', 'Float', 'Expected mutant allele copies'),
    ]
    for prefix, model in [('ccf_Mcopies', 'major allele copies'),
                          ('ccf_1copy', 'one mutant copy'),
                          ('ccf_expected_copies', 'expected mutant copies')]:
        for suffix, detail in [('', 'estimate'), ('_lower', 'lower bound'),
                               ('_upper', 'upper bound'), ('_prob95', 'prob95'),
                               ('_prob90', 'prob90')]:
            columns.append((prefix + suffix, 'Float', 'CCF {} assuming {}'.format(detail, model)))
    columns.append(('clonality', 'String', 'Mutation clonality classification'))
    return [dict(source_col=col, info_field=col, type=kind,
                 description='FACETS ' + description, source='facets-suite', version=None)
            for col, kind, description in columns]


def normalize_chromosome(chrom):
    chrom = chrom[3:] if chrom.startswith('chr') else chrom
    return 'MT' if chrom in ('M', 'MT') else chrom


def vcf_to_maf_key(chrom, pos, ref, alt):
    """Reproduce vcf2maf 1.6.22 prefix removal and insertion coordinates.

    Keep the original VCF position as well as the MAF representation. This is
    intentionally not reference-based left alignment or coordinate-only matching.
    """
    ref, alt = ref.upper(), alt.upper()
    if not re.fullmatch('[ACGTN]+', ref) or not re.fullmatch('[ACGTN]+', alt) or ref == alt:
        return None  # Symbolic/breakend/spanning-deletion alleles have no small-variant key.
    start = pos
    while ref and alt and ref[0] == alt[0] and ref != alt:
        ref, alt = ref[1:], alt[1:]
        start += 1
    if not ref:
        start -= 1
    return normalize_chromosome(chrom), pos, start, ref or '-', alt or '-'


def annotation_value(value, kind):
    if value.strip().lower() in ('', '.', 'na', 'nan', 'inf', '+inf', '-inf'):
        return None
    value = coerce_value(value, kind)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if kind in ('String', 'Character'):
        if re.search(r'[\s;,=]', value):
            raise ValueError('String annotations cannot contain whitespace, semicolons, commas or equals signs')
        if kind == 'Character' and len(value) != 1:
            raise ValueError('Character annotations must have length 1')
    return value


def load_maf(path, mappings, sample=None):
    required = {'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele',
                'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 'vcf_pos'}
    required.update(m['source_col'] for m in mappings)
    values, samples = {}, set()
    with open_maybe_gz(path) as stream:
        reader = csv.DictReader((line for line in stream if line.strip() and not line.startswith('#')),
                                delimiter='\t')
        header = reader.fieldnames or []
        if len(header) != len(set(header)):
            raise ValueError('MAF contains duplicate column names')
        missing = required - set(header)
        if missing:
            raise ValueError('MAF missing columns: ' + ', '.join(sorted(missing)))
        for row_number, row in enumerate(reader, 1):
            if None in row or any(value is None for value in row.values()):
                raise ValueError('Malformed MAF data row {}'.format(row_number))
            barcode = row['Tumor_Sample_Barcode']
            if sample is not None and barcode != sample:
                continue
            if not barcode or barcode in ('.', 'NA'):
                raise ValueError('MAF tumor sample is missing')
            samples.add(barcode)
            if len(samples) > 1:
                raise ValueError('MAF contains multiple tumors; select one with --tumor-sample-name')
            try:
                pos, start, end = (int(row[c]) for c in ('vcf_pos', 'Start_Position', 'End_Position'))
                ref, alt = row['Reference_Allele'].upper(), row['Tumor_Seq_Allele2'].upper()
                if not all(re.fullmatch(r'(?:[ACGTN]+|-)', a) for a in (ref, alt)) or ref == alt:
                    raise ValueError('Expected distinct sequence alleles or a MAF indel dash')
                expected_end = start + 1 if ref == '-' else start + len(ref) - 1
                if min(pos, start) < 1 or end != expected_end or start < pos:
                    raise ValueError('Inconsistent vcf2maf coordinates')
                if row.get('Strand', '+') != '+':
                    raise ValueError('Only positive-strand MAF alleles are supported')
                key = normalize_chromosome(row['Chromosome']), pos, start, ref, alt
                if key in values:
                    raise ValueError('Duplicate MAF variant key {}'.format(key))
                values[key] = {m['info_field']: annotation_value(row[m['source_col']], m['type'])
                               for m in mappings}
            except ValueError as error:
                raise ValueError('MAF data row {}: {}'.format(row_number, error)) from error
    if not samples:
        raise ValueError('No MAF rows found for the selected tumor')
    return values, next(iter(samples))


def define_parser():
    parser = argparse.ArgumentParser('vcf-facets-maf-annotator', description=__doc__)
    parser.add_argument('input_vcf', help='Original VCF or a descendant preserving POS/REF/ALT; may be gzipped')
    parser.add_argument('maf', help='FACETS-annotated vcf2maf MAF, including vcf_pos; may be gzipped')
    parser.add_argument('-o', '--output-vcf', help='Output path; default: <input>.facets.vcf[.gz]')
    parser.add_argument('-s', '--tumor-sample-name', help='Select a MAF tumor; otherwise require exactly one')
    parser.add_argument('-m', '--column-mappings', type=column_mappings_type(),
                        help='Override the 22 FACETS mappings with source_col:info_field:type:description'
                             '[:source[:version]],...; Flag is not supported for per-ALT fields')
    parser.add_argument('-w', '--overwrite', action='store_true',
                        help='Replace mapped INFO headers and clear their old values on ALL records')
    parser.add_argument('--require-all-maf-matched', action='store_true',
                        help='Fail without publishing output if any selected MAF row has no VCF allele match')
    return parser


def main(args_input=None):
    parser = define_parser()
    args = parser.parse_args(args_input)
    mappings = args.column_mappings if args.column_mappings is not None else facets_mappings()
    ids = [m['info_field'] for m in mappings]
    if len(ids) != len(set(ids)) or any(not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_.]*', name) for name in ids):
        parser.error('Mapped INFO IDs must be valid and unique')
    if any(m['type'] == 'Flag' for m in mappings):
        parser.error('Flag cannot represent per-ALT values (Number=A)')
    output = args.output_vcf
    if output is None:
        stem = re.sub(r'\.vcf(?:\.gz)?$', '', args.input_vcf)
        output = stem + '.facets.vcf' + ('.gz' if args.input_vcf.endswith('.gz') else '')
    for source in (args.input_vcf, args.maf):
        if Path(output).resolve() == Path(source).resolve() or (os.path.exists(output) and os.path.samefile(output, source)):
            parser.error('Output must differ from both input files')
    values, sample = load_maf(args.maf, mappings, args.tumor_sample_name)
    matched, seen = set(), set()
    counts = dict(records=0, matched_records=0, matched_alleles=0, unmatched_alleles=0)
    temporary = None
    try:
        with vcfpy.Reader.from_path(args.input_vcf) as reader:
            if sample not in reader.header.samples.names:
                raise ValueError('MAF tumor {} is not a sample in the VCF'.format(sample))
            header = reader.header.copy()
            conflicts = set(ids) & set(header.info_ids())
            if conflicts and not args.overwrite:
                raise ValueError('INFO fields already exist: {}; use --overwrite to replace them'.format(', '.join(sorted(conflicts))))
            if conflicts:
                header.lines = [line for line in header.lines
                                if not (line.key == 'INFO' and line.mapping['ID'] in conflicts)]
                header._indices = header._build_indices()
            for mapping in mappings:
                info = dict(ID=mapping['info_field'], Number='A', Type=mapping['type'],
                            Description=mapping['description'] + ' (tumor: ' + sample + ')')
                if mapping['source']:
                    info['Source'] = mapping['source']
                if mapping['version']:
                    info['Version'] = mapping['version']
                header.add_info_line(info)
            suffix = '.vcf.gz' if output.endswith('.gz') else '.vcf'
            fd, temporary = tempfile.mkstemp(prefix='.vcf-maf-', suffix=suffix,
                                              dir=Path(output).resolve().parent)
            os.close(fd)
            with vcfpy.Writer.from_path(temporary, header) as writer:
                for record in reader:
                    counts['records'] += 1
                    rows = []
                    for alt in record.ALT:
                        key = vcf_to_maf_key(record.CHROM, record.POS, record.REF, alt.serialize())
                        row = values.get(key)
                        if row is not None:
                            if key in seen:
                                raise ValueError('Multiple VCF alleles match the same MAF key: {}'.format(key))
                            seen.add(key)
                            matched.add(key)
                            counts['matched_alleles'] += 1
                        else:
                            counts['unmatched_alleles'] += 1
                        rows.append(row)
                    if any(row is not None for row in rows):
                        counts['matched_records'] += 1
                    for field in ids:
                        if args.overwrite:
                            record.INFO.pop(field, None)
                        annotations = [row[field] if row is not None else None for row in rows]
                        if any(value is not None for value in annotations):
                            record.INFO[field] = annotations
                    write_record(record, writer)
        counts['unmatched_maf_rows'] = len(values) - len(matched)
        if args.require_all_maf_matched and counts['unmatched_maf_rows']:
            raise ValueError('{} MAF rows have no VCF allele match'.format(counts['unmatched_maf_rows']))
        if not matched:
            print('WARNING: no MAF variants matched the VCF; check sample, assembly and allele representation', file=sys.stderr)
        os.replace(temporary, output)
        temporary = None
    finally:
        if temporary is not None:
            os.unlink(temporary)
    print('Tumor {}: {}'.format(sample, ', '.join('{}={}'.format(k, v) for k, v in counts.items())), file=sys.stderr)


if __name__ == '__main__':
    main()
