import argparse
import sys
import os
import re
import vcfpy
import tempfile
import csv
from collections import OrderedDict, namedtuple
import logging
from vatools.utils import open_maybe_gz, write_record

# parse the input params
def define_parser():
    parser = argparse.ArgumentParser(
        'vcf-readcount-annotator',
        description="A tool that will add the data from bam-readcount files to the VCF sample column."
    )
    parser.add_argument(
        "input_vcf",
        help="A VCF file"
    )
    parser.add_argument(
        "bam_readcount_file",
        help="A bam-readcount output file",
    )
    parser.add_argument(
        "data_type",
        choices=['DNA', 'RNA'],
        help="The type of data in the bam_readcount_file. "
            +"If `DNA` is chosen, the readcounts will be written to the AD, AF, and DP fields. "
            +"If `RNA` is chosen, the readcounts will be written to the RAD, RAF, and RDP fields."
    )
    parser.add_argument(
        "-s", "--sample-name",
        help="If the input_vcf contains multiple samples, the name of the sample to annotate."
    )
    parser.add_argument(
        "-o", "--output-vcf",
        help="Path to write the output VCF file. If not provided, the output VCF file will be "
            +"written next to the input VCF file with a .readcount.vcf file ending."
    )
    parser.add_argument(
        "-t", "--variant-type",
        help="The type of variant to process. `snv` will only annotate SNVs. `indel` will only annotate InDels. `all` will annotate all variant types. `snv` and `indel` mode currently do not support multi-allelic VCF entries that contain both SNVs and InDels. It is recommended to split multi-allelic sites before running in `snv` or `indel` mode.",
        choices=['snv', 'indel', 'all'],
        default='all'
    )
    extra = parser.add_argument_group('extra bam-readcount fields')
    extra.add_argument(
        '-a', '--all-fields', action='store_true', default=False,
        help='Append all extra bam-readcount fields to the output.'
    )
    extra.add_argument(
        '-q', '--avg-mapping-quality', action='store_true', default=False,
        help='Append avg mapping quality of variant-supporting reads (FORMAT tag: VAMQ).'
    )
    extra.add_argument(
        '-b', '--avg-basequality', action='store_true', default=False,
        help='Append avg base quality of variant-supporting reads (FORMAT tag: VABQ).'
    )
    extra.add_argument(
        '-e', '--avg-se-mapping-quality', action='store_true', default=False,
        help='Append avg SE mapping quality of variant-supporting reads (FORMAT tag: VASEMQ).'
    )
    extra.add_argument(
        '-r', '--strand-counts', action='store_true', default=False,
        help='Append ref and var forward/reverse strand read counts (FORMAT tags: ADF, ADR). '
             'In DNA mode ADF/ADR are already written by default, so this flag is a no-op.'
    )
    extra.add_argument(
        '-f', '--avg-pos-fraction', action='store_true', default=False,
        help='Append avg position of variant reads as fraction of read length (FORMAT tag: VAPF).'
    )
    extra.add_argument(
        '-m', '--avg-mismatches', action='store_true', default=False,
        help='Append avg mismatches per variant-supporting read as fraction (FORMAT tag: VAMF).'
    )
    extra.add_argument(
        '-k', '--sum-mismatch-qual', action='store_true', default=False,
        help='Append avg sum of mismatch base qualities for variant reads (FORMAT tag: VAMQS).'
    )
    extra.add_argument(
        '-2', '--num-q2-reads', action='store_true', default=False,
        help='Append number of variant-supporting reads containing a Q2 base (FORMAT tag: VAQ2).'
    )
    extra.add_argument(
        '-d', '--avg-q2-distance', action='store_true', default=False,
        help='Append avg distance to Q2 start in Q2-containing reads (FORMAT tag: VAQD).'
    )
    extra.add_argument(
        '-c', '--avg-clipped-length', action='store_true', default=False,
        help='Append avg clipped read length for variant-supporting reads (FORMAT tag: VACL).'
    )
    extra.add_argument(
        '-3', '--avg-3p-distance', action='store_true', default=False,
        help="Append avg distance to effective 3' end for variant reads (FORMAT tag: VA3P)."
    )
    return parser

# create an object to hold all of the accessory information about the bam-readcount columns:
#   brct_col: column name in bam-readcount per-base output (None for composite fields)
#   arg_name: argument name
#   tag:      VCF FORMAT tag, or 'strand_counts' for the ADF/ADR pair
#   vcf_type, number, desc: VCF header values (empty strings for the strand_counts)
ExtraField = namedtuple('ExtraField', ['brct_col', 'arg_name', 'tag', 'vcf_type', 'number', 'desc'])

EXTRA_FIELDS = [
    ExtraField('avg_mapping_quality',                 'avg_mapping_quality',    'VAMQ',         'Float',   '1', 'Avg mapping quality for variant-supporting reads'),
    ExtraField('avg_basequality',                     'avg_basequality',        'VABQ',         'Float',   '1', 'Avg base quality for variant-supporting reads'),
    ExtraField('avg_se_mapping_quality',              'avg_se_mapping_quality', 'VASEMQ',       'Float',   '1', 'Avg SE mapping quality for variant-supporting reads'),
    ExtraField(None,                                  'strand_counts',          'strand_counts','',        '',  ''),
    ExtraField('avg_pos_as_fraction',                 'avg_pos_fraction',       'VAPF',         'Float',   '1', 'Avg position of variant reads as fraction of read length'),
    ExtraField('avg_num_mismatches_as_fraction',      'avg_mismatches',         'VAMF',         'Float',   '1', 'Avg mismatches per variant-supporting read as fraction'),
    ExtraField('avg_sum_mismatch_qualities',          'sum_mismatch_qual',      'VAMQS',        'Float',   '1', 'Avg sum of mismatch base qualities for variant reads'),
    ExtraField('num_q2_containing_reads',             'num_q2_reads',           'VAQ2',         'Integer', '1', 'Number of variant reads containing a Q2 base'),
    ExtraField('avg_distance_to_q2_start_in_q2_reads','avg_q2_distance',        'VAQD',         'Float',   '1', 'Avg distance to Q2 start in Q2-containing reads'),
    ExtraField('avg_clipped_length',                  'avg_clipped_length',     'VACL',         'Float',   '1', 'Avg clipped read length for variant-supporting reads'),
    ExtraField('avg_distance_to_effective_3p_end',    'avg_3p_distance',        'VA3P',         'Float',   '1', "Avg distance to effective 3' end for variant reads"),
]

def get_requested_extra_fields(args):
    if getattr(args, 'all_fields', False):
        return list(EXTRA_FIELDS)
    return [f for f in EXTRA_FIELDS if getattr(args, f.arg_name, False)]

# parse the fields out for the detailed metrics
def parse_brct_field(brcts):
    # Column names match bam-readcount per-base output order
    quality_cols = (
        'avg_mapping_quality', 'avg_basequality', 'avg_se_mapping_quality',
        'num_plus_strand', 'num_minus_strand',
        'avg_pos_as_fraction', 'avg_num_mismatches_as_fraction',
        'avg_sum_mismatch_qualities', 'num_q2_containing_reads',
        'avg_distance_to_q2_start_in_q2_reads',
        'avg_clipped_length', 'avg_distance_to_effective_3p_end',
    )
    counts = {}
    forward_counts = {}
    reverse_counts = {}
    qualities = {}
    for brct in brcts:
        parts = brct.split(':')
        base = parts[0].upper()
        counts[base] = parts[1]
        qualities[base] = dict(zip(quality_cols, parts[2:]))
        forward_counts[base] = qualities[base].get('num_plus_strand', '0')
        reverse_counts[base] = qualities[base].get('num_minus_strand', '0')
    return counts, forward_counts, reverse_counts, qualities

# read the bam readcount file
def parse_bam_readcount_file(args):
    coverage = {}
    with open_maybe_gz(args.bam_readcount_file) as reader:
        coverage_tsv_reader = csv.reader(reader, delimiter='\t')
        for row in coverage_tsv_reader:
            chromosome     = row[0]
            position       = row[1]
            reference_base = row[2].upper()
            depth          = row[3]
            brct           = row[4:]
            counts, forward_counts, reverse_counts, qualities = parse_brct_field(brct)
            parsed_brct = {
                'counts': counts,
                'forward_counts': forward_counts,
                'reverse_counts': reverse_counts,
                'qualities': qualities,
                'depth': depth
            }
            if (chromosome, position, reference_base) in coverage and parsed_brct != coverage[(chromosome,position,reference_base)]:
                prev_brct = coverage[(chromosome, position, reference_base)]
                if prev_brct["depth"] == depth:
                    coverage[(chromosome, position, reference_base)] = {"depth" : depth}
                    logging.warning("Duplicate bam-readcount entry for chr {} pos {} ref {}. Both depths match, so this field will be written, but count and frequency fields will be skipped. Offending entries:\n{}\n{}".format(chromosome, position, reference_base, parsed_brct, prev_brct))
                else:
                    coverage[(chromosome, position, reference_base)] = [prev_brct, parsed_brct]
                    logging.warning("Duplicate bam-readcount entry for chr {} pos {} ref {}. Depths are discrepant, so neither entry will be included in the output vcf. Offending entries:\n{}\n{}".format(chromosome, position, reference_base, parsed_brct, prev_brct))
            else:
                coverage[(chromosome,position,reference_base)] = parsed_brct
    return coverage

def has_snv(entry):
    for alt in entry.ALT:
        alt = alt.serialize()
        ref = entry.REF
        if len(ref) == 1 and len(alt) == 1:
            return True
    return False

def has_indel(entry):
    for alt in entry.ALT:
        alt = alt.serialize()
        ref = entry.REF
        if len(ref) > 1 or len(alt) > 1:
            return True
    return False

def is_insertion(ref, alt):
    return len(alt) > len(ref)

def is_deletion(ref, alt):
    return len(alt) < len(ref)

def has_complex_variant(entry):
    if has_indel(entry):
        for alt in entry.ALT:
            alt = alt.serialize()
            ref = entry.REF
            if len(ref) == len(alt) or (len(ref) > 1 and len(alt) > 1):
                return True
    return False

def simplify_indel_allele(ref, alt):
    while len(ref)> 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref = ref[0:-1]
        alt = alt[0:-1]
    while len(ref)> 0 and len(alt) > 0 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
    return ref, alt

def calculate_coverage(ref, var):
    return ref + var

def calculate_vaf(var, depth):
    if int(depth) > 0:
        return format(var / int(depth), '.5f')
    else:
        return 0

def parse_to_bam_readcount(start, reference, alt, position):
    if len(alt) != len(reference):
        if is_deletion(reference, alt):
            bam_readcount_position = str(start + 2)
            (simplified_reference, simplified_alt) = simplify_indel_allele(reference, alt)
            ref_base = reference[1:2]
            var_base = '-' + simplified_reference
        elif is_insertion(reference, alt):
            bam_readcount_position = str(start)
            (simplified_reference, simplified_alt) = simplify_indel_allele(reference, alt)
            ref_base = reference
            var_base = '+' + simplified_alt
    else:
        bam_readcount_position = str(position)
        ref_base = reference
        var_base = alt
    return (bam_readcount_position, ref_base, var_base)

def create_vcf_reader(args):
    vcf_reader = vcfpy.Reader.from_path(args.input_vcf)
    if len(vcf_reader.header.samples.names) > 1:
        if args.sample_name is None:
            vcf_reader.close()
            raise Exception("ERROR: VCF {} contains more than one sample. Please use the -s option to specify which sample to annotate.".format(args.input_vcf))
        elif args.sample_name not in vcf_reader.header.samples.names:
            vcf_reader.close()
            raise Exception("ERROR: VCF {} does not contain a sample column for sample {}.".format(args.input_vcf, args.sample_name))
        else:
            sample_name = args.sample_name
    else:
        sample_name = vcf_reader.header.samples.names[0]
    return vcf_reader, sample_name

def create_vcf_writer(args, vcf_reader, extra_fields=None):
    if extra_fields is None:
        extra_fields = []
    if args.output_vcf:
        output_file = args.output_vcf
    else:
        (head, sep, tail) = args.input_vcf.rpartition('.vcf')
        output_file = ('').join([head, '.readcount.vcf', tail])
    new_header = vcfpy.Header(samples = vcf_reader.header.samples)
    if args.data_type == 'DNA':
        for line in vcf_reader.header.lines:
            if not (line.key == 'FORMAT' and line.id in ['DP', 'AD', 'ADF', 'ADR','AF']):
                new_header.add_line(line)
        new_header.add_format_line(OrderedDict([('ID', 'DP'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Read depth')]))
        new_header.add_format_line(OrderedDict([('ID', 'AD'), ('Number', 'R'), ('Type', 'Integer'), ('Description', 'Allelic depths for the ref and alt alleles in the order listed')]))
        new_header.add_format_line(OrderedDict([('ID', 'ADF'), ('Number', 'R'), ('Type', 'Integer'), ('Description', 'Allelic depths on the forward strand (high-quality bases)')]))
        new_header.add_format_line(OrderedDict([('ID', 'ADR'), ('Number', 'R'), ('Type', 'Integer'), ('Description', 'Allelic depths on the reverse strand (high-quality bases)')]))
        new_header.add_format_line(OrderedDict([('ID', 'AF'), ('Number', 'A'), ('Type', 'Float'), ('Description', 'Variant-allele frequency for the alt alleles')]))
    if args.data_type == 'RNA':
        for line in vcf_reader.header.lines:
            if not (line.key == 'FORMAT' and line.id in ['RDP', 'RAD', 'RADF', 'RADR', 'RAF']):
                new_header.add_line(line)
        new_header.add_format_line(OrderedDict([('ID', 'RDP'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'RNA Read depth')]))
        new_header.add_format_line(OrderedDict([('ID', 'RAD'), ('Number', 'R'), ('Type', 'Integer'), ('Description', 'RNA Allelic depths for the ref and alt alleles in the order listed')]))
        new_header.add_format_line(OrderedDict([('ID', 'RADF'), ('Number', 'R'), ('Type', 'Integer'), ('Description', 'RNA Allelic depths on the forward strand (high-quality bases)')]))
        new_header.add_format_line(OrderedDict([('ID', 'RADR'), ('Number', 'R'), ('Type', 'Integer'), ('Description', 'RNA Allelic depths on the reverse strand (high-quality bases)')]))
        new_header.add_format_line(OrderedDict([('ID', 'RAF'), ('Number', 'A'), ('Type', 'Float'), ('Description', 'RNA Variant-allele frequency for the alt alleles')]))
    for f in extra_fields:
        if f.tag == 'strand_counts':
            if args.data_type != 'DNA':
                new_header.add_format_line(OrderedDict([('ID', 'ADF'), ('Number', 'R'), ('Type', 'Integer'), ('Description', 'Allelic depths on the forward strand')]))
                new_header.add_format_line(OrderedDict([('ID', 'ADR'), ('Number', 'R'), ('Type', 'Integer'), ('Description', 'Allelic depths on the reverse strand')]))
        else:
            new_header.add_format_line(OrderedDict([
                ('ID', f.tag), ('Number', f.number), ('Type', f.vcf_type), ('Description', f.desc),
            ]))
    return vcfpy.Writer.from_path(output_file, new_header)

def add_format_value(entry, sample_name, field, value):
    if field not in entry.FORMAT:
        entry.FORMAT += [field]
    entry.call_for_sample[sample_name].data[field] = value

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    read_counts = parse_bam_readcount_file(args)
    extra_fields = get_requested_extra_fields(args)
    (vcf_reader, sample_name) = create_vcf_reader(args)
    vcf_writer = create_vcf_writer(args, vcf_reader, extra_fields)

    if args.data_type == 'DNA':
        depth_field = 'DP'
        count_field = 'AD'
        forward_count_field = 'ADF'
        reverse_count_field = 'ADR'
        frequency_field = 'AF'
    elif args.data_type == 'RNA':
        depth_field = 'RDP'
        count_field = 'RAD'
        forward_count_field = 'RADF'
        reverse_count_field = 'RADR'
        frequency_field = 'RAF'

    for entry in vcf_reader:
        chromosome = entry.CHROM
        start      = entry.affected_start
        stop       = entry.affected_end
        reference  = entry.REF
        alts       = entry.ALT

        #If we limit the annotations to only SNVs and the entry contains an InDel, skip it
        if args.variant_type == 'snv' and has_indel(entry):
            if has_snv(entry):
                logging.warning("Running in `snv` variant type mode but VCF entry for chr {} pos {} ref {} alts {} contains both SNVs and InDels. Skipping.".format(chromosome, entry.POS, reference, alts))
            write_record(entry, vcf_writer)
            continue

        #If we limit the annotations to only InDels and the entry contains a SNV, skip it
        if args.variant_type == 'indel' and has_snv(entry):
            if has_indel(entry):
                logging.warning("Running in `indel` variant type mode but VCF entry for chr {} pos {} ref {} alts {} contains both SNVs and InDels. Skipping.".format(chromosome, entry.POS, reference, alts))
            write_record(entry, vcf_writer)
            continue

        #If the entry contains a complex variant, skip it
        if has_complex_variant(entry):
            write_record(entry, vcf_writer)
            continue

        (bam_readcount_position, ref_base, var_base) = parse_to_bam_readcount(start, reference, alts[0].serialize(), entry.POS)
        brct = read_counts.get((chromosome,bam_readcount_position,ref_base), None)
        if brct is None:
            add_format_value(entry, sample_name, depth_field, 0)
            if frequency_field not in entry.FORMAT:
                entry.FORMAT += [frequency_field]
            vafs = [0] * len(alts)
            entry.call_for_sample[sample_name].data[frequency_field] = vafs
            if count_field not in entry.FORMAT:
                entry.FORMAT += [count_field]
            ads = [0] * (len(alts) + 1)
            entry.call_for_sample[sample_name].data[count_field] = ads
            write_record(entry, vcf_writer)
            continue

        #Discrepant bam-readcount entries; none of the fields should be written.
        if isinstance(brct, list):
            write_record(entry, vcf_writer)
            continue

        if 'depth' not in brct:
            raise Exception("Error: Malformed bam-readcount output for chromosome {} reference {} position {}. Missing depth field: {}".format(chromosome, reference, start, brct))

        #DP - read depth
        depth = brct['depth']
        add_format_value(entry, sample_name, depth_field, depth)

        #If `depth` is the only key in this hash, then this must have
        #been a duplicate bam-readcount entry where only the depths matched.
        #The only field to write is depth; frequency and count fields should not be written.
        if len(brct.keys()) == 1 and list(brct.keys())[0] == 'depth':
            write_record(entry, vcf_writer)
            continue

        primary_brct = brct
        primary_ref_base = ref_base
        primary_var_base = var_base

        #AF - variant allele frequencies
        if frequency_field not in entry.FORMAT:
            entry.FORMAT += [frequency_field]
        vafs = []
        for alt in alts:
            alt = alt.serialize()
            (bam_readcount_position, ref_base, var_base) = parse_to_bam_readcount(start, reference, alt, entry.POS)
            brct = read_counts.get((chromosome,bam_readcount_position,ref_base), None)
            if brct is not None:
                if var_base not in brct['counts']:
                    print("Warning: variant base {} is not present in the bam-readcount entry for variant {} {}. This might indicate that the bam-readcount file doesn't match the VCF.".format(var_base, chromosome, start))
                    vafs.append(0)
                else:
                    vafs.append(calculate_vaf(int(brct['counts'][var_base]), depth))
            else:
                vafs.append(0)
        entry.call_for_sample[sample_name].data[frequency_field] = vafs

        #AD/ADF/ADR - ref, var1..varN counts
        for (field_name, value_name) in zip([count_field, forward_count_field, reverse_count_field], ['counts', 'forward_counts', 'reverse_counts']):
            if field_name not in entry.FORMAT:
                entry.FORMAT += [field_name]
            (bam_readcount_position, ref_base, var_base) = parse_to_bam_readcount(start, reference, alts[0].serialize(), entry.POS)
            brct = read_counts.get((chromosome,bam_readcount_position,ref_base), None)
            ads = []
            ads.append(brct[value_name][ref_base])
            for alt in alts:
                alt = alt.serialize()
                (bam_readcount_position, ref_base, var_base) = parse_to_bam_readcount(start, reference, alt, entry.POS)
                brct = read_counts.get((chromosome,bam_readcount_position,ref_base), None)
                if brct is not None:
                    if var_base not in brct[value_name]:
                        print("Warning: variant base {} is not present in the bam-readcount entry for variant {} {}. This might indicate that the bam-readcount file doesn't match the VCF.".format(var_base, chromosome, start))
                        ads.append(0)
                    else:
                        ads.append(brct[value_name][var_base])
                else:
                    ads.append(0)
            entry.call_for_sample[sample_name].data[field_name] = ads

        # Skip the whole block if no extra fields were requested, or if the bam-readcount entry doesn't have per-base quality data
        if extra_fields and 'qualities' in primary_brct:
            for f in extra_fields:
                # strand_counts maps to two VCF tags (ADF/ADR) instead of one, so it can't follow the same path
                # as the other extra fields. It's also skipped in DNA mode because ADF/ADR are already written
                # earlier in the function as standard fields. In RNA mode they aren't written by default, so 
                # this is where they get added.
                if f.tag == 'strand_counts':
                    if args.data_type != 'DNA':
                        ref_q = primary_brct['qualities'].get(primary_ref_base, {})
                        var_q = primary_brct['qualities'].get(primary_var_base, {})
                        add_format_value(entry, sample_name, 'ADF', [
                            int(float(ref_q.get('num_plus_strand', '0'))),
                            int(float(var_q.get('num_plus_strand', '0'))),
                        ])
                        add_format_value(entry, sample_name, 'ADR', [
                            int(float(ref_q.get('num_minus_strand', '0'))),
                            int(float(var_q.get('num_minus_strand', '0'))),
                        ])
                else:
                    var_q = primary_brct['qualities'].get(primary_var_base, {})
                    raw = var_q.get(f.brct_col, None)
                    if raw is not None:
                        val = int(float(raw)) if f.vcf_type == 'Integer' else float(raw)
                    else:
                        val = None
                    add_format_value(entry, sample_name, f.tag, val)

        write_record(entry, vcf_writer)

    vcf_writer.close()
    vcf_reader.close()

if __name__ == '__main__':
    main()
