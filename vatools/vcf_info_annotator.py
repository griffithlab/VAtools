import argparse
import sys
import vcfpy
import csv
from collections import OrderedDict
from vatools.utils import open_maybe_gz

def to_array(dictionary):
    array = []
    for key, value in dictionary.items():
        array.append("{}|{}".format(key, value))
    return sorted(array)

def parse_tsv_file(args):
    values={}
    with open_maybe_gz(args.values_file) as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        header = next(tsvin)
        for row in tsvin:
            if any(x.strip() for x in row):
                key = row[0] + ":" + row[1]
                values[key] = {col: row[header.index(col)] for col in [m['source_col'] for m in mappings]}
    return values

def create_vcf_reader(args, mappings):
    vcf_reader = vcfpy.Reader.from_path(args.input_vcf)
    for m in mappings:
        if m['info_field'] in vcf_reader.header.info_ids() and not args.overwrite:
            vcf_reader.close()
            raise Exception("INFO already contains a {} field. Choose a different label, or use the --overwrite flag to retain this field and overwrite values".format(m['info_field']))
    return vcf_reader

def create_vcf_writer(args, vcf_reader, mappings):
    if args.output_vcf:
        output_file = args.output_vcf
    else:
        (head, sep, tail) = args.input_vcf.rpartition('.vcf')
        output_file = ('').join([head, '.info.vcf', tail])

    new_header = vcf_reader.header.copy()

    for m in mappings:
        if args.overwrite and m['info_field'] in new_header.info_ids():
            print("WARNING: overwriting existing INFO header line for '{}'".format(m['info_field']), file=sys.stderr)
            new_header.lines = [
                l for l in new_header.lines
                if not (hasattr(l, 'key') and l.key == 'INFO' and l.mapping.get('ID') == m['info_field'])
            ]
            new_header._indices.pop(('INFO', m['info_field']), None)
        od = OrderedDict([('ID', m['info_field']), ('Number', '1'), ('Type', m['type']), ('Description', m['description'])])
        if m['source']:
            od['Source'] = m['source']
        if m['version']:
            od['Version'] = m['version']
        new_header.add_info_line(od)
    return vcfpy.Writer.from_path(output_file, new_header)

def define_parser():
    parser = argparse.ArgumentParser(
        "vcf-info-annotator",
        description = "A tool that will add data from a tab-delimited file to a user-specified " +
                      "field in the VCF INFO column."
    )
    parser.add_argument(
        "input_vcf",
        help="A VCF file"
    )
    parser.add_argument(
        "values_file",
        help="A TSV file with a header row. The first two columns must be chromosome and position. "
            +"One or more additional columns can be mapped to VCF INFO fields via --column-mappings."
    )
    parser.add_argument(
        "-m", "--column-mappings",
        required=True,
        help="Comma-separated list of mappings from TSV columns to VCF INFO fields. "
            +"Required fields: source_col:info_field:type:description. "
            +"Optional fields: :source and :version appended after description. "
            +"Example: \"mq:MQ0:Integer:Mapping quality:dbSNP:1.0,class:CVCLASS:String:ClinVar Classification\""
    )
    parser.add_argument(
        "-o", "--output-vcf",
        help="Path to write the output VCF file. If not provided, the output VCF file will be "
            +"written next to the input VCF file with a .info.vcf file ending."
    )
    parser.add_argument('-w', "--overwrite", action='store_true',
        help="By default, this tool will raise an exception if a mapped field already exists in the VCF. "
            +"This flag allows existing fields to be overwritten."
    )
    parser.add_argument('--clear-existing', action='store_true',
        help="When overwriting, remove preexisting values for mapped INFO fields from every record, "
            +"not just those with a TSV match. Requires --overwrite."
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    mappings = parse_column_mappings(args.column_mappings)
    if args.clear_existing and not args.overwrite:
        parser.error("--clear-existing requires --overwrite")
    vcf_reader = create_vcf_reader(args, mappings)
    vcf_writer = create_vcf_writer(args, vcf_reader, mappings)
    values = parse_tsv_file(args, mappings)

    for entry in vcf_reader:
        key = entry.CHROM + ":" + str(entry.POS)
        if args.clear_existing:
            for m in mappings:
                entry.INFO.pop(m['info_field'], None)
        if key in values:
            row_vals = values[key]
            for m in mappings:
                if m['source_col'] in row_vals:
                    entry.INFO[m['info_field']] = coerce_value(row_vals[m['source_col']], m['type'])
        vcf_writer.write_record(entry)

    vcf_reader.close()
    vcf_writer.close()

if __name__ == '__main__':
    main()
