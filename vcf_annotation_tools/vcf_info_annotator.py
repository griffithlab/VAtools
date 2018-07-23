import argparse
import sys
import vcfpy
import csv
from collections import OrderedDict

def to_array(dictionary):
    array = []
    for key, value in dictionary.items():
        array.append("{}|{}".format(key, value))
    return sorted(array)

def parse_tsv_file(args):
    values={}
    with open(args.values_file,'r') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        for row in tsvin:
            values[(row[0] + ":" + row[1])] = row[2]
    return values


def create_vcf_reader(args):
    vcf_reader = vcfpy.Reader.from_path(args.input_vcf)
    if args.info_field in vcf_reader.header.info_ids():
        vcf_reader.close()
        raise Exception("INFO already contains a {} field. Choose a different label".format(args.info_field))
    return vcf_reader

def create_vcf_writer(args, vcf_reader):
    new_header = vcf_reader.header.copy()
    od = OrderedDict([('ID', args.info_field), ('Number', '.'), ('Type', args.value_format), ('Description', args.description)])
    if args.source:
        od['Source'] = args.source
    if args.version:
        od['Version'] = args.version
    new_header.add_info_line(od)
    output_file = args.output_vcf
    return vcfpy.Writer.from_path(output_file, new_header)

def define_parser():
    parser = argparse.ArgumentParser("vcf-info-annotator")

    parser.add_argument(
        "input_vcf",
        help="A VCF file"
    )
    parser.add_argument(
        "values_file",
        help="A TSV file containing three columns: chromosome, position, value"
    )
    parser.add_argument(
        "info_field",
        help="The INFO field name to add"
    )
    parser.add_argument(
        "description",
        help="The description of the INFO field to add to the header"
    )
    parser.add_argument(
        "value_format",
        choices=['Integer', 'Float', 'Flag', 'Character', 'String'],
        help="The format of the values to be placed into the info field. ",
    )
    parser.add_argument(
        "output_vcf",
        help="Path to write the output VCF file"
    )
    parser.add_argument(
        "-s", "--source",
        help="The string to put in the \"source\" section of the INFO header line - optional "
    )
    parser.add_argument(
        "-v", "--version",
        help="The string to put in the \"version\" section of the INFO header line - optional "
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)
    
    vcf_reader  = create_vcf_reader(args)
    vcf_writer = create_vcf_writer(args, vcf_reader)

    values = parse_tsv_file(args)

    for entry in vcf_reader:
        if entry.CHROM + ":" + str(entry.POS) in values:
            entry.INFO[args.info_field] = [values[entry.CHROM + ":" + str(entry.POS)]]
        vcf_writer.write_record(entry)

    vcf_reader.close()
    vcf_writer.close()

if __name__ == '__main__':
    main()
