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
            if any(x.strip() for x in row): #skip blank lines
                values[(row[0] + ":" + row[1])] = row[2]
    return values


def create_vcf_reader(args):
    vcf_reader = vcfpy.Reader.from_path(args.input_vcf)
    if args.info_field in vcf_reader.header.info_ids() and not args.overwrite:
        vcf_reader.close()
        raise Exception("INFO already contains a {} field. Choose a different label, or use the --overwrite flag to retain this field and overwrite values".format(args.info_field))
    if args.overwrite and not args.info_field in vcf_reader.header.info_ids():
        vcf_reader.close()
        raise Exception("INFO field {} does not exist and thus cannot be overwritten!".format(args.info_field))
    return vcf_reader

def create_vcf_writer(args, vcf_reader):
    if args.output_vcf:
        output_file = args.output_vcf
    else:
        (head, sep, tail) = args.input_vcf.rpartition('.vcf')
        output_file = ('').join([head, '.info.vcf', tail])

    new_header = vcf_reader.header.copy()

    if args.overwrite:
        if args.value_format or args.description:
            print("Warning: --overwrite flag is set, so existing header for {} field will be retained and --value_format and --description inputs will be ignored".format(args.info_field))
    else:
        od = OrderedDict([('ID', args.info_field), ('Number', '1'), ('Type', args.value_format), ('Description', args.description)])
        if args.source:
            od['Source'] = args.source
        if args.version:
            od['Version'] = args.version
        new_header.add_info_line(od)
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
        help="The INFO field name to use"
    )
    parser.add_argument(
        "-d", "--description",
        help="The description of the INFO field to add to the header"
    )
    parser.add_argument(
        "-f", "--value_format",
        choices=['Integer', 'Float', 'Flag', 'Character', 'String'],
        help="The format of the values to be placed into the info field.",
    )
    parser.add_argument(
        "-o", "--output-vcf",
        help="Path to write the output VCF file. If not provided, the output VCF file will be "
            +"written next to the input VCF file with a .info.vcf file ending."
    )
    parser.add_argument(
        "-s", "--source",
        help="The string to put in the \"source\" section of the INFO header line - optional"
    )
    parser.add_argument('-w', "--overwrite", action='store_true',
        help="by default, ths tool will raise an exception if the field specified already exists in the VCF. This flag allows existing fields to be overwritten."
    )
    parser.add_argument(
        "-v", "--version",
        help="The string to put in the \"version\" section of the INFO header line - optional"
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    #if we're not overwriting an existing field, then these are required
    if not args.overwrite:
        if args.description is None or args.value_format is None:
            raise Exception("the --description and --value_format arguments are required unless updating/overwriting an existing field (with flag --overwrite)")

    vcf_reader  = create_vcf_reader(args)
    vcf_writer = create_vcf_writer(args, vcf_reader)

    values = parse_tsv_file(args)

    for entry in vcf_reader:
        if entry.CHROM + ":" + str(entry.POS) in values:
            entry.INFO[args.info_field] = values[entry.CHROM + ":" + str(entry.POS)]
        vcf_writer.write_record(entry)

    vcf_reader.close()
    vcf_writer.close()

if __name__ == '__main__':
    main()
