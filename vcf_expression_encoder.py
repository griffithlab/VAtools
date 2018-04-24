import argparse
import sys
import vcfpy
import re

def create_vcf_reader(args):
    vcf_reader = vcfpy.Reader.from_path(args.input_vcf)
    is_multi_sample = len(vcf_reader.header.samples.names) > 1
    if is_multi_sample and args.sample_name is None:
        vcf_reader.close()
        raise Exception("ERROR: VCF {} contains more than one sample. Please use the -s option to specify which sample to annotate.".format(args.input_vcf))
    elif is_multi_sample and args.sample_name not in vcf_reader.header.samples.names:
        vcf_reader.close()
        raise Exception("ERROR: VCF {} does not contain a sample column for sample {}.".format(args.input_vcf, args.sample_name))
    if 'CSQ' not in vcf_reader.header.info_ids():
        vcf_reader.close()
        raise Exception("ERROR: VCF {} is not VEP-annotated. Please annotate the VCF with VEP before running this tool.".format(args.input_vcf))
    if args.mode == 'gene' and 'GX' in vcf_reader.header.format_ids():
        vcf_reader.close()
        raise Exception("ERROR: VCF {} is already gene expression annotated. GX format header already exists.".format(args.input_vcf))
    elif args.mode == 'transcript' and 'TX' in vcf_reader.header.format_ids():
        vcf_reader.close()
        raise Exception("ERROR: VCF {} is already transcript expression annotated. TX format header already exists.".format(args.input_vcf))
    return vcf_reader, is_multi_sample

def define_parser():
    parser = argparse.ArgumentParser("vcf-expression-encoder")

    parser.add_argument(
        "input_vcf",
        help="A VEP-annotated VCF file"
    )
    parser.add_argument(
        "expression_file",
        help="A TSV file containing expression estimates"
    )
    parser.add_argument(
        "format",
        choices=['kallisto','stringtie','cufflinks','custom'],
        help="The file format of the expression file to process. "
            +"Use `custom` to process file formats not explicitly supported. "
            +"The `custom` option requires the use of the --id-column and --expression-column arguments."
    )
    parser.add_argument(
        "mode",
        choices=['gene', 'transcript'],
        help="The type of expression data in the expression_file"
    )
    parser.add_argument(
        "-i", "--id-column",
        help="The column header in the expression_file for the column containing gene names/transcript ids. Required when using the `custom` format."
    )
    parser.add_argument(
        "-e", "--expression-column",
        help="The column header in the expression_file for the column containing expression data. Required when using the `custom` format."
    )
    parser.add_argument(
        "-s", "--sample-name",
        help="If the input_vcf contains multiple samples, the name of the sample to annotate."
    )

    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    (vcf_reader, is_multi_sample) = create_vcf_reader(args)
    format_pattern = re.compile('Format: (.*)')
    csq_format = format_pattern.search(vcf_reader.header.get_info_field_info('CSQ').description).group(1).split('|')

if __name__ == '__main__':
    main()
