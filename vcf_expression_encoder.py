import argparse
import sys

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

if __name__ == '__main__':
    main()
