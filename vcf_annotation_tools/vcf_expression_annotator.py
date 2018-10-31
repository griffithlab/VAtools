import argparse
import sys
import vcfpy
import pandas as pd
import re
from collections import OrderedDict
from gtfparse import read_gtf
import logging

def resolve_id_column(args):
    if args.format == 'cufflinks':
        return 'tracking_id'
    elif args.format == 'kallisto':
        if args.mode == 'gene':
            return 'gene_name'
        elif args.mode == 'transcript':
            return 'target_id'
    elif args.format == 'stringtie':
        if args.mode == 'gene':
            return 'Gene ID'
    elif args.format == 'custom':
        if args.id_column is None:
            raise Exception("ERROR: `--id-column` option is required when using the `custom` format")
        else:
            return args.id_column

def resolve_stringtie_id_column(args, headers):
    if 'reference_id' in headers:
        return 'reference_id'
    else:
        return 'transcript_id'

def resolve_expression_column(args):
    if args.format == 'cufflinks':
        return 'FPKM'
    elif args.format == 'kallisto':
        if args.mode == 'gene':
            return 'abundance'
        elif args.mode == 'transcript':
            return 'tpm'
    elif args.format == 'stringtie':
        return 'TPM'
    elif args.format == 'custom':
        if args.expression_column is None:
            raise Exception("ERROR: `--id-column` option is required when using the `custom` format")
        else:
            return args.expression_column

def to_array(dictionary):
    array = []
    for key, value in dictionary.items():
        array.append("{}|{}".format(key, value))
    return sorted(array)

def parse_expression_file(args, vcf_reader, vcf_writer):
    if args.format == 'stringtie' and args.mode == 'transcript':
        df_all = read_gtf(args.expression_file)
        df = df_all[df_all["feature"] == "transcript"]
        id_column = resolve_stringtie_id_column(args, df.columns.values)
    else:
        id_column = resolve_id_column(args)
        df = pd.read_csv(args.expression_file, sep='\t')
    expression_column = resolve_expression_column(args)
    if expression_column not in df.columns.values:
        vcf_reader.close()
        vcf_writer.close()
        raise Exception("ERROR: expression_column header {} does not exist in expression_file {}".format(expression_column, args.expression_file))
    if id_column not in df.columns.values:
        vcf_reader.close()
        vcf_writer.close()
        raise Exception("ERROR: id_column header {} does not exist in expression_file {}".format(id_column, args.expression_file))
    return df, id_column, expression_column

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

def create_vcf_writer(args, vcf_reader):
    (head, sep, tail) = args.input_vcf.rpartition('.vcf')
    new_header = vcf_reader.header.copy()
    if args.mode == 'gene':
        new_header.add_format_line(OrderedDict([('ID', 'GX'), ('Number', '.'), ('Type', 'String'), ('Description', 'Gene Expressions')]))
        output_file = ('').join([head, '.gx.vcf', tail])
    elif args.mode == 'transcript':
        new_header.add_format_line(OrderedDict([('ID', 'TX'), ('Number', '.'), ('Type', 'String'), ('Description', 'Transcript Expressions')]))
        output_file = ('').join([head, '.tx.vcf', tail])
    if args.output_vcf:
        output_file = args.output_vcf
    return vcfpy.Writer.from_path(output_file, new_header)

def add_expressions(entry, is_multi_sample, sample_name, df, items, tag, id_column, expression_column, ignore_transcript_version, missing_expressions_count, entry_count):
    expressions = {}
    for item in items:
        entry_count += 1
        if tag == 'TX' and ignore_transcript_version:
            df['transcript_without_version'] = df[id_column].apply(lambda x: re.sub(r'\.[0-9]+$', '', x))
            subset = df.loc[df['transcript_without_version'] == re.sub(r'\.[0-9]+$', '', item)]
        else:
            subset = df.loc[df[id_column] == item]
        if len(subset) > 0:
            expressions[item] = subset[expression_column].sum()
        else:
            missing_expressions_count += 1
    if is_multi_sample:
        entry.FORMAT += [tag]
        entry.call_for_sample[sample_name].data[tag] = to_array(expressions)
    else:
        entry.add_format(tag, to_array(expressions))
    return (entry, missing_expressions_count, entry_count)

def define_parser():
    parser = argparse.ArgumentParser("vcf-expression-annotator")

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
    parser.add_argument(
        "-o", "--output-vcf",
        help="Path to write the output VCF file. If not provided, the output VCF file will be "
            +"written next to the input VCF file with a .tx.vcf or .gx.vcf file ending."
    )
    parser.add_argument(
        "--ignore-transcript-version",
        help='Assumes that the final period and number denotes the transcript version and ignores it (i.e. for "ENST00001234.3" - ignores the ".3").',
        action="store_true"
    )

    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    if args.format == 'custom':
        if args.id_column is None:
            raise Exception("--id-column is not set. This is required when using the `custom` format.")
        if args.expression_column is None:
            raise Exception("--expression-column is not set. This is required when using the `custom` format.")

    (vcf_reader, is_multi_sample) = create_vcf_reader(args)
    format_pattern = re.compile('Format: (.*)')
    csq_format = format_pattern.search(vcf_reader.header.get_info_field_info('CSQ').description).group(1).split('|')

    vcf_writer = create_vcf_writer(args, vcf_reader)

    (df, id_column, expression_column) = parse_expression_file(args, vcf_reader, vcf_writer)

    missing_expressions_count = 0
    entry_count = 0
    for entry in vcf_reader:
        transcript_ids = set()
        genes = set()
        if 'CSQ' not in entry.INFO:
            logging.warning("Variant is missing VEP annotation. INFO column doesn't contain CSQ field for variant {}".format(entry))
            vcf_writer.write_record(entry)
            continue
        for transcript in entry.INFO['CSQ']:
            for key, value in zip(csq_format, transcript.split('|')):
                if key == 'Feature' and value != '' and not value.startswith('ENSR'):
                    transcript_ids.add(value)
                if args.format == 'kallisto':
                    if key == 'SYMBOL' and value != '':
                        genes.add(value)
                else:
                    if key == 'Gene' and value != '':
                        genes.add(value)

        if args.mode == 'gene':
            genes = list(genes)
            if len(genes) > 0:
                (entry, missing_expressions_count, entry_count) = add_expressions(entry, is_multi_sample, args.sample_name, df, genes, 'GX', id_column, expression_column, False, missing_expressions_count, entry_count)
        elif args.mode == 'transcript':
            transcript_ids = list(transcript_ids)
            if len(transcript_ids) > 0:
                (entry, missing_expressions_count, entry_count) = add_expressions(entry, is_multi_sample, args.sample_name, df, transcript_ids, 'TX', id_column, expression_column, args.ignore_transcript_version, missing_expressions_count, entry_count)
        vcf_writer.write_record(entry)

    vcf_reader.close()
    vcf_writer.close()

    if missing_expressions_count > 0:
        logging.warning("{} of {} transcripts did not have an expression entry for their {} id.".format(missing_expressions_count, entry_count, args.mode))

if __name__ == '__main__':
    main()
