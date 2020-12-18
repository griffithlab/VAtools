#!/usr/bin/env python

import sys
import os
import argparse
import re
import vcfpy
import tempfile
import csv
import binascii
from statistics import mean, median, stdev
import logging

def define_parser():
    parser = argparse.ArgumentParser(
        'transform-split-values',
        description = "A tool that extracts and manipulates values from existing sample fields " +
                      "and outputs the results to a TSV file."
    )
    parser.add_argument(
        "input_vcf",
        help="The VCF file from which to extract information. Multi-allelic sites must be decomposed."
    )
    parser.add_argument(
        "format_field",
        help="The multi-value format field to report.",
    )
    parser.add_argument(
        "operations",
        help="The operation to execute on the chosen field.\n"
            +"ref:  Extract the first value in a R-number field (the reference value).\n"
            +"alt:  Extract the second value in a R-number field (the alt value).\n"
            +"sum:  Calculate the sum of all the numbers in the field.\n"
            +"min:  Calculate the minimum of all the numbers in the field.\n"
            +"max:  Calculate the maximum of all the numbers in the field.\n"
            +"mean: Calculate the mean of all the numbers in the field.\n"
            +"median: Calculate the median of all the numbers in the field.\n"
            +"stdev: Calculate the standard deviation of all the numbers in the field.\n"
            +"ref_ratio:  The first value in a R-number field divided by the sum of all the numbers (the reference ratio).\n"
            +"alt_ratio:  The second value in a R-number field divided by the sum of all the numbers (the alt ratio).\n",
        choices=['ref', 'alt', 'sum', 'min', 'max', 'mean', 'median', 'stdev', 'ref_ratio', 'alt_ratio'],
        nargs='+',
    )
    parser.add_argument(
        "-t", "--input_tsv",
        help="A TSV report file to add information to. Required columns are CHROM, POS, REF, ALT. "
            +"These are used to match each TSV entry to a VCF entry. Must be tab-delimited."
    )
    parser.add_argument(
        "-s", "--sample-name",
        help="If the input_vcf contains multiple samples, the name of the sample to extract information for."
    )
    parser.add_argument(
        "-o", "--output-tsv",
        help="Path to write the output report TSV file. If not provided, the output TSV will be written "
            +"next to the input VCF with a .tsv file ending."
    )
    return parser

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
        if args.sample_name != sample_name:
            logging.warn("VCF doesn't contain requested sample {}. Using sample {}".format(args.sample_name, sample_name))
    if args.format_field not in vcf_reader.header.format_ids():
        vcf_reader.close()
        raise Exception("ERROR: VCF {} does not contain a format field {}.".format(args.input_vcf, args.format_field))
    format_field_info = vcf_reader.header.get_format_field_info(args.format_field)
    for operation in args.operations:
        if operation in ['ref', 'alt', 'ref_ratio', 'alt_ratio'] and format_field_info.number != 'R':
            vcf_reader.close()
            raise Exception("ERROR: VCF {} format field {} incompatible with operation {}. Not of Number R.".format(args.input_vcf, args.format_field, operation))
        if operation in ['sum', 'min', 'max', 'mean', 'median', 'stdev', 'ref_ratio', 'alt_ratio'] and format_field_info.type not in ['Integer', 'Float']:
            vcf_reader.close()
            raise Exception("ERROR: VCF {} format field {} incompatible with operation {}. Not of Type Integer or Float.".format(args.input_vcf, args.format_field, operation))
    return vcf_reader, sample_name

def create_tsv_reader(input_filehandle):
    tsv_reader = csv.DictReader(input_filehandle, delimiter = "\t")
    for field in ['CHROM', 'POS', 'REF', 'ALT']:
        if field not in tsv_reader.fieldnames:
            raise Exception("ERROR: Input TSV {} doesn't contain required column '{}'.".format(input_filehandle.name, field))
    return tsv_reader

def ref(field_values):
    return handle_empty_array(field_values, "field_values[0]")

def alt(field_values):
    return handle_empty_array(field_values, "field_values[1]")

def ref_ratio(field_values):
    return handle_empty_array(field_values, "field_values[0] / sum(field_values)")

def alt_ratio(field_values):
    return handle_empty_array(field_values, "field_values[1] / sum(field_values)")

def handle_empty_array(field_values, action):
    if len(field_values) == 0:
        return '-'
    else:
        return eval(action)

def extract_field_values(args):
    (vcf_reader, sample_name) = create_vcf_reader(args)
    values = {}
    for variant in vcf_reader:
        chr = str(variant.CHROM)
        pos = str(variant.POS)
        reference = str(variant.REF)
        alts = variant.ALT

        if chr not in values:
            values[chr] = {}
        if pos not in values[chr]:
            values[chr][pos] = {}
        if reference not in values[chr][pos]:
            values[chr][pos][reference] = {}
        alts_string = ','.join(a.serialize() for a in alts)
        if alts_string not in values[chr][pos][reference]:
            values[chr][pos][reference][alts_string] = {}

        for operation in args.operations:
            if operation in ['ref', 'alt', 'ref_ratio', 'alt_ratio'] and len(alts) > 1:
                logging.warn("Multi-allelic sites are not supported for operation {}. Skipping entry {} {}.".format(operation, chr, pos))
                values[chr][pos][reference][alts_string][operation] = '-'
                continue

            if args.format_field in variant.call_for_sample[sample_name].data:
                field_value = variant.call_for_sample[sample_name].data[args.format_field]
                result = eval("{}({})".format(operation, field_value))
                if result != '-':
                    result = round(result, 4)
                values[chr][pos][reference][alts_string][operation] = result
            else:
                values[chr][pos][reference][alts_string][operation] = '-'
    vcf_reader.close()
    return (values, sample_name)

def field_names(args, sample_name):
    return list(map(lambda operation: "{}-{}-{}".format(sample_name, args.format_field, operation), args.operations))

def add_field_values_to_row(args, row, field_values, sample_name):
    field_annotations = []
    if row['CHROM'] in field_values and row['POS'] in field_values[row['CHROM']] and row['REF'] in field_values[row['CHROM']][row['POS']] and row['ALT'] in field_values[row['CHROM']][row['POS']][row['REF']]:
        value = field_values[row['CHROM']][row['POS']][row['REF']][row['ALT']]
        for operation in args.operations:
            if value is not None and operation in value:
                row["{}-{}-{}".format(sample_name, args.format_field, operation)] = value[operation]
            else:
                row["{}-{}-{}".format(sample_name, args.format_field, operation)] = '-'
    else:
        for operation in args.operations:
            row["{}-{}-{}".format(sample_name, args.format_field, operation)] = '-'
    return row

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    (field_values, sample_name) = extract_field_values(args)

    if args.output_tsv:
        output_file = args.output_tsv
    else:
        (head, sep, tail) = args.input_vcf.rpartition('.vcf')
        output_file = "{}.tsv".format(head)

    if args.input_tsv:
        with open(args.input_tsv, 'r') as input_filehandle:
            tsv_reader = create_tsv_reader(input_filehandle)
            output_filehandle = open(output_file, 'w')
            writer = csv.DictWriter(output_filehandle, fieldnames = tsv_reader.fieldnames + field_names(args, sample_name), delimiter = "\t")
            writer.writeheader()
            for entry in tsv_reader:
                row = entry
                row = add_field_values_to_row(args, row, field_values, sample_name)
                writer.writerow(row)
            output_filehandle.close()
    else:
        (vcf_reader, sample_name) = create_vcf_reader(args)
        with open(output_file, 'w') as output_filehandle:
            writer = csv.DictWriter(output_filehandle, fieldnames = ['CHROM', 'POS', 'REF', 'ALT'] + field_names(args, sample_name), delimiter = "\t")
            writer.writeheader()
            for variant in vcf_reader:
                row = {
                    'CHROM': str(variant.CHROM),
                    'POS'  : str(variant.POS),
                    'REF'  : variant.REF,
                    'ALT'  : ','.join(a.serialize() for a in variant.ALT)
                }
                row = add_field_values_to_row(args, row, field_values, sample_name)
                writer.writerow(row)

if __name__ == '__main__':
    main()
