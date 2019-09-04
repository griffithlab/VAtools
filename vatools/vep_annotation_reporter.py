#!/usr/bin/env python

import sys
import os
import argparse
import re
import vcfpy
import tempfile
import csv
import binascii

def define_parser():
    parser = argparse.ArgumentParser('vep-annotation-reporter')
    parser.add_argument(
        "input_vcf",
        help="The VCF file with VEP annotations to report."
    )
    parser.add_argument(
        "vep_fields",
        help="The VEP fields to report. Takes a space-separated list of fields. "
            +"Example:  Consequence SYMBOL Feature",
        nargs='+',
    )
    parser.add_argument(
        "-t", "--input_tsv",
        help="A TSV report file to add VEP annotations to. Required columns are CHROM, POS, REF, ALT. "
            +"These are used to match each TSV entry to a VCF entry. Must be tab-delimited."
    )
    parser.add_argument(
        "-o", "--output-tsv",
        help="Path to write the output report TSV file. If not provided, the output TSV will be written "
            +"next to the input VCF with a .tsv file ending."
    )
    return parser

def create_vcf_reader(args):
    vcf_reader = vcfpy.Reader.from_path(args.input_vcf)
    if 'CSQ' not in vcf_reader.header.info_ids():
        vcf_reader.close()
        raise Exception("ERROR: VCF {} is not VEP-annotated. Please annotate the VCF with VEP before running this tool.".format(args.input_vcf))
    return vcf_reader

def create_tsv_reader(input_filehandle):
    tsv_reader = csv.DictReader(input_filehandle, delimiter = "\t")
    for field in ['CHROM', 'POS', 'REF', 'ALT']:
        if field not in tsv_reader.fieldnames:
            raise Exception("ERROR: Input TSV {} doesn't contain required column '{}'.".format(input_filehandle.name, field))
    return tsv_reader

def parse_csq_header(vcf_reader):
    format_pattern = re.compile('Format: (.*)')
    return format_pattern.search(vcf_reader.header.get_info_field_info('CSQ').description).group(1).split('|')

def parse_csq_entries(csq_entries, csq_fields):
    transcripts = {}
    for entry in csq_entries:
        values = entry.split('|')
        transcript = {}
        for key, value in zip(csq_fields, values):
            transcript[key] = value
        if transcript['Allele'] not in transcripts.keys():
            transcripts[transcript['Allele']] = []
        transcripts[transcript['Allele']].append(transcript)
    return transcripts

def is_sv(entry):
    return 'SVTYPE' in entry.INFO

def resolve_alleles(entry, csq_alleles):
    alleles = {}
    if is_sv(entry):
        for alt in entry.ALT:
            alt = alt.serialize()
            if len(alt) > len(entry.REF) and 'insertion' in csq_alleles:
                alleles[alt] = 'insertion'
            elif len(alt) < len(entry.REF) and 'deletion' in csq_alleles:
                alleles[alt] = 'deletion'
            elif len(csq_alleles) == 1:
                alleles[alt] = list(csq_alleles)[0]
            else:
                print("Warning: Unsupported alleles {} in VEP annotation for SV {}".format(csq_alleles, entry))
    else:
        for alt in entry.ALT:
            alt = alt.serialize()
            if len(alt) != len(entry.REF):
                if alt[0:1] != entry.REF[0:1]:
                    csq_allele = alt
                elif alt[1:] == "":
                    csq_allele = '-'
                else:
                    csq_allele = alt[1:]
                alleles[alt] = csq_allele
            else:
                alleles[alt] = alt
    return alleles

def transcript_for_alt(transcripts, alt):
    for transcript in transcripts[alt]:
        if 'PICK' in transcript and transcript['PICK'] == '1':
            return transcript
    merged_transcripts = {}
    for key in transcripts[alt][0].keys():
        merged_transcripts[key] = ",".join([transcript[key] for transcript in transcripts[alt]])
    return merged_transcripts

def decode_hex(match_string):
    hex_string = match_string.group(0).replace('%', '')
    return binascii.unhexlify(hex_string).decode('utf-8')

def extract_vep_fields(args):
    vcf_reader = create_vcf_reader(args)
    csq_fields = parse_csq_header(vcf_reader)
    vep = {}
    for variant in vcf_reader:
        chr = str(variant.CHROM)
        pos = str(variant.POS)
        ref = str(variant.REF)
        alts = variant.ALT

        if chr not in vep:
            vep[chr] = {}

        if pos not in vep[chr]:
            vep[chr][pos] = {}

        if ref not in vep[chr][pos]:
            vep[chr][pos][ref] = {}

        if 'CSQ' not in variant.INFO:
            for alt in alts:
                vep[chr][pos][ref][alt.serialize()] = None
            continue
        else:
            transcripts = parse_csq_entries(variant.INFO['CSQ'], csq_fields)

        alleles_dict = resolve_alleles(variant, transcripts.keys())
        for alt in alts:
            alt = alt.serialize()
            if alt not in vep[chr][pos][ref]:
                if alleles_dict[alt] in transcripts:
                    vep[chr][pos][ref][alt] = transcript_for_alt(transcripts, alleles_dict[alt])
                else:
                    vep[chr][pos][ref][alt] = None
            else:
                sys.exit("VEP entry for at CHR %s, POS %s, REF %s , ALT % already exists" % (chr, pos, ref, alt) )
    vcf_reader.close()
    return vep

def add_vep_fields_to_row(args, row, vep):
    for field in args.vep_fields:
        field_annotations = []
        for alt in row['ALT'].split(','):
            vep_annotations = vep[row['CHROM']][row['POS']][row['REF']][alt]
            if vep_annotations is not None and field in vep_annotations:
                annotation = vep_annotations[field]
                decoded_annotation = re.sub(r'%[0-9|A-F][0-9|A-F]', decode_hex, annotation)
                field_annotations.append(decoded_annotation)
            else:
                field_annotations.append('-')
        row[field] = ','.join(field_annotations)
    return row

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    vep = extract_vep_fields(args)

    if args.output_tsv:
        output_file = args.output_tsv
    else:
        (head, sep, tail) = args.input_vcf.rpartition('.vcf')
        output_file = "{}.tsv".format(head)

    if args.input_tsv:
        with open(args.input_tsv, 'r') as input_filehandle:
            tsv_reader = create_tsv_reader(input_filehandle)
            output_filehandle = open(output_file, 'w')
            writer = csv.DictWriter(output_filehandle, fieldnames = tsv_reader.fieldnames + args.vep_fields, delimiter = "\t")
            writer.writeheader()
            for entry in tsv_reader:
                row = entry
                row = add_vep_fields_to_row(args, row, vep)
                writer.writerow(row)
            output_filehandle.close()
    else:
        vcf_reader = create_vcf_reader(args)
        with open(output_file, 'w') as output_filehandle:
            writer = csv.DictWriter(output_filehandle, fieldnames = ['CHROM', 'POS', 'REF', 'ALT'] + args.vep_fields, delimiter = "\t")
            writer.writeheader()
            for variant in vcf_reader:
                row = {
                    'CHROM': str(variant.CHROM),
                    'POS'  : str(variant.POS),
                    'REF'  : variant.REF,
                    'ALT'  : ','.join(map(lambda a: a.serialize(), variant.ALT)),
                }
                row = add_vep_fields_to_row(args, row, vep)
                writer.writerow(row)

if __name__ == '__main__':
    main()
