import argparse
import sys
import vcfpy
import re
import csv
from collections import OrderedDict
import logging

def resolve_consequence(consequence_string, allele):
    if '&' in consequence_string:
        consequences = {consequence.lower() for consequence in consequence_string.split('&')}
    elif '.' in consequence_string:
        consequences = {consequence.lower() for consequence in consequence_string.split('.')}
    else:
        consequences = [consequence_string.lower()]

    if 'start_lost' in consequences:
        consequence = None
    elif 'stop_retained_variant' in consequences:
        consequence = None
    elif 'frameshift_variant' in consequences:
        consequence = 'FS'
    elif 'missense_variant' in consequences:
        consequence = 'missense'
    elif 'inframe_insertion' in consequences:
        consequence = 'inframe_ins'
    elif 'inframe_deletion' in consequences:
        consequence = 'inframe_del'
    elif 'protein_altering_variant' in consequences:
        if '-' in allele:
            allele = allele.replace('-', '')
            if len(allele) % 3 == 0:
                consequence = 'inframe_del'
            else:
                consequence = None
        elif len(allele) % 3 == 0:
            consequence = 'inframe_ins'
        else:
            consequence = None
    else:
        consequence = None
    return consequence

def create_vcf_reader(args):
    vcf_reader = vcfpy.Reader.from_path(args.input_vcf)
    if 'CSQ' not in vcf_reader.header.info_ids():
        vcf_reader.close()
        raise Exception("ERROR: VCF {} is not VEP-annotated. Please annotate the VCF with VEP before running this tool.".format(args.input_vcf))
    return vcf_reader

def create_vcf_writer(args, vcf_reader):
    (head, sep, tail) = args.input_vcf.rpartition('.vcf')
    new_header = vcf_reader.header.copy()
    if args.filter == 'soft':
        new_header.add_filter_line(OrderedDict([('ID', 'CSQ_MISMATCH'), ('Description', 'Variant has a mismatch between the reference bases (REF) and the reference transcript bases (CSQ) at the mutation position')]))
    if args.output_vcf:
        output_file = args.output_vcf
    else:
        output_file = ('').join([head, '.filtered.vcf', tail])
    return vcfpy.Writer.from_path(output_file, new_header)

def report_fieldnames():
    return [
        'CHROM',
        'POS',
        'Transcript',
        'Mutation Position',
        'Reference Genome Wildtype Amino Acid(s)',
        'Reference Transcript Wildtype Amino Acid(s)',
    ]

def define_parser():
    parser = argparse.ArgumentParser("ref-transcript-mismatch-reporter", description=
        "A tool to identify variants in a VCF where the reference genome used to " +
        "align and call variants doesn't match the Ensembl reference transcript " +
        "used by VEP for variant consequence annotations."
    )

    parser.add_argument(
        "input_vcf",
        help="A VEP-annotated VCF file with Wildtype plugin annotation"
    )
    parser.add_argument(
        "-f", "--filter",
        help='soft: Write a soft-filtered VCF file which identifies variants with mismatched VEP annotations with the CSQ_MISMATCH filter.\n' + 
             'hard: Write a hard-filtered VCF file that removes variants with mismatched VEP annotations.',
        choices=['soft', 'hard']
    )
    parser.add_argument(
        "-o", "--output-vcf",
        help="Path to write the output VCF file to if a --filter is chosen. If not provided, the output VCF file will be "
            +"written next to the input VCF file with a .filtered.vcf file ending."
    )

    return parser

def main(args_input = sys.argv[1:]):
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    parser = define_parser()
    args = parser.parse_args(args_input)

    vcf_reader = create_vcf_reader(args)
    format_pattern = re.compile('Format: (.*)')
    csq_format = format_pattern.search(vcf_reader.header.get_info_field_info('CSQ').description).group(1).split('|')
    if 'WildtypeProtein' not in csq_format:
        raise Exception("VEP was not run with the Wildtype plugin.")

    if args.filter is not None:
        vcf_writer = create_vcf_writer(args, vcf_reader)

    (head, sep, tail) = args.input_vcf.rpartition('.vcf')
    report = ('').join([head, '.mismatches.tsv', tail])
    with open(report, 'w') as output_filehandle:
        writer = csv.DictWriter(output_filehandle, fieldnames = report_fieldnames(), delimiter = "\t")
        writer.writeheader()

        variant_count = 0
        processable_variant_count = 0
        filtered_variant_count = 0
        transcript_count = 0
        processable_transcript_count = 0
        filtered_transcript_count = 0
        for entry in vcf_reader:
            variant_count += 1
            filtered = False
            processable_variant = False
            if 'CSQ' not in entry.INFO:
                logging.warning("Variant is missing VEP annotation. INFO column doesn't contain CSQ field for variant {}".format(entry))
                if args.filter is not None:
                    vcf_writer.write_record(entry)
                continue
            for transcript in entry.INFO['CSQ']:
                transcript_count += 1
                for key, value in zip(csq_format, transcript.split('|')):
                    if key == 'Allele':
                        allele = value
                    if key == 'WildtypeProtein':
                        full_wildtype_sequence = value
                    if key == 'Amino_acids':
                        wildtype_amino_acid = value.split('/')[0]
                        if '*' in wildtype_amino_acid:
                            wildtype_amino_acid = wildtype_amino_acid.split('*')[0]
                        if 'X' in wildtype_amino_acid:
                            wildtype_amino_acid = wildtype_amino_acid.split('X')[0]
                    if key == 'Protein_position':
                        protein_position = value
                        if '/' in value:
                            protein_position = value.split('/')[0]
                            if protein_position == '-':
                                protein_position = value.split('/')[1]
                    if key == 'Consequence':
                        variant_type = resolve_consequence(value, allele)
                    if key == 'Feature':
                        transcript_id = value

                if protein_position == '':
                    continue

                if '*' in full_wildtype_sequence:
                    continue

                if variant_type == 'missense' or variant_type == 'inframe_ins':
                    if '-' in protein_position:
                        position = int(protein_position.split('-', 1)[0]) - 1
                    else:
                        position = int(protein_position) - 1
                elif variant_type == 'inframe_del':
                    position = int(protein_position.split('-', 1)[0]) - 1
                elif variant_type == 'FS':
                    position = int(protein_position.split('-', 1)[0]) - 1
                else:
                    continue

                if position == '-':
                    continue

                if wildtype_amino_acid != '-':
                    processable_transcript_count += 1
                    processable_variant = True
                    mutation_end_position = position + len(wildtype_amino_acid)
                    ref_transcript_wildtype_amino_acid = full_wildtype_sequence[position:mutation_end_position]
                    if wildtype_amino_acid != ref_transcript_wildtype_amino_acid:
                        report_entry = {
                            'CHROM': entry.CHROM,
                            'POS': entry.POS,
                            'Transcript': transcript_id,
                            'Mutation Position': protein_position,
                            'Reference Genome Wildtype Amino Acid(s)': wildtype_amino_acid,
                            'Reference Transcript Wildtype Amino Acid(s)': ref_transcript_wildtype_amino_acid,
                        }
                        writer.writerow(report_entry)
                        filtered = True
                        filtered_transcript_count += 1

            if filtered == True:
                filtered_variant_count += 1

            if processable_variant == True:
                processable_variant_count += 1

            if args.filter == 'soft':
                if filtered == True:
                    entry.add_filter('CSQ_MISMATCH')
                vcf_writer.write_record(entry)
            if args.filter == 'hard':
                if filtered == False:
                    vcf_writer.write_record(entry)

    vcf_reader.close()
    if args.filter is not None:
        vcf_writer.close()

    processable_variant_percentage = "N/A" if processable_variant_count == 0 else round(float(filtered_variant_count)/float(processable_variant_count)*100, 2)
    variant_percentage = "N/A" if variant_count == 0 else round(float(filtered_variant_count)/float(variant_count)*100, 2)
    processable_transcript_percentage = "N/A" if processable_transcript_count == 0 else round(float(filtered_transcript_count)/float(processable_transcript_count)*100, 2)
    transcript_percentage = "N/A" if transcript_count == 0 else round(float(filtered_transcript_count)/float(transcript_count)*100, 2)
    logging.info(
        "\nTotal number of variants: {}\n".format(variant_count) +
        "Total number of processable variants (at least one missense, inframe indels, or frameshift transcript): {}\n".format(processable_variant_count) +
        "Total number of variants with mismatched annotations: {}\n".format(filtered_variant_count) +
        "Percentage of processable variants with mismatched annotations: {}%\n".format(processable_variant_percentage) +
        "Percentage of variants with mismatched annotations: {}%\n".format(variant_percentage) +
        "\n" +
        "Total number of transcripts: {}\n".format(transcript_count) +
        "Total number of processable transcripts (missense, inframe indels, frameshifts): {}\n".format(processable_transcript_count) +
        "Total number of transcripts with mismatched annotations: {}\n".format(filtered_transcript_count) +
        "Percentage of processable transcripts with mismatched annotations: {}%\n".format(processable_transcript_percentage) +
        "Percentage of all transcripts with mismatched annotations: {}%\n".format(transcript_percentage)
    )

if __name__ == '__main__':
    main()
