import argparse
import sys
import os
import re
import vcfpy
import tempfile
import csv
from collections import OrderedDict

def define_parser():
    parser = argparse.ArgumentParser('vcf-readcount-annotator')
    parser.add_argument(
        "input_vcf",
        help="A VCF file"
    )
    parser.add_argument(
        "bam_readcount_file",
        help="A bam-readcount output file",
        nargs='+',
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
    return parser

def parse_brct_field(brcts):
    parsed_brct = {}
    for brct in brcts:
        (base, count, rest) = brct.split(':', 2)
        parsed_brct[base.upper()] = count
    return parsed_brct

def parse_bam_readcount_file(args):
    coverage = {}
    for bam_readcount_file in args.bam_readcount_file:
        with open(bam_readcount_file, 'r') as reader:
            coverage_tsv_reader = csv.reader(reader, delimiter='\t')
            for row in coverage_tsv_reader:
                chromosome     = row[0]
                position       = row[1]
                reference_base = row[2].upper()
                depth          = row[3]
                brct           = row[4:]
                parsed_brct = parse_brct_field(brct)
                parsed_brct['depth'] = depth
                coverage[(chromosome,position,reference_base)] = parsed_brct
    return coverage

def is_insertion(ref, alt):
    return len(alt) > len(ref)

def is_deletion(ref, alt):
    return len(alt) < len(ref)

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
    is_multi_sample = len(vcf_reader.header.samples.names) > 1
    if is_multi_sample:
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

def create_vcf_writer(args, vcf_reader):
    if args.output_vcf:
        output_file = args.output_vcf
    else:
        (head, sep, tail) = args.input_vcf.rpartition('.vcf')
        output_file = ('').join([head, '.readcount.vcf', tail])
    new_header = vcfpy.Header(samples = vcf_reader.header.samples)
    if args.data_type == 'DNA':
        for line in vcf_reader.header.lines:
            if not (line.key == 'FORMAT' and line.id in ['DP', 'AD', 'AF']):
                new_header.add_line(line)
        new_header.add_format_line(OrderedDict([('ID', 'DP'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Read depth')]))
        new_header.add_format_line(OrderedDict([('ID', 'AD'), ('Number', 'R'), ('Type', 'Integer'), ('Description', 'Allelic depths for the ref and alt alleles in the order listed')]))
        new_header.add_format_line(OrderedDict([('ID', 'AF'), ('Number', 'A'), ('Type', 'Float'), ('Description', 'Variant-allele frequency for the alt alleles')]))
    if args.data_type == 'RNA':
        for line in vcf_reader.header.lines:
            if not (line.key == 'FORMAT' and line.id in ['RDP', 'RAD', 'RAF']):
                new_header.add_line(line)
        new_header.add_format_line(OrderedDict([('ID', 'RDP'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'RNA Read depth')]))
        new_header.add_format_line(OrderedDict([('ID', 'RAD'), ('Number', 'R'), ('Type', 'Integer'), ('Description', 'RNA Allelic depths for the ref and alt alleles in the order listed')]))
        new_header.add_format_line(OrderedDict([('ID', 'RAF'), ('Number', 'A'), ('Type', 'Float'), ('Description', 'RNA Variant-allele frequency for the alt alleles')]))
    return vcfpy.Writer.from_path(output_file, new_header)

def write_depth(entry, sample_name, field, value):
    if field not in entry.FORMAT:
        entry.FORMAT += [field]
    entry.call_for_sample[sample_name].data[field] = value

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    read_counts = parse_bam_readcount_file(args)
    (vcf_reader, sample_name) = create_vcf_reader(args)
    vcf_writer = create_vcf_writer(args, vcf_reader)

    if args.data_type == 'DNA':
        depth_field = 'DP'
        count_field = 'AD'
        frequency_field = 'AF'
    elif args.data_type == 'RNA':
        depth_field = 'RDP'
        count_field = 'RAD'
        frequency_field = 'RAF'

    for entry in vcf_reader:
        chromosome = entry.CHROM
        start      = entry.affected_start
        stop       = entry.affected_end
        reference  = entry.REF
        alts       = entry.ALT

        genotype_alts = [a for a in entry.call_for_sample[sample_name].gt_bases if a != reference]
        if len(list(set(genotype_alts))) == 1:
            genotype_alt = genotype_alts[0]
        else:
            vcf_writer.write_record(entry)
            continue

        (bam_readcount_position, ref_base, var_base) = parse_to_bam_readcount(start, reference, genotype_alt, entry.POS)
        brct = read_counts.get((chromosome,bam_readcount_position,ref_base), None)
        if brct is None:
            write_depth(entry, sample_name, depth_field, 0)
            if frequency_field not in entry.FORMAT:
                entry.FORMAT += [frequency_field]
            vafs = [0] * len(alts)
            entry.call_for_sample[sample_name].data[frequency_field] = vafs
            if count_field not in entry.FORMAT:
                entry.FORMAT += [count_field]
            ads = [0] * (len(alts) + 1)
            entry.call_for_sample[sample_name].data[count_field] = ads
            vcf_writer.write_record(entry)
            continue

        #DP - read depth
        depth = brct['depth']
        write_depth(entry, sample_name, depth_field, depth)

        #AF - variant allele frequencies
        if frequency_field not in entry.FORMAT:
            entry.FORMAT += [frequency_field]
        vafs = []
        for alt in alts:
            alt = alt.serialize()
            (bam_readcount_position, ref_base, var_base) = parse_to_bam_readcount(start, reference, alt, entry.POS)
            brct = read_counts.get((chromosome,bam_readcount_position,ref_base), None)
            if brct is not None:
                if var_base not in brct:
                    print("Warning: variant base {} is not present in the bam-readcount entry for variant {} {}. This might indicate that the bam-readcount file doesn't match the VCF.".format(var_base, chromosome, start))
                    vafs.append(0)
                else:
                    vafs.append(calculate_vaf(int(brct[var_base]), depth))
            else:
                vafs.append(0)
        entry.call_for_sample[sample_name].data[frequency_field] = vafs

        #AD - ref, var1..varN counts
        if count_field not in entry.FORMAT:
            entry.FORMAT += [count_field]
        ads = []
        ads.append(brct[ref_base])
        for alt in alts:
            alt = alt.serialize()
            (bam_readcount_position, ref_base, var_base) = parse_to_bam_readcount(start, reference, alt, entry.POS)
            brct = read_counts.get((chromosome,bam_readcount_position,ref_base), None)
            if brct is not None:
                if var_base not in brct:
                    print("Warning: variant base {} is not present in the bam-readcount entry for variant {} {}. This might indicate that the bam-readcount file doesn't match the VCF.".format(var_base, chromosome, start))
                    ads.append(0)
                else:
                    ads.append(brct[var_base])
            else:
                ads.append(0)
        entry.call_for_sample[sample_name].data[count_field] = ads

        vcf_writer.write_record(entry)

    vcf_writer.close()
    vcf_reader.close()

if __name__ == '__main__':
    main()
