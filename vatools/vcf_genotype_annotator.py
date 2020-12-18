import argparse
import sys
import vcfpy
import csv
from collections import OrderedDict

def create_vcf_reader(args):
    vcf_reader = vcfpy.Reader.from_path(args.input_vcf)
    if args.sample_name in vcf_reader.header.samples.names:
        if 'GT' in vcf_reader.header.format_ids():
            vcf_reader.close()
            raise Exception("VCF already contains a sample column for sample {} with a GT field.".format(args.sample_name))
    return vcf_reader

def create_vcf_writer(args, vcf_reader):
    if args.output_vcf:
        output_file = args.output_vcf
    else:
        (head, sep, tail) = args.input_vcf.rpartition('.vcf')
        output_file = ('').join([head, '.genotype.vcf', tail])
    sample_info = vcf_reader.header.samples
    if args.sample_name in sample_info.names:
        append_to_existing_sample = True
    else:
        append_to_existing_sample = False
        sample_info.names.append(args.sample_name)
        sample_info.name_to_idx[args.sample_name] = len(sample_info.names)-1
    new_header = vcfpy.Header(samples = sample_info)
    for line in vcf_reader.header.lines:
        if not (line.key == 'FORMAT' and line.id == 'GT'):
            new_header.add_line(line)
    new_header.add_format_line(OrderedDict([('ID', 'GT'), ('Number', '1'), ('Type', 'String'), ('Description', 'Genotype')]))
    return ( vcfpy.Writer.from_path(output_file, new_header), append_to_existing_sample )

def define_parser():
    parser = argparse.ArgumentParser(
        "vcf-genotype-annotator",
        description = "A tool to add a new sample to an existing VCF file or fill in the GT field " +
                      "for an existing sample in a VCF."
    )

    parser.add_argument(
        "input_vcf",
        help="A VCF file"
    )
    parser.add_argument(
        "sample_name",
        help="The name of the sample to add",
    )
    parser.add_argument(
        "genotype_value",
        choices=['0/1', '1/1', '0/0', '.'],
        default='0/1',
        help="The genotype value to add to the GT field"
    )
    parser.add_argument(
        "-o", "--output-vcf",
        help="Path to write the output VCF file. If not provided, the output VCF file will be "
            +"written next to the input VCF file with a .genotype.vcf file ending."
    )
    return parser

def main(args_input = sys.argv[1:]):
    parser = define_parser()
    args = parser.parse_args(args_input)

    vcf_reader = create_vcf_reader(args)
    (vcf_writer, append_to_existing_sample) = create_vcf_writer(args, vcf_reader)

    for entry in vcf_reader:
        if "GT" not in entry.FORMAT:
            if isinstance(entry.FORMAT, tuple):
                entry.FORMAT = ["GT"]
            else:
                entry.FORMAT.insert(0, 'GT')
        if append_to_existing_sample:
            entry.call_for_sample[args.sample_name].data['GT'] = args.genotype_value
        else:
            new_sample_call = vcfpy.Call(args.sample_name, data={'GT': args.genotype_value})
            if entry.calls:
                entry.calls.append(new_sample_call)
            else:
                entry.calls = [new_sample_call]
            entry.call_for_sample = {call.sample: call for call in entry.calls}
        vcf_writer.write_record(entry)

    vcf_reader.close()
    vcf_writer.close()

if __name__ == '__main__':
    main()
