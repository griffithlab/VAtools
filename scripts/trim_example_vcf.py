#!/usr/bin/env python3
"""
Trim a VCF down to the header lines its data rows actually use.

Used to build the small, readable example VCFs shipped in
vatools-examples.tar.gz (see scripts/build_example_data.sh) from the much
noisier multi-caller fixtures in tests/test_data/. Keeps ##fileformat, a
single PASS ##FILTER line, and only the ##INFO/##FORMAT header lines whose
ID is actually referenced in some record's INFO/FORMAT. Every record's
FILTER is rewritten to PASS. Everything else (contig lines, unused
caller-specific FORMAT/FILTER/INFO declarations, etc.) is dropped.

Usage: trim_example_vcf.py INPUT.vcf > OUTPUT.vcf
"""
import sys
from collections import OrderedDict
import vcfpy


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: trim_example_vcf.py INPUT.vcf > OUTPUT.vcf")

    reader = vcfpy.Reader.from_path(sys.argv[1])
    records = list(reader)
    reader.close()

    info_ids = {key for record in records for key in record.INFO}
    format_ids = {field_id for record in records for field_id in record.FORMAT}

    header = reader.header.copy()
    fileformat_line = next(l for l in header.lines if l.key == 'fileformat')
    seen_info, seen_format = set(), set()
    kept_lines = [fileformat_line]
    for line in header.lines:
        if line.key == 'INFO' and line.mapping['ID'] in info_ids and line.mapping['ID'] not in seen_info:
            kept_lines.append(line)
            seen_info.add(line.mapping['ID'])
        elif line.key == 'FORMAT' and line.mapping['ID'] in format_ids and line.mapping['ID'] not in seen_format:
            kept_lines.append(line)
            seen_format.add(line.mapping['ID'])
    header.lines = kept_lines
    header._indices = header._build_indices()
    header.add_filter_line(OrderedDict([('ID', 'PASS'), ('Description', 'Passed all filters')]))

    writer = vcfpy.Writer.from_stream(sys.stdout, header)
    for record in records:
        record.FILTER = ['PASS']
        writer.write_record(record)
    writer.close()


if __name__ == "__main__":
    main()
