#!/usr/bin/env python3
"""
Trim a VCF down to the header lines its data rows actually use.

Used to build the small, readable example VCFs shipped in
vatools-examples.tar.gz (see scripts/build_example_data.sh) from the much
noisier multi-caller fixtures in tests/test_data/. Keeps ##fileformat, a
single PASS ##FILTER line, and only the ##INFO/##FORMAT header lines whose
ID is actually referenced in some data row's INFO/FORMAT column. Every
record's FILTER is rewritten to PASS. Everything else (contig lines,
unused caller-specific FORMAT/FILTER/INFO declarations, etc.) is dropped.

Usage: trim_example_vcf.py INPUT.vcf > OUTPUT.vcf
"""
import re
import sys

ID_RE = re.compile(r"ID=([^,>]+)")


def used_ids(data_lines, column_index):
    ids = set()
    for line in data_lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) <= column_index:
            continue
        if column_index == 7:  # INFO
            for entry in fields[7].split(";"):
                ids.add(entry.split("=", 1)[0])
        else:  # FORMAT
            ids.update(fields[column_index].split(":"))
    return ids


def header_id(meta_line):
    match = ID_RE.search(meta_line)
    return match.group(1) if match else None


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: trim_example_vcf.py INPUT.vcf > OUTPUT.vcf")

    with open(sys.argv[1]) as f:
        lines = f.readlines()

    meta_lines = [l for l in lines if l.startswith("##")]
    chrom_line = next(l for l in lines if l.startswith("#CHROM"))
    data_lines = [l for l in lines if l and not l.startswith("#")]

    info_ids = used_ids(data_lines, 7)
    format_ids = used_ids(data_lines, 8)

    fileformat_line = next(l for l in meta_lines if l.startswith("##fileformat="))

    out_meta = [fileformat_line, '##FILTER=<ID=PASS,Description="Passed all filters">\n']
    seen_info, seen_format = set(), set()
    for line in meta_lines:
        line_id = header_id(line)
        if line.startswith("##INFO=") and line_id in info_ids and line_id not in seen_info:
            out_meta.append(line)
            seen_info.add(line_id)
        elif line.startswith("##FORMAT=") and line_id in format_ids and line_id not in seen_format:
            out_meta.append(line)
            seen_format.add(line_id)

    sys.stdout.writelines(out_meta)
    sys.stdout.write(chrom_line)
    for line in data_lines:
        fields = line.rstrip("\n").split("\t")
        fields[6] = "PASS"
        sys.stdout.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    main()
