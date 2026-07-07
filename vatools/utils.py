import gzip

# handle opening files that may or may not be gzipped
# check the magic bytes at the beginning of the file to determine gzip status
# which is more reliable than looking at the file extension
def open_maybe_gz(path, mode='r'):
    with open(path, 'rb') as f:
        magic = f.read(2)
    if magic == b'\x1f\x8b':
        return gzip.open(path, mode + 't' if 'b' not in mode else mode)
    return open(path, mode)

# vcfpy>=0.14 requires every FORMAT field with Number != 1 to have a list
# value (even when data is missing) for every sample in the record. Some tools
# only populate a field for one sample, so these need to backfill an
# empty list for other samples, which will ultimately be inserted in the 
# VCF as '.'. If a sample has a pre-existing scalar value for a field that 
# used to be declared Number=1 in the input header (e.g. AF) but is now 
# being redeclared as multi-value here, we should keep that value, not discard it.
def fill_missing_multi_value_format_fields(entry, header):
    for field in entry.FORMAT:
        field_info = header.get_format_field_info(field)
        if field_info is None or field_info.number == 1:
            continue
        for call in entry.calls:
            if field not in call.data:
                call.data[field] = []
            elif not isinstance(call.data[field], list):
                call.data[field] = [call.data[field]]

# make sure that missing format fields (in all samples) are filled correctly 
# before writing the VCF record.
def write_record(entry, vcf_writer):
    fill_missing_multi_value_format_fields(entry, vcf_writer.header)
    vcf_writer.write_record(entry)
