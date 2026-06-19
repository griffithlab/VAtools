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
