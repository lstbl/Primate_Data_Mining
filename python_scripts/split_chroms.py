# this will split entries from a fasta file into seperate files (useful for parsing VCF files)
# use as python split_chroms.py <fasta_file>

import gzip
import sys
def split_chroms(fasta_file):
    if fasta_file[-3:] == '.gz':
        openfn = gzip.open
    else:
        openfn = open
    with openfn(fasta_file) as fh:
        for line in fh:
            if line.startswith('>'):
                current_chr = line.strip()[1:]
                wfile = open(current_chr+'.fa','w')
                wfile.write(">"+current_chr+'\n')
            elif line.strip():
                wfile.write(line)
    wfile.close()

if __name__ == '__main__':
    filename = sys.argv[1:][0]
    split_chroms(filename)

