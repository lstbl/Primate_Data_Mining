# this function will use a constantly-updating hash function to compare the chromosomes of different hash
# functions. The purpose is to figure out if two genome files are identical, or if they are not, what parts are
# the same and which are different. This will allow for more easy comparisons and no ambiguity when performing liftover
# and mapping
# NOTE: this function will take stdin and output stdout. The chromosomes will be ordered by length.
# To run you'd do something like:
# cat fasta_file.fa | python hash_fasta.py > hashed_chromosomes.txt    or, if the fasta file was gziped:
# gzip -cd fasta_file.fa.gz | python hash_fasta.py > hashed_chromosomes.txt

import sys
import hashlib
from operator import itemgetter
def hash_fasta_entries():
    chr_info = []
    chr_name = ''
    chr_len = 0
    chrhash = hashlib.md5()
    for line in sys.stdin:
        if line.startswith('>'):
            if chr_len > 0:
                chr_info.append((chr_name,chr_len,chrhash.hexdigest()))
            chr_name = line.strip()[1:]
            chrhash = hashlib.md5()
            chr_len = 0
        else:
            line = line.strip()
            chr_len += len(line)
            chrhash.update(line.upper())
    chr_info = sorted(chr_info,key=itemgetter(1))[::-1]
    sys.stdout.write('chr_name\tchr_len\tchr_hash\n')
    for entry in chr_info:
        sys.stdout.write('{0}\t{1}\t{2}\n'.format(*entry))

if __name__ == '__main__':
    hash_fasta_entries()
