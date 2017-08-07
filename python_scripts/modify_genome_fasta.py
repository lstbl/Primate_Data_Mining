# this will take a fasta file and replace the canonical chromosome names and remove all trailing information past the first space.
# the purpose of this is to make all chromosome names consistent and as short as possible. The "chrXX" identifiers
# are used by UCSC, and I am using the gtf files downloaded from their table browser. I have put the "mappings" from
# their nomenclature to NCBI's nomenclature in the file "<organism_name>_chr.names" in each GENOME directory

import sys, os, gzip

def helpmessage():
    sys.stderr.write('python modify_genome_fasta.py <genome fasta file>\n\n')

def change_fasta(filename):
    matched = False
    if filename[-3:] == '.gz':
        openfn = gzip.open
        new_name = filename[:-3]
    else:
        openfn = open
        new_name = filename
    with openfn(filename) as fh, open(new_name+'.tmp', 'w') as wfile:
        for line in fh:
            if line.startswith('>'):
                line = line.strip().split()
                wfile.write(line[0]+'\n')
            else:
                wfile.write(line)
        wfile.close()
        fh.close()
    os.system('rm {0}; mv {1}.tmp {1}'.format(filename,new_name))

if __name__ == '__main__':
    try:
        filename = sys.argv[1:][0]
        change_fasta(filename)
    except:
        helpmessage()
