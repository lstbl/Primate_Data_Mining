# this will take a fasta file and replace the canonical chromosome names with just the number indicators to be consistent
# with the VCF file format on NCBI
import re, sys, os


def change_fasta(filename):
    matched = False
    with open(filename) as fh, open(filename+'.tmp', 'w') as wfile:
        for line in fh:
            if line.startswith('>NC_'):
                match = re.findall('chromosome ([0-9ABXY]+),', line)
                if match:
                    wfile.write('>'+match[0]+'\n')
                    matched = True
                else:
                    matched = False
            elif line.startswith('>'):
                matched = False
            elif matched:
                wfile.write(line)
        wfile.close()
        fh.close()
    os.system('rm {0}; mv {0}.tmp {0}'.format(filename))

if __name__ == '__main__':
    filename = sys.argv[1:][0]
    change_fasta(filename)




