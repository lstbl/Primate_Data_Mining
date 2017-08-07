# this will take information from "hash_fasta.py" from the genome which short reads were mapped to and
# the genome that the annotation file is in (usually the NCBI genome). It will first create a dictionary between the NCBI
# chromosome nomenclature and the reference nomenclature and then alter the gff file to reflect this change.
# run as such:
# python convert_NCBI_gff.py <hash1> <hash2> <gff_file>

import sys
import os
import gzip

def help_message():
    sys.stderr.write('python convert_NCBI_gff.py <hash1> <hash2> <gff_file>')

def parsehashfile(open_file):
    hash_dic = {}
    for line in open_file:
        if line.startswith('chr_name'): continue
        line = line.strip().split('\t')
        if line:
            chrname, length, chrhash = line
            hash_dic[chrhash] = chrname
    return hash_dic

def mapping_dic(hash1,hash2):
    with open(hash1) as fh1, open(hash2) as fh2:
        hash1_dic, hash2_dic = parsehashfile(fh1), parsehashfile(fh2)
    mapping_dictionary = {}
    for entry in hash1_dic:
        if entry in hash2_dic:
            mapping_dictionary[hash1_dic[entry]] = hash2_dic[entry]
            mapping_dictionary[hash2_dic[entry]] = hash1_dic[entry]
    return mapping_dictionary

def modify_gff(gff_file, hash1, hash2):
    chr_map = mapping_dic(hash1,hash2)
    if gff_file[-3:] == '.gz':
        openfn = gzip.open
        newname = gff_file[:-3]
    else:
        openfn = open
        newname = gff_file
    with openfn(gff_file) as fh, open(newname+'.tmp','w') as wfile:
        for line in fh:
            if line.startswith('#'):
                wfile.write(line)
            else:
                line = line.split('\t')
                if line[0] in chr_map:
                    line[0] = chr_map[line[0]]
                    wfile.write('\t'.join(_ for _ in line))
    os.system('rm {0} && mv {1}.tmp {1}'.format(gff_file,newname))


if __name__ == "__main__":
    try:
        hash1, hash2, gff_file = sys.argv[1:]
    except:
        help_message()
    modify_gff(hash1=hash1, hash2=hash2, gff_file=gff_file)


