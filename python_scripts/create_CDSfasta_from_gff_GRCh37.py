# This will take a gff file and genome fasta file and create a fasta file with all CDSs. This will also put the CDS
# boarders in the name. The fasta file will be built as such:
# ><gene_name>|<NCBI gene ID>|<refseq accession number>|<chromosome>|<strand>|<CDS boundaries>
# for example:
# >Apobec3G|449577|NP_001009001|chr22|+|25610-25615,256640-256724,256801-256950
# ATG...
# NOTE: bondaries are 1-based and inclusive
# NOTE: this is modified to handle the GRCh37 annotation, which does not include the "gene=" designation in the CDS
# entries


import gzip
import re
from pyfaidx import Fasta
import sys

def helpmessage():
    sys.stderr.write('\nUSE AS FOLLOWS:\npython <genome_fasta_file> > <output_file>\n NOTE: there must be a indexed fasta file and gff file with same prefix in the folder\n\n')

def RevComp(sequence):
    RCDict = {'A':'T','G':'C','C':'G','T':'A','N':'N'}
    sequence = [RCDict[_] for _ in sequence]
    return ''.join(_ for _ in sequence[::-1])

def collect_CDSs(gff_file):
    gene_name_dic = {}
    gff_info = {}
    gene_order = []
    geneID_set = set()
    if gff_file[-3:] == '.gz':
        openfn = gzip.open
    else:
        openfn = open
    with openfn(gff_file) as fh:
        for line in fh:
            if line.startswith('#'): continue
            line = line.strip().split('\t')
            if line[2] == 'gene':
                gene = re.findall('gene=([^;\n]+)',line[-1])[0]
                GeneID = re.findall('GeneID:([0-9]+)', line[-1])[0]
                gene_name_dic[GeneID] = gene
            if line[2] == 'CDS':
                parent = re.findall('Parent=([^;]+)',line[-1])
                genbank = re.findall('Genbank:([^;,]+)',line[-1])
                GeneID = re.findall('GeneID:([0-9]+)',line[-1])
                if GeneID:
                    assert GeneID[0] in gene_name_dic
                gene = gene_name_dic[GeneID[0]]
                start, end = int(line[3]), int(line[4])
                strand = line[6]
                chrom = line[0]
                if parent and genbank and GeneID and gene:
                    parent, genbank, GeneID = parent[0], genbank[0], GeneID[0]
                    if GeneID not in geneID_set:
                        gene_order.append(GeneID)
                        geneID_set.add(GeneID)
                    if GeneID not in gff_info:
                        gff_info[GeneID] = {'CDS':{genbank:{'strand':strand,'chrom':chrom,'parent':parent,'boundaries':[(start,end),]}},'gene':gene}
                    elif genbank not in gff_info[GeneID]['CDS']:
                        assert gene == gff_info[GeneID]['gene']
                        gff_info[GeneID]['CDS'][genbank] = {'strand':strand,'chrom':chrom,'parent':parent,'boundaries':[(start,end),]}
                    else:
                        if parent != gff_info[GeneID]['CDS'][genbank]['parent']: continue
                        assert gene == gff_info[GeneID]['gene'] and strand == gff_info[GeneID]['CDS'][genbank]['strand'] and chrom == gff_info[GeneID]['CDS'][genbank]['chrom']
                        gff_info[GeneID]['CDS'][genbank]['boundaries'].append((start,end))
    return gff_info, gene_order

def write_CDS_fasta(fasta_file, gff_file):
    gff_info, gene_order = collect_CDSs(gff_file)
    GenomeReader = Fasta(fasta_file,sequence_always_upper=True)
    for GeneID in gene_order:
        for genbank in gff_info[GeneID]['CDS']:
            chrom, strand = gff_info[GeneID]['CDS'][genbank]['chrom'],gff_info[GeneID]['CDS'][genbank]['strand']
            seq = ''
            gff_info[GeneID]['CDS'][genbank]['boundaries'] = sorted(gff_info[GeneID]['CDS'][genbank]['boundaries'])
            for boundary in gff_info[GeneID]['CDS'][genbank]['boundaries']:
                seq += GenomeReader[chrom][boundary[0]-1:boundary[1]].seq
            if strand == '-':
                seq = RevComp(seq)
            if len(seq)%3 != 0:
                sys.stderr.write('genebank accession {0} of gene {1} is not divisible by 3, skipping\n'.format(genbank,gff_info[GeneID]['gene']))
                continue
            if seq[:3] != 'ATG':
                sys.stderr.write('genebank accession {0} of gene {1} has a non-canonical start site ({2})\n'.format(genbank,gff_info[GeneID]['gene'],seq[:3]))
            # >Apobec3G|449577|NP_001009001|chr22|+|25610-25615,256640-256724,256801-256950
            sys.stdout.write('>{0}|{1}|{2}|{3}|{4}|{5}\n'.format(gff_info[GeneID]['gene'],GeneID,genbank,chrom,strand,','.join(str(_[0])+'-'+str(_[1]) for _ in gff_info[GeneID]['CDS'][genbank]['boundaries'])))
            for pos in xrange(0,len(seq),50):
                sys.stdout.write(seq[pos:pos+50]+'\n')

if __name__ == '__main__':
    try:
        fasta_file = sys.argv[1:][0]
    except:
        helpmessage()
    if fasta_file[-3:] != '.fa':
        sys.stderr.write('must be a ".fa" file')
    prefix = fasta_file[:-3]
    write_CDS_fasta(fasta_file,prefix+'.gff')


