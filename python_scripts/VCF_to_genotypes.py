# This file will take as input the VCF file, a ".fa.vcfcds" file, and a file that lists the NCBI genbank protein
# accession number of the genes you are interested in and will output a fasta file with genotype information.
# NOTE: the .fa.vcfcds file MUST have the same prefix as the fasta file that corresponds to the genome assembly.

# USAGE: python VCF_to_genotypes.py <VCF_file> <genome Fasta file> <file with genbank accession numbers>

# NOTE: file with genbank accession numbers must have a hard return between the entries. For example:
# XP_31135213.1
# NP_44821784.4
#
# etc
import sys
import os
from collections import defaultdict

def helpstatement():
    sys.stderr.write('python VCF_to_genotypes.py <VCF_file> <genome Fasta file> <file with genbank accession numbers>\n\n')

def get_genbank_accessions(genbank_accession_file):
    accession_numbers = set()
    with open(genbank_accession_file) as fh:
        for line in fh:
            line = line.strip()
            if line:
                accession_numbers.add(line.split('.')[0])
    return accession_numbers



def find_accessions(vcfcds_file, accession_number):
    seq = ''
    GeneID = ''
    genbank_info = {}
    Genbank = ''
    with open(vcfcds_file) as fh:
        for line in fh:
            if line.startswith('genebank accession'): continue
            line = line.strip()
            if (not line and seq) or line.startswith('>'):
                if seq and Genbank == accession_number:
                    genbank_info[Genbank] = {'chrom':chrom,'strand':strand,'GeneName':GeneName[1:],'GeneID':GeneID,'boarders':boarders,'seq':seq}
                    break
                    #seq = ''
                    #GeneID = ''
                if line:
                    line = line.split('|')
                    GeneName, GeneID, Genbank, chrom, strand, boarders = line
                    Genbank = Genbank.split('.')[0]
                    seq = ''
                    if Genbank != accession_number: continue
                    boarders = [tuple(int(j) for j in _.split('-')) for _ in boarders.split(',')]
                    boarders = [x for _ in boarders for x in range(_[0],_[1]+1)]
            elif line:
                seq += line
        if seq and Genbank == accession_number:
            genbank_info[Genbank] = {'chrom': chrom, 'strand': strand, 'GeneName': GeneName[1:],'GeneID':GeneID,'boarders': boarders,'seq':seq}
    return genbank_info

#generate a indexed VCF file for quick lookups. This file will be organized by chromosome.
# The format is as follows:
# >chr1
# 0 <offset in bytes from start of file>
# 1000000 <offset in bytes for all entries >= 1000000>
# etc
def make_VCF_index(VCF_file):
    current_pos = 0
    current_chr = ''
    VCF_index = VCF_file+'.pos.dat'
    if not os.path.isfile(VCF_index):
        with open(VCF_file) as fh, open(VCF_index,'w') as wfile:
            while True:
                filepos = str(fh.tell())
                line = fh.readline().strip()
                if not line: break
                if line.startswith('#'): continue
                line = line.split('\t')[:2]
                if line[0] != current_chr:
                    current_pos = 0
                    wfile.write('>'+line[0]+'\n')
                    wfile.write(str(current_pos)+' '+filepos+'\n')
                    current_chr = line[0]
                    current_pos = 100000
                while current_pos <= int(line[1]):
                    wfile.write(str(current_pos)+' '+filepos+'\n')
                    current_pos += 100000
    VCFINDEX = defaultdict(list)
    with open(VCF_index) as fh:
        for line in fh:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                current_chr = line[1:]
            else:
                line = line.split()
                VCFINDEX[current_chr].append(int(line[-1]))
    return VCFINDEX


#function to get variant information from VCF file:
def get_variant_info(VCF_file, vcfcds_file, accession_number):
    VCFINDEX = make_VCF_index(VCF_file)
    gene_info = find_accessions(vcfcds_file, accession_number)
    if accession_number not in gene_info:
        gene_info = find_accessions('/'.join(_ for _ in vcfcds_file.split('/')[:-1])+'/vcfcds.stderr',accession_number)
        if accession_number in gene_info:
            sys.stderr.write('gene {1} ({0}) is out of frame due to issue with genome fasta file, keep in mind when interpreting results\n'.format(accession_number,new_gene_info[accession_number]['GeneName']))
        else:
            sys.stderr.write('gene {0} not found in vcfcds file or stderr file\n'.format(accession_number))
            return None
    PASS = False
    with open(VCF_file) as fh:
        for line in fh:
            if line.startswith('#CHR'):
                BIOSAMPS = line.strip().split('\t')[9:]
                break
        sys.stderr.write('testing gene {0}\n'.format(gene_info[accession_number]['GeneName']))
        fh.seek(0, 0)
        chr, boarders, strand, start, end, seq = gene_info[accession_number]['chrom'], gene_info[accession_number]['boarders'], \
                                                 gene_info[accession_number]['strand'], int(gene_info[accession_number]['boarders'][0]), \
                                                 int(gene_info[accession_number]['boarders'][-1]), gene_info[accession_number]['seq']
        #check to see if this SNP could even be in the VCF file (if the position is greater than anything found in the VCF)
        if start/100000 >= len(VCFINDEX[chr]):
            sys.stderr.write('gene {0} not in VCF file'.format(gene_info[accession_number]['GeneName']))
            return None

        #check to see if gene being tested is on a chromosome that was called in the VCF file
        if chr not in VCFINDEX:
            sys.stderr.write('chr {0} not found in VCF file, skipping gene {1}'.format(chr,accession_number))
            return None

        #if the gene is on reverse strand, revcomp the gene (this will be corrected later
        if strand == '-':
            seq = revcomp(seq)
        sequences = [[_ for _ in seq] for __ in range(len(BIOSAMPS))]

        #seek to correct spot in VCF file:
        fh.seek(VCFINDEX[chr][start/100000],0)

        while not PASS:
            line = fh.readline().strip().split('\t')
            if line[6] != '.':
                PASS = True
        CHROM,POS,ID,REF,ALT,QUAL  = line[:6]
        GENOTYPES = line[9:]
        POS = int(POS)
        pos, current = 0, boarders[0]
        while POS <= end:
            while POS < current:
                PASS = False
                while not PASS:
                    line = fh.readline().strip().split('\t')
                    assert line
                    if line[6] != '.':
                        PASS = True
                CHROM, POS, ID, REF, ALT, QUAL = line[:6]
                POS = int(POS)
                GENOTYPES = line[9:]
            while POS > current and pos<len(boarders)-1:
                pos += 1
                current = boarders[pos]
            if POS == current:
                GENOTYPES = get_genotypes(GENOTYPES,REF,ALT)
                for N,sequence in enumerate(sequences):
                    sequence[pos] = GENOTYPES[N]
                PASS = False
                while not PASS:
                    line = fh.readline().strip().split('\t')
                    assert line
                    if line[6] != '.':
                        PASS = True
                CHROM, POS, ID, REF, ALT, QUAL = line[:6]
                POS = int(POS)
                GENOTYPES = line[9:]

        sequences = [''.join(_ for _ in sequence) for sequence in sequences]
        if strand == '-':
            sequences = [revcomp(seq) for seq in sequences]
        gene_info[accession_number]['variants'] = sequences
    return gene_info, BIOSAMPS

def link_biosamples(BIOSAMPLE_INFO_FILE):
    SUBSPECIES = defaultdict(list)
    with open(BIOSAMPLE_INFO_FILE) as fh:
        for line in fh:
            line = line.strip().split(',')
            if line:
                SUBSPECIES[line[1]].append(line[0])
    return SUBSPECIES

def write_variant_info(gene_info,BIOSAMPS,SUBSPECIES_DESIGNATIONS,genome_fasta_file):
    if not gene_info:
        return None
    #this should only loop through once (changed code so each iteration is only a single gene)
    for gene in gene_info:
        with open(gene_info[gene]['GeneName']+'_'+gene+'_reference.fa','w') as wfile:
            boarders = collapse_boundaries(gene_info[gene]['boarders'])
            wfile.write('>{0}|{1}|{2}|{3}|{4}|{5}\n'.format(gene_info[gene]['GeneName'], gene_info[gene]['GeneID'], gene, gene_info[gene]['chrom'], genome_fasta_file.split('/')[-1][:-3],boarders))
            sequence = gene_info[gene]['seq']
            wfile.write('\n'.join(sequence[_:_ + 50] for _ in range(0, len(sequence), 50)) + '\n\n')
        for species in SUBSPECIES_DESIGNATIONS:
            with open(gene_info[gene]['GeneName']+'_'+gene+'_'+species.replace(' ','_')+'.fa','w') as wfile:
                for individual in SUBSPECIES_DESIGNATIONS[species]:
                    index = BIOSAMPS.index(individual)
                    sequence = gene_info[gene]['variants'][index]
                    wfile.write('>{0}|{1}|{2}|{3}|{4}\n'.format(gene_info[gene]['GeneName'],gene_info[gene]['GeneID'],gene,individual,species))
                    wfile.write('\n'.join(sequence[_:_+50] for _ in range(0,len(sequence),50))+'\n\n')
    return True

def collapse_boundaries(boudary_list):
    new_list = [boudary_list[0]]
    for pos,entry in enumerate(boudary_list[1:]):
        if entry-1 != boudary_list[pos]:
            new_list.append(boudary_list[pos])
            new_list.append(entry)
    new_list.append(entry)
    assert len(new_list)%2 == 0
    boarders = ''
    for n in range(0,len(new_list),2):
        boarders += '-'.join(str(_) for _ in new_list[n:n+2])+','
    return boarders[:-1]

#><assembly (ucsc designation)>|<gene name (NCBI annotation designation)>|<genbank accession>|<location of gene in genome assembly (inclusive)>


def revcomp(seq):
    return ''.join(REVCOMP[_] for _ in seq[::-1])

REVCOMP = {'A':'T','T':'A','C':'G','G':'C','Y':'R','R':'Y','W':'W','S':'S','K':'M','M':'K','N':'N','-':'-'}
AMBIGUOUS = {('C', 'T'): 'Y', ('A', 'G'): 'R', ('A', 'T'): 'W', ('G', 'C'): 'S', ('T', 'G'): 'K', ('C', 'A'): 'M',
             ('T', 'C'): 'Y', ('G', 'A'): 'R', ('T', 'A'): 'W', ('C', 'G'): 'S', ('G', 'T'): 'K', ('A', 'C'): 'M',
             ('A','*'):'A',('*','A'):'A',('C','*'):'C',('*','C'):'C',('T','*'):'T',('*','T'):'T',('G','*'):'G',('*','G'):'G',
             ('A','A'):'A',('C','C'):'C',('T','T'):'T',('G','G'):'G',('*','*'):'-'}


def get_genotypes(pos_9_to_end,REF,ALT):
    ALT = ALT.split(',')
    bases = [REF] + ALT
    genotypes = [_.split(':')[0].split('/') for _ in pos_9_to_end]

    for pos,genotype in enumerate(genotypes):
        if '.' in genotype:
            genotypes[pos] = 'N'
        else:
            BASE1, BASE2 = bases[int(genotypes[pos][0])],bases[int(genotypes[pos][1])]
            base_call = AMBIGUOUS[(BASE1,BASE2)]
            if base_call == '-':
                base_call = REF
            genotypes[pos] = base_call
    return genotypes


if __name__ == '__main__':
    exit = False
    VCF_file, GENOME_Fasta, Genbank_Accession_Number = sys.argv[1:]
    sys.stderr.write('{} {} {}\n'.format(VCF_file,GENOME_Fasta,Genbank_Accession_Number))
    vcfcds_file = GENOME_Fasta[:-3]+'.fa.vcfcds'
    gene_info, BIOSAMPS = get_variant_info(VCF_file, vcfcds_file, Genbank_Accession_Number)
    SUBSPECIES_DESIGNATIONS = link_biosamples(VCF_file[:-4]+'.INDIVIDUALS')
    exit = write_variant_info(gene_info, BIOSAMPS, SUBSPECIES_DESIGNATIONS,GENOME_Fasta)
    if not exit:
        helpstatement()


# python VCF_to_genotypes.py <VCF_file> <genome Fasta file> <file with genbank accession numbers>
