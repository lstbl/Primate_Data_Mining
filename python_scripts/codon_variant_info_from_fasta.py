# this will take variant info from fasta files produced from the "VCF_to_genotypes.py" program and calculate
# basic statistics such as those found when using the "search_vcf.py" program
# the output will be as such:
#
# #AA position  ref   mutant    alt_allele_count(AC);alt_allele_freq(AF);total_num_alleles(AN)  subspecies1 allele count    subpecies2 allele count ...etc

import glob
from collections import defaultdict
from itertools import izip
import sys

def helpstatement():
    sys.stderr.write('python codon_variant_info_from_fastas.py <file_prefix>')

code = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
        'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
        'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
        'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
        'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
        'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
        'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
        'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
        'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
        'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
        'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
        'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
        'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
        'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
        'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
        'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
        }

REVCOMP = {'A':'T','T':'A','C':'G','G':'C','Y':'R','R':'Y','W':'W','S':'S','K':'M','M':'K','N':'N','-':'-'}
AMBIGUOUS = {('C', 'T'): 'Y', ('A', 'G'): 'R', ('A', 'T'): 'W', ('G', 'C'): 'S', ('T', 'G'): 'K', ('C', 'A'): 'M',
             ('T', 'C'): 'Y', ('G', 'A'): 'R', ('T', 'A'): 'W', ('C', 'G'): 'S', ('G', 'T'): 'K', ('A', 'C'): 'M',
             ('A','*'):'A',('*','A'):'A',('C','*'):'C',('*','C'):'C',('T','*'):'T',('*','T'):'T',('G','*'):'G',('*','G'):'G',
             ('A','A'):'A',('C','C'):'C',('T','T'):'T',('G','G'):'G',('*','*'):'-'}
REV_AMBIGUOUS = {'A': ('A', 'A'), 'C': ('C', 'C'), 'G': ('G', 'G'), 'K': ('T', 'G'), 'M': ('C', 'A'),
                 'R': ('G', 'A'), 'S': ('C', 'G'), 'T': ('T', 'T'), 'W': ('T', 'A'), 'Y': ('C', 'T'),
                 'N': ('N','N')}


def SUMMARY_STATS(genbank_accession):
    files = glob.glob(genbank_accession+'*')
    assert genbank_accession+'_reference.fa' in files
    files.remove(genbank_accession+'_reference.fa')
    with open(genbank_accession+'_reference.fa') as fh:
        refseq = ''
        for line in fh:
            if line.startswith('>'):
                line = line.strip()[1:].split('|')
                boarders = [tuple(int(__) for __ in _.split('-')) for _ in line[-1].split(',')]
                GeneName = line[0]
                GenbankAccession = line[2]
                GeneID = line[1]
            else:
                refseq += line.strip()
    Subspecies_Alleles = defaultdict(list)
    for file in files:
        with open(file) as fh:
            current_seq = ''
            for line in fh:
                if line.startswith('#'): break
                if line.startswith('>'):
                    line = line.strip()[1:].split('|')
                    current_species = line[-1]
                    if current_seq:
                        Subspecies_Alleles[current_species].append(current_seq)
                        Subspecies_Alleles['ALL'].append(current_seq)
                        current_seq = ''
                else:
                    current_seq += line.strip()
            if current_seq:
                Subspecies_Alleles[current_species].append(current_seq)
                Subspecies_Alleles['ALL'].append(current_seq)
    # create a dictionary with summary statistics
    SummaryStats = {}
    for pos, bases in enumerate(izip(*[refseq]+Subspecies_Alleles['ALL'])):
        if TEST_FOR_DIFFERENCES(bases):
            codon_pos = pos/3
            refcodon = refseq[3*(codon_pos):3*(codon_pos)+3]
            SummaryStats[codon_pos] = dict((_,{}) for _ in Subspecies_Alleles.keys())
            for sample in Subspecies_Alleles:
                SummaryStats[codon_pos][sample] = GET_STATS([_[pos] for _ in  Subspecies_Alleles[sample]],refcodon,pos%3)

    SpeciesList = [_ for _ in Subspecies_Alleles if _ != "ALL"]
    sys.stdout.write('#AA_position\tREF\tALT\tAlleleCount(AC),AlleleFrequence(AF),AlleleNumber(AN)\t{0}\n'.format(
        '\t'.join(_.replace(' ', '_') for _ in SpeciesList)))

    Variable_Positions = sorted(SummaryStats.keys())
    for position in Variable_Positions:
        refcodon = refseq[position*3:position*3+3]
        for altcodon in SummaryStats[position]['ALL']:
            if altcodon != 'AN':
                sys.stdout.write('{0}\t{1}({2})\t{3}({4})\t{5}\t{6}\n'.
                                 format(position+1,refcodon,code.get(refcodon,''),altcodon,code.get(altcodon,''),
                                        PRINT_STATS(altcodon,SummaryStats[position]['ALL']),
                                        '\t'.join(PRINT_STATS(altcodon,SummaryStats[position][_]) for _ in SpeciesList)))
    return True


def PRINT_STATS(AltCodon,InfoDic):
    if AltCodon in InfoDic:
        return 'AC={0},AF={1},AN={2}'.format(InfoDic[AltCodon]['AC'],round(InfoDic[AltCodon]['AF'],3),InfoDic['AN'])
    else:
        return 'AC={0},AF={1},AN={2}'.format(0,0,InfoDic['AN'])


def GET_STATS(bases, refcodon, position):
    refseq = refcodon[position]
    AN = 0
    AC = defaultdict(int)
    for base in bases:
        if base != 'N':
            AN += 2
            altcodon = refcodon[:position]+base+refcodon[position+1:]
            genotypes = SPLIT_CODONS(altcodon, refcodon)
            for genotype in genotypes:
                if genotype != refcodon:
                    AC[genotype] += 1
    SummaryStats = {'AN':AN}

    for genotype in AC:
        SummaryStats[genotype] = {'AC':AC[genotype],'AF':ALLELE_FREQ(AC[genotype],AN)}
    return SummaryStats

def ALLELE_FREQ(AC,AN):
    if AN != 0:
        if AC <= AN:
            return float(AC)/AN
    return 'NA'


def TEST_FOR_DIFFERENCES(bases):
    refseq = bases[0]
    for i in bases[1:]:
        if i != refseq:
            return True
    return False


def SPLIT_CODONS(codon, ref):
    differences = []
    for pos, (base1, base2) in enumerate(zip(codon,ref)):
        if base1 != base2:
            differences.append(pos)
    if len(differences) >1:
        # can't figure out exact genotype for this codon given this data, just return variant with ambigous bases (or homozygous bases)
        return (codon,codon)
    elif len(differences) == 0:
        return (ref,ref)
    else:
        assert len(differences) == 1
        difpos = differences[0]
        differences = REV_AMBIGUOUS[codon[differences[0]]]
        codon1, codon2 = [_ for _ in codon], [_ for _ in codon]
        codon1[difpos], codon2[difpos] = differences[0], differences[1]
        return (''.join(_ for _ in codon1), ''.join(_ for _ in codon2))


    # function that takes codon and returns a 2-tuple with the variants of that codon (if it's homozygous, return a 2-tuple of the same entry)

if __name__ == '__main__':
    arg = sys.argv[1:][0]
    exit = SUMMARY_STATS(arg)
    if not exit:
        helpstatement()



