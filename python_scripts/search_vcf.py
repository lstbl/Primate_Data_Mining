# this will take coordinates and a chromosome and output a file with variant information


import gzip, sys, getopt, re
from collections import defaultdict
def help():
    return 'usage: \n python search_vcf.py -g <gene name (NM_XXXX refseq designation)> <VCF_file> <GFF_file> <chromosome_fasta_file_path> > outfile.vcf'

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

def get_variant_info(gff_file = 'hg18_refseq.gtf.gz',gene_id = None, vcf_file = 'Pan_paniscus', fasta_path = './chr_fastas'):
#    variant_files = ['Gorilla','Homo','Pan_paniscus','Pan_troglodytes','Pongo_abelii','Pongo_pygmaeus']
#    if vcf_file not in variant_files:
#        print 'this is not a valid organism name, please input one of the following: {0}'.format(', '.join(_ for _ in variant_files))
#        return help()
    CDSboundaries = []
    chrom = ''
    strand = ''
    print 'searching gff file...\n'
    if gff_file[-3:] == '.gz':
        openfn = gzip.open
    else:
        openfn = open
    with openfn(gff_file) as gff:
        for line in gff:
            if re.findall(gene_id,line):
                line = line.strip().split()
                if line[2] == 'CDS':
                    CDSboundaries.append((int(line[3]),int(line[4])))
                if not chrom:
                    chrom = line[0].strip()
                else:
                    assert line[0] == chrom
                if not strand:
                    strand = line[6]
                else:
                    assert strand == line[6]
    CDSboundaries = sorted(CDSboundaries)
    CDS_positions = [x for _ in CDSboundaries for x in range(_[0],_[1]+1)]
    sys.stderr.write('found gene location on {0}\n'.format(chrom))
    start, end = CDSboundaries[0][0], CDSboundaries[-1][-1]
#    genomic_positions = range(start,end+1)
    genomic_seq = get_genomic_sequence(boundaries= CDSboundaries, chrom= chrom, fasta_path = fasta_path)
    sys.stderr.write('searching {0} fasta file...\n'.format(chrom))
    coding_seq = ''.join(genomic_seq[_[0]-start:_[1]-start+1] for _ in CDSboundaries)
    sys.stderr.write('found coding sequence for {0}\n'.format(gene_id))
    assert len(coding_seq) == len(CDS_positions)
#    vcf_file = '{0}/{1}.vcf.gz'.format(vcf_file,chrom)
    sys.stderr.write('searching {0} for variants... this may take a while'.format(vcf_file))
    variant_dic = find_variants(chr= chrom, CDS= CDS_positions, vcf_file= vcf_file, gene_name= gene_id)
    sys.stdout.write('#{0}\t{1}\t{2}\t{3}\n'.format('AA position','human_ref','mutant','alt_allele_count(AC);alt_allele_freq(AF);total_num_allels(AN)'))
    for pos, basenum in enumerate(CDS_positions):
        if basenum in variant_dic:
            for entry in variant_dic[basenum]:
                codon, mutant = coding_seq[3*(pos/3):3*(pos/3)+3],[_ for _ in coding_seq[3*(pos/3):3*(pos/3)+3]]
                mutant[pos%3] = entry[0][1]
                mutant = ''.join(_ for _ in mutant)
                if strand == '-':
                    codon, mutant, codon_pos = revcomp(codon), revcomp(mutant), ((len(coding_seq)-pos)/3)+1 # add one because of zero-indexing in python
                else:
                    codon_pos = (pos/3) + 1 # add one because of zero-indexing in python
                sys.stdout.write('{0}\t{1}({2})\t{3}({4})\t{5}\n'.format(codon_pos,codon,code[codon],mutant,code[mutant],'AC={0};AF={1};AN={2}'.format(*entry[1])))

    print '\ndone! stored intron and exon variant information (hg18 coordinate system) in {0}.vcf\n'.format(gene_id)

# this will not return sequence reverse-complemented!
def get_genomic_sequence(fasta_path = "hg18_chroms", boundaries = None, chrom = None):
    start, end = int(boundaries[0][0]), int(boundaries[-1][-1])
    with open(fasta_path+'/'+chrom+'.fa') as ofile:
        current_chr = ofile.readline().strip()[1:]
        pos = 0
        seq = ''
        # find the start of the sequence
        for line in ofile:
            line = line.strip()
            if pos + len(line) >= start:
                start = start - pos - 1
                # make sure the end isn't in the same line (most likely never will be)
                if pos + len(line) >= end:
                    end = end - pos
                    return line[start:end]
                seq += line[start:]
                pos += len(line)
                # break out of loop once we've found the start and restart the iteration to find the end
                break
            else:
                pos += len(line)
        # find the end of the sequence
        for line in ofile:
            line = line.strip()
            if pos + len(line) >= end:
                seq += line[:end - pos]
                return seq.upper()
            seq += line
            pos += len(line)
#    coding_seq = ''.join(seq[_[0] - start:_[1]-start+1] for _ in boundaries).upper()
#    return coding_seq


def revcomp(seq):
    revdic = {'A':'T','T':'A','C':'G','G':'C'}
    return ''.join(revdic[_] for _ in seq[::-1])


def find_variants(chr = 'chr5', CDS = [0,0], vcf_file = 'vcf', gene_name = 'out'):
    fivep, threep = CDS[0], CDS[-1]
    CDS = set(CDS)
    variant_dic = defaultdict(list)
    if vcf_file[-3:] == '.gz':
        openfn = gzip.open
    else:
        openfn = open
    with openfn(vcf_file,'rb') as vcfF,  open(gene_name+'.vcf','w') as wfile:
        wfile.write('## This file was produced using "search_vcf.py"\n##found variants in {0}: {1} {2} - {3}\n'.format(vcf_file,chr,fivep,threep))
        wfile.write('#CHROM\tPOS\tID\tREF\tALT\tAlt_N;Alt_freq;N_alleles\n')
        for line in vcfF:
            if line.startswith('#'):
                continue
            elif line.startswith(chr):
                line = line.strip().split('\t')
                line[1] = int(line[1])
                if line[1] >= fivep:
                    if line[1] <= threep:
                        ACAFAN = re.findall('A[CFN]=([^;]+)', line[7])
                        wfile.write('{0}\t{1}\t{2}\t{3}\t{4}\tAC={5};AF={6};AN={7}\n'.format(*line[0:5]+ACAFAN))
                        if int(line[1]) in CDS:
                            variant_dic[int(line[1])].append(((line[3],line[4]), (int(ACAFAN[0]),float(ACAFAN[1]),int(ACAFAN[2]))))
                    elif line[1] > threep: break
    return variant_dic

if __name__ == '__main__':
    args = sys.argv[1:]
    sys.stderr.write(str(args)+'\n\n')
    opts, args = getopt.getopt(args,'g:')
    optdict = dict(opts)
    if not '-g' in optdict:
        print help()
    else:
        get_variant_info(gene_id= optdict['-g'], vcf_file= args[0] ,gff_file=args[1], fasta_path = args[2])


