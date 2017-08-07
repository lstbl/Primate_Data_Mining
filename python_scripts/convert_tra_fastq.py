# this function will take 3 files, the fasta file with the trace reads, the file with quality scores, and
# the info file that indicates where the 5' and 3' trimming should be. I will initiall put all capillary
# reads into a single fastq file and map that to the reference ape genome.

##### NOTE: THIS PROGRAM IS DESIGNED TO TAKE IN GZIPPED FILES; it will output non-gzipped fastq files #####


import sys, getopt, gzip
def help():
    print "Usage: \n python convert_tra-fastq.py -o output_filename <fasta_file> <quality_score_file> <clip_file> \n\n"

# this function will take the string values corresponding to the clipping categories and output either the integer-
# converted values or 0 (for left boarder) and 9999999999999 (for right boarder) if there are no values assigned.
def get_clipping(qual):
    best_interval = (0,0)
    current_interval = []
    for n in range(len(qual)):
        if ord(qual[n]) >= 20:
            current_interval.append(n)
        if ord(qual[n]) < 20:
            if len(current_interval)>0:
                if best_interval[1]-best_interval[0] < current_interval[-1]-current_interval[0]:
                    best_interval = (current_interval[0],current_interval[-1])
            current_interval = []
    return best_interval


# create fastq (PHRED33 scaled) from fasta file, qual file, and info (anc) file:
def convert_to_fastq(output_filename, fasta_file, qual_file, clip_file):
    with gzip.open(fasta_file,'rb') as fastaF, gzip.open(qual_file,'rb') as qualF, gzip.open(clip_file, 'rb') as clipF, \
            open(output_filename,'wb') as wfile:
        seqname = fastaF.readline().strip()
        qualname = qualF.readline().strip()
        clipF.readline().strip().split('\t')
        clipping = clipF.readline().strip().split('\t')
        while seqname and qualname:
            try:
                assert seqname == qualname
            except:
                print "Sequence header not identical to quality header"
                raise AssertionError
            sequence = ''
            qual = ''
            seqbuffer = True
            qualbuffer = True
            while seqbuffer:
                seqbuffer = fastaF.readline().strip()
                if seqbuffer.startswith('>'):
                    break
                else:
                    sequence += seqbuffer
            newseqname = seqbuffer
            while qualbuffer:
                qualbuffer = qualF.readline().strip()
                if qualbuffer.startswith('>'):
                    break
                else:
                    qual += ''.join(chr(33+int(_)) for _ in qualbuffer.split())
            newqualname = qualbuffer
            assert len(sequence) == len(qual)
            if clipping:
                if clipping[0] == seqname.split()[0].split('|')[-1]:
                    clip_left, clip_right = int(clipping[1]), int(clipping[2])
                else:
                    clip_left, clip_right = get_clipping(qual)
            sequence, qual = sequence[clip_left:clip_right],qual[clip_left:clip_right]
            wfile.write('@{0}\n{1}\n+\n{2}\n'.format(seqname[1:],sequence,qual))
            seqname, qualname = newseqname, newqualname
            clipping = clipF.readline().strip().split('\t')
    return

if __name__ == '__main__':
    args = sys.argv[1:]
    opts, files = getopt.getopt(args,'o:')
    optsdict = dict(opts)
    convert_to_fastq(optsdict['-o'],*files)

