# this will take coordinates and a chromosome and stdout variant information


import gzip, sys, getopt
def help():
    return 'usage: \n python search_vcf_coords.py -c <chromosome_name> -f <five-prime_position> -t <three-prime_position> <VCF_file> > <out_file>'

def find_variants(chr = 'chr5', fivep = '0', threep = '100', vcf_file = 'vcf'):
    fivep, threep = int(fivep), int(threep)
    if vcf_file[-3:] == '.gz':
        openfn = gzip.open
    else:
        openfn = open
    with openfn(vcf_file,'rb') as vcfF:
        sys.stdout.write('## This file was produced using "search_vcf_coords.py"\n##found variants in {0}: {1} {2} - {3}\n'.format(vcf_file,chr,fivep,threep))
        sys.stdout.write('#CHROM\tPOS\tID\tREF\tALT\n')
        for line in vcfF:
            if line.startswith('#'):
                continue
            elif line.startswith(chr):
                line = line.strip().split('\t')
                line[1] = int(line[1])
                if line[1] >= fivep:
                    if line[1] <= threep:
                        sys.stdout.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(*line[0:5]))
                    elif line[1] > threep: break
    sys.stderr.write("finished\n\n")

if __name__ == '__main__':
    args = sys.argv[1:]
    opts, args = getopt.getopt(args,'c:f:t:')
    optdict = dict(opts)
    if not '-c' in optdict or not '-f' in optdict or not '-t' in optdict:
        print help()
    else:
        find_variants(chr = optdict['-c'], fivep = optdict['-f'], threep = optdict['-t'], vcf_file = args[0])


