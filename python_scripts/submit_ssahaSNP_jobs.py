# This python script will modify pbs scripts and submit jobs using the 'convert_tra-fastq.py' program.
# use in the directory with all of the fasta, qual, and anc files downloaded from the trace read archive.
# it will use glob to grab all the appropriate filenames and modify the pbs script accordingly

import glob
import re
import os
def modify_script(prefix):
    pbs_script = """#PBS -m ae
#PBS -M alst5940@colorado.edu
#PBS -N ssahaSNP{0}
#PBS -q long
#PBS -l walltime=200:00:00
#PBS -e localhost:$MNTDIR/PBS/eo_files/ssahaSNP.err
#PBS -o localhost:$MNTDIR/PBS/eo_files/ssahaSNP.log
#PBS -l nodes=1:ppn=2
#PBS -l mem=16gb
MNTDIR="/mnt/scratch_nobackup/alst5940/Primate_data_mining"
FASTADIR="$MNTDIR/GENOMES/Pan_troglodytes/GCF_000001515.5.hash"
FASTQDIR="$MNTDIR/FASTQ/tra_reads/Pan_troglodytes/{0}.fastq"
OUTFILEDIR="$MNTDIR/knownSNPs/Pan_troglodytes"
EXEDIR="$MNTDIR/ssaha2/ssahaSNP"
cd $OUTFILEDIR
$EXEDIR -save $FASTADIR $FASTQDIR > {0}_ssahaSNP.ssaha2 2> {0}.err""".format(prefix)
    with open('/mnt/scratch_nobackup/alst5940/Primate_data_mining/PBS/ssahaSNP.pbs','w') as wfile:
        wfile.write(pbs_script)
        wfile.close()
    return

if __name__ == '__main__':
    prefixes = sorted([re.findall('(.*)\.fastq',entry)[0] for entry in glob.glob('*.fastq')])
    for prefix in prefixes:
        modify_script(prefix)
        os.system('qsub /mnt/scratch_nobackup/alst5940/Primate_data_mining/PBS/ssahaSNP.pbs')

