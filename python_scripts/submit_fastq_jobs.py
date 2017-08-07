# This python script will modify pbs scripts and submit jobs using the 'convert_tra-fastq.py' program.
# use in the directory with all of the fasta, qual, and anc files downloaded from the trace read archive.
# it will use glob to grab all the appropriate filenames and modify the pbs script accordingly

import glob
import re
import os
def modify_script(suffix):
    pbs_script = """#PBS -m ae
#PBS -M alst5940@colorado.edu
#PBS -N {0}
#PBS -q short
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=8gb
MNTDIR='/mnt/scratch_nobackup/alst5940/Primate_data_mining'
cd $MNTDIR/FASTQ/tra_reads/Pan_troglodytes
module load python_2.7.3
python /mnt/scratch_nobackup/alst5940/Primate_data_mining/python_scripts/convert_tra_fastq.py -o {0}.fastq fasta.{0} qual.{0} clip.{0} 2> {0}.err""".format(suffix)
    with open('/mnt/scratch_nobackup/alst5940/Primate_data_mining/PBS/fastq_convert.pbs','w') as wfile:
        wfile.write(pbs_script)
        wfile.close()
    return

if __name__ == '__main__':
#    suffixes = sorted(list(set([re.findall('(anc|fasta|qual)\.(.*.gz)',entry)[0][1] for entry in glob.glob('*pan_troglodytes_schweinfurthii*')])))
#    suffixes = sorted(list(set([re.findall('(clip|fasta|qual)\.(.*\.gz)',entry)[0][1] for entry in glob.glob('*')])))
    suffixes = ['pan_troglodytes.002.gz','pan_troglodytes.003.gz','pan_troglodytes.050.gz','pan_troglodytes_schweinfurthii.001.gz']
    for suffix in suffixes:
        modify_script(suffix)
        os.system('qsub /mnt/scratch_nobackup/alst5940/Primate_data_mining/PBS/fastq_convert.pbs')
