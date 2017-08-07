# This program will create a makefile pair for each biosample
# run as a for loop, e.g. for i in $BIOSAMPLES; do python create_makefiles.py $i; done

import sys, os
def write_makefiles_and_pbs(input):
    with open("../makefiles/"+input+"_fastqs.mk",'w') as wfile:
        wfile.write('''include gmsl
QUERY_INPUT={0}
QUERY_SRC=./pipeline.py
QUERY_EXE=python $(QUERY_SRC) --query $(QUERY_INPUT)
FASTQDUMP_EXE=/opt/sra/2.5.2/bin/fastq-dump --split-3 -F -I --gzip
SRADOWNLOAD_EXE=wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/
SCI_NAME=$(shell $(QUERY_EXE) --other_info=ScientificName)
FASTQ_DIR=../FASTQ/$(SCI_NAME)/$(QUERY_INPUT)
XRRs=$(shell $(QUERY_EXE) --XRRs)
SRA_FILES=$(addprefix $(FASTQ_DIR)/, $(addsuffix .sra, $(XRRs)))


## get_fastqs   : get all fastq files for particular query:
.PHONY : get_fastqs
.INTERMEDIATE : $(SRA_FILES)
get_fastqs : $(SRA_FILES)
%.sra :
	@mkdir -p $(dir $*); cd $(dir $*); if ls|grep -q $(notdir $*)[_12]*\.'fastq\|sra'; then : ; else $(SRADOWNLOAD_EXE)$(call substr,$(notdir $*),1,3)/$(call substr,$(notdir $*),1,6)/$(notdir $*)/$(notdir $*).sra &&  $(FASTQDUMP_EXE) $(notdir $*).sra; fi
'''.format(input))

    with open("../makefiles/"+input+"_Makefile.mk",'w') as wfile2:
        wfile2.write('''QUERY_INPUT={0}
QUERY_SRC=./pipeline.py
QUERY_EXE=python $(QUERY_SRC) --query $(QUERY_INPUT)
FASTQDUMP_EXE=/opt/sra/2.5.2/bin/fastq-dump --split-3 -F -I --gzip
SRADOWNLOAD_EXE=wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/
SCI_NAME=$(shell $(QUERY_EXE) --other_info=ScientificName)
FASTQ_DIR=../FASTQ/$(SCI_NAME)/$(QUERY_INPUT)
XRRs=$(shell $(QUERY_EXE) --XRRs)
SRA_FILES=$(addprefix $(FASTQ_DIR)/, $(addsuffix .sra, $(XRRs)))

GATK_EXE=../jdk1.8.0_91/bin/java -jar -Djava.io.tmpdir=`pwd`/tmp ../GenomeAnalysisTK.jar
PICARD_EXE=../jdk1.8.0_91/bin/java -jar -Djava.io.tmpdir=`pwd`/tmp ../picard-tools/picard.jar
SORTSAM_EXE=../jdk1.8.0_91/bin/java -Xmx10g -Djava.io.tmpdir=`pwd`/tmp -jar ../picard-tools/picard.jar SortSam
BWA_EXE=../bwa/bwa
SAMTOOLS_EXE=../samtools/samtools
BCFTOOLS_EXE=../bcftools/bcftools
NPROC=64
NPROC_BWA=8
REF_GENOME = $(shell $(QUERY_EXE) --ref_genome)
REF_GENOME_PATH=../GENOMES/$(REF_GENOME)/*.fa
PAIRED=$(shell ls $(FASTQ_DIR) | grep '_[12].fastq.gz' | sed 's/[_][12].fastq.gz//'|sort|uniq)
UNPAIRED=$(shell ls $(FASTQ_DIR) | grep '[^_12].fastq.gz'|sed 's/.fastq.gz//')
ALL_FASTQs = $(addprefix $(FASTQ_DIR)/,$(addsuffix _1.fastq.gz,$(PAIRED)) $(addsuffix _2.fastq.gz,$(PAIRED)) $(addsuffix .fastq.gz,$(UNPAIRED)))
FASTQs = $(XRRs)
RG=\\tSM:$(QUERY_INPUT)\\tPL:$(shell $(QUERY_EXE) --other_info=Platform)
SAMFILE_TEMPDIR=../SAMFILES/$(SCI_NAME)/$(QUERY_INPUT)
TEMP_SAMFILES=$(addprefix $(SAMFILE_TEMPDIR)/,$(addsuffix _paired.sam,$(PAIRED)) $(addsuffix _unpaired.sam,$(UNPAIRED)))
SAMFILES=$(addprefix $(SAMFILE_TEMPDIR)/,$(addsuffix _sorted.bam, $(addsuffix _paired,$(PAIRED)) $(addsuffix _unpaired,$(UNPAIRED))))
BAMFILE_DIR=../BAMFILES/$(SCI_NAME)
BAMFILE=$(BAMFILE_DIR)/$(QUERY_INPUT)
KNOWN_SNP_NAME=$(shell $(QUERY_EXE) --other_info=known_SNPs_genome_name)
KNOWNVARIANTS_DIR=../knownSNPs/$(KNOWN_SNP_NAME)
KNOWNVARIANTS=$(shell ls $(KNOWNVARIANTS_DIR)/*.vcf)
VCF_DIR=../VCFs/$(SCI_NAME)
VCFs=$(VCF_DIR)/$(QUERY_INPUT)
DIRs_to_clean=$(FASTQ_DIR)/* $(SAMFILE_TEMPDIR)/* $(BAMFILE_DIR)/*

#.INTERMEDIATE : $(SAMFILES) $(BAMFILE)_dedup_reads.bam $(BAMFILE)_target_intervals.list $(BAMFILE)_mpileup.vcf $(BAMFILE)_recal_data.grp $(BAMFILE)_recal_reads.bam
## NOTE 	: FASTQ files are obtained through the makefile 'fastqs.mk', which must be executed before this makefile executes
## all		: execute all functions in this makefile
.PHONY : all
all: $(VCFs)_raw_variants.vcf

## HaplotypeCaller	: call first round of haplotypes on analysis-ready reads
.PHONY : HaplotypeCaller
HaplotypeCaller : $(VCFs)_raw_variants.vcf
$(VCFs)_raw_variants.vcf : $(BAMFILE)_recal_reads.bam
\tmkdir -p $(VCF_DIR) && ../jdk1.8.0_91/bin/java -jar -Djava.io.tmpdir=`pwd`/tmp ../GenomeAnalysisTK_3.5.jar -T HaplotypeCaller -nct $(NPROC) -variant_index_type LINEAR -variant_index_parameter 128000 -R $(REF_GENOME_PATH) -I $(BAMFILE)_recal_reads.bam -o $(VCFs)_raw_variants.vcf -ERC GVCF && rm $(BAMFILE)_recal_reads.ba*

## recalibrate_scores	: use the bqsr recalibration to recalibrate sequence data
.PHONY : recalibrate_scores
.SECONDARY : $(BAMFILE)_recal_reads.bam
recalibrate_scores : $(BAMFILE)_recal_reads.bam $(BAMFILE)_realigned_reads.bam
$(BAMFILE)_recal_reads.bam : $(BAMFILE)_realigned_reads.bam $(BAMFILE)_recal_data.grp
\t$(GATK_EXE) -T PrintReads -nct $(NPROC) -R $(REF_GENOME_PATH) -I $(BAMFILE)_realigned_reads.bam -BQSR $(BAMFILE)_recal_data.grp -o $(BAMFILE)_recal_reads.bam && rm $(BAMFILE)_realigned_reads.ba*

$(BAMFILE)_recal_data.grp : $(VCF_FILES) $(BAMFILE)_realigned_reads.bam
\t$(GATK_EXE) -T BaseRecalibrator -nct $(NPROC) -R $(REF_GENOME_PATH) -I $(BAMFILE)_realigned_reads.bam $(addprefix --knownSites , $(KNOWNVARIANTS)) -o $(BAMFILE)_recal_data.grp && $(GATK_EXE) -T BaseRecalibrator -nct $(NPROC) -R $(REF_GENOME_PATH) -I $(BAMFILE)_realigned_reads.bam $(addprefix --knownSites , $(KNOWNVARIANTS)) -BQSR $(BAMFILE)_recal_data.grp -o $(BAMFILE)_post_recal_data.grp && $(GATK_EXE) -T AnalyzeCovariates -R $(REF_GENOME_PATH) -before $(BAMFILE)_recal_data.grp -after $(BAMFILE)_post_recal_data.grp -plots $(BAMFILE)_BQSR.pdf


## realign_indels	: realign raw sequencing reads around indels using GATK's RealignerTargetCreator and IndelRealiner
.PHONY : realign_indels
.SECONDARY : $(BAMFILE)_realigned_reads.bam
realign_indels : $(BAMFILE)_realigned_reads.bam
$(BAMFILE)_realigned_reads.bam : $(BAMFILE)_target_intervals.list $(BAMFILE)_dedup_reads.bam
\t$(GATK_EXE) -T IndelRealigner -R $(REF_GENOME_PATH) -I $(BAMFILE)_dedup_reads.bam -targetIntervals $(BAMFILE)_target_intervals.list -o $(BAMFILE)_realigned_reads.bam && rm $(BAMFILE)_dedup_reads.ba*

$(BAMFILE)_target_intervals.list : $(BAMFILE)_dedup_reads.bam
\t$(GATK_EXE) -T RealignerTargetCreator -nt $(NPROC) -R $(REF_GENOME_PATH) -I $(BAMFILE)_dedup_reads.bam -o $(BAMFILE)_target_intervals.list

## remove_dups	: remove duplicate reads from SAMfiles, sort, and convert to a single merged BAM file
.PHONY : remove_dups
.SECONDARY : $(BAMFILE)_dedup_reads.bam
remove_dups : $(BAMFILE)_dedup_reads.bam
$(BAMFILE)_dedup_reads.bam : $(SAMFILES)
\tmkdir -p $(BAMFILE_DIR)/metrics; $(PICARD_EXE) MarkDuplicates $(addprefix I=, $(SAMFILES)) OUTPUT=$(BAMFILE)_dedup_reads.bam M=$(BAMFILE_DIR)/metrics/$(QUERY_INPUT).metrics && $(PICARD_EXE) BuildBamIndex I=$(BAMFILE)_dedup_reads.bam && rm -r $(SAMFILE_TEMPDIR)/* $(FASTQ_DIR)/*

#mark samfiles and fastq files as secondary:
## map_fastqs	: map fastq files to appropriate genome
.PHONY : map_fastqs sort_SAM
.SECONDARY : $(TEMP_SAMFILES)
sort_SAM : $(SAMFILES)
%_paired_sorted.bam : $(TEMP_SAMFILES)
\tmkdir -p ../tmp ; $(SORTSAM_EXE) I=$*_paired.sam O=../tmp/$(notdir $*)_paired_sorted.bam SORT_ORDER=coordinate && rm $*_paired.sam && mv ../tmp/$(notdir $*)_paired_sorted.bam $(SAMFILE_TEMPDIR)

%_unpaired_sorted.bam : $(TEMP_SAMFILES)
\tmkdir -p ../tmp ; $(SORTSAM_EXE) I=$*_unpaired.sam O=../tmp/$(notdir $*)_unpaired_sorted.bam SORT_ORDER=coordinate && rm $*_unpaired.sam && mv ../tmp/$(notdir $*)_unpaired_sorted.bam $(SAMFILE_TEMPDIR)

map_fastqs : $(TEMP_SAMFILES)
%_paired.sam : $(ALL_FASTQs)
\tmkdir -p $(dir $*) ; $(BWA_EXE) mem -t $(NPROC_BWA) -M -R '@RG\\tID:$(notdir $*)$(RG)' $(REF_GENOME_PATH) $(FASTQ_DIR)/$(notdir $*)_1.fastq.gz $(FASTQ_DIR)/$(notdir $*)_2.fastq.gz > $*_paired.sam 2> $(notdir $*)_bwa.stderr
%_unpaired.sam : $(ALL_FASTQs)
\tmkdir -p $(dir $*) ; $(BWA_EXE) mem -t $(NPROC_BWA) -M -R '@RG\\tID:$(notdir $*)$(RG)' $(REF_GENOME_PATH) $(FASTQ_DIR)/$(notdir $*).fastq.gz > $*_unpaired.sam 2> $(notdir $*)_bwa.stderr

## print_fastqs	: print full path fastq files to be written
.PHONY : print_fastqs
print_fastqs :
\t@echo $(ALL_FASTQs)
## give_org	: give the scientific name of organism attached to this biosample
give_org :
\t@echo $(SCI_NAME)

## print_SRRs	: print the SRR files that will be downloaded
.PHONY : print_SRRs
print_SRRs :
\t@printf '\\nPAIRED = $(PAIRED)\\nUNPAIRED = $(UNPAIRED)\\n'

## print_known_SNPs	: print the known SNP vcf files obtained from different sources
.PHONY : print_known_SNPs
print_known_SNPs :
\t@echo $(KNOWN_SNP_NAME)
## clean		: clean up temporary files
.PHONY : clean
clean :
\trm -r $(DIRs_to_clean)

## help		: print help statement
.PHONY : help
help : Makefile
\t@sed -n 's/^##//p' $<
'''.format(input))

    with open("../PBS/RheMac_1/"+input+".pbs",'w') as wfile3:
        wfile3.write('''#PBS -N {0}
#PBS -q long8gb
#PBS -S /bin/bash
#PBS -m ae
#PBS -M alst5940@colorado.edu
#PBS -e localhost:/scratch/Users/alst5940/Primate_data_mining/PBS/eo_files/{0}.err
#PBS -o localhost:/scratch/Users/alst5940/Primate_data_mining/PBS/eo_files/{0}.out
#PBS -l nodes=1:ppn=64
#PBS -l mem=500gb
#PBS -l walltime=239:59:59
cd /scratch/Users/alst5940/Primate_data_mining/makefiles
make -j -f {0}_Makefile.mk all
'''.format(input))

if __name__ == "__main__":
    input = sys.argv[1:][0]
    write_makefiles_and_pbs(input)
#    os.system('qsub ../PBS/GAGP_jobs/{0}.pbs'.format(input))


