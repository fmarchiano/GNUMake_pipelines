ifndef CONFIG_INC


############################################################
########## OPTIONS MOST LIKELY NEED TO BE CUSTOMIZED ##########
############################################################


# reference genome
# e.g.: hg19, b37, hg19_ionref, b37_hbv_hcv etc
REF ?= b37

# possible values: ILLUMINA, IONTORRENT
SEQ_PLATFORM ?= ILLUMINA

# possible values: NONE, BAITS, PCR, RNA, CHIP
CAPTURE_METHOD ?= NONE

# e.g: NONE, AGILENT_CLINICAL_EXOME, CCP, AGILENT_ALLEXON_MOUSE, HCC, WXS
# For WGS, RNA-seq and ChiP-seq, this should be NONE.
PANEL ?= NONE

# Single-end or paired-end, set to false if single-end
PAIRED_END ?= true

# for pair-end RNA-seq, [NONE|FIRST_READ_TRANSCRIPTION_STRAND|SECOND_READ_TRANSCRIPTION_STRAND]
STRAND_SPECIFICITY ?= NONE

# possible values: SOMATIC, GERMLINE
ANALYSIS_TYPE ?= SOMATIC

# if you anticipate a complex project requiring analyzing subsets of samples
# this is useful for distinguishing them
PROJECT_PREFIX ?=

INCLUDE_CHR_Y ?= true

# This parameter determines which aligned BAM files get put into the bam/ directory
# This is useful when multiple aligners are used and you want to specify the
# set that goes into bam/, since most downstream tools get the data from bam/
ifeq ($(strip $(CAPTURE_METHOD)),RNA)
PRIMARY_ALIGNER ?= star
else
PRIMARY_ALIGNER ?= bwamem
endif

# Specify the sample sheets
SAMPLE_SET_FILE ?= sample_sets.txt
SAMPLE_FILE ?= samples.txt
SAMPLE_SPLIT_FILE ?= samples.split.txt

# Normal samples to use for pooled_norm and panel_of_normals: 
# particularly useful for when samples are being added to specify a fixed set of normals to use.
# Otherwise, if one single new normal sample is added, everything re-runs.
SAMPLE_PON_FILE ?= samples.pon.txt
SAMPLE_POOLEDNORM_FILE ?= samples.poolednorm.txt

MERGE_SAMPLE_FILE ?= merge_samples.txt


ifeq ($(HPC),ubelix)
  GENCODE_VER ?= v48
else
  GENCODE_VER ?= v33
endif

# is it a tumor-only dataset
TUMOR_ONLY ?= false

############################################################
########### END MOST FREQUENTLY CUSTOMIZED PARAMS ##########
############################################################

# this grabs the parameters for the genome and target panels
include usb-modules-v2/genome_inc/$(REF).inc


#############################################
############### Pre-alignment ###############
#############################################

## These parameters are used when you align FASTQs to BAMs
## If MERGE_SPLIT_FASTQ=true, this happens before FASTQ_FILTER
## If set, <FASTQ_FILTER> is run before alignment

# If multiple sets of fastqs per sample, should they be merged before alignment
# For transcriptome, true. For genomic, false.
# There are performance advantages to setting this to false.
ifeq ($(CAPTURE_METHOD),RNA)
MERGE_SPLIT_FASTQ ?= true 
endif
ifeq ($(wildcard $(SAMPLE_SPLIT_FILE)),)
MERGE_SPLIT_FASTQ ?= true
endif
MERGE_SPLIT_FASTQ ?= false

# Pre-processing the FASTQs, trimming bad cycles or adaptors
# [<NULL>|trim|trimgalore|cutadapt]
# RNA-seq needs the adaptors trimmed off so need either trimgalore or cutadapt
# trim: uniforming trim reads down to TRIM_LENGTH (see below)
# trimgalore: trims adaptors, and optionally clip reads to CLIP_FASTQ_R1/2 base pairs
# cutadapt: trims adaptors, but more flexible options than trimgalore
ifeq ($(CAPTURE_METHOD),RNA)
FASTQ_FILTER ?= cutadapt
endif
FASTQ_FILTER ?=

# Uniforming trim to this length if FASTQ_FILTER=trim
TRIM_LENGTH ?= 150

# setting the --clip_R1 and --clip_R2 parameters if FASTQ_FILTER=trimgalore
CLIP_FASTQ_R1 ?=
CLIP_FASTQ_R2 ?=
TRIM_GALORE_OPTS ?=


# options for FASTQ_FILTER=cutadapt
# ADAPTOR_POS =  [3prime|5prime]
# NEXTSEQ_TRIM = true may help for sequencing done with NextSeq
ADAPTOR_POS ?= 3prime
ADAPTOR_SEQ ?= AGATCGGAAGAGC
NEXTSEQ_TRIM ?= false

ifeq ($(ADAPTOR_POS),3prime)
CUTADAPT_OPTS ?= -a $(ADAPTOR_SEQ) $(if $(findstring true,$(PAIRED_END)),-A $(ADAPTOR_SEQ))
endif
ifeq ($(ADAPTOR_POS),5prime)
CUTADAPT_OPTS ?= -g $(ADAPTOR_SEQ) $(if $(findstring true,$(PAIRED_END)),-G $(ADAPTOR_SEQ))
endif
ifeq ($(NEXTSEQ_TRIM),true)
CUTADAPT_OPTS += --nextseq-trim=20
endif
CUTADAPT_OPTS += -q 20,20 --trim-n -m 35 --max-n=3


#############################################
############### Post-alignment ###############
#############################################

# After alignment with the favourite aligner, 
# the BAM files are post-processed in the following order
# For RNA-seq, if mutation calls will eventually be needed, 
# then post-process, otherwise just dont, it takes a long time and creates hard-clipped
# reads that may be problematic. Use true with caution
POST_PROCESS_RNA_BAM ?= false 

# this runs Picard ReorderSam
BAM_REORDER ?= true 

# this runs Picard SortSam
BAM_SORT ?= true 

# this runs Picard AddOrReplaceReadGroups
BAM_FIX_RG ?= false 

# this filters out reads with BAM_FILTER_FLAGS
# default 768 means "not primary alignment" and "read fails platform/vendor quality checks"
BAM_FILTER ?= true 
BAM_FILTER_FLAGS ?= 768

# this removes inter-chromosomal chimeric paired reads
BAM_INTRACHR ?= false 

# this runs GATK SplitNCigarReads, default to false except for RNA-seq
ifeq ($(CAPTURE_METHOD),RNA)
ifeq ($(POST_PROCESS_RNA_BAM),true)
BAM_SPLIT_N_TRIM ?= true 
endif
endif
BAM_SPLIT_N_TRIM ?= false 

# this splits bams into individual chroms betfore BAM_DUP_TYPE, 
# BAM_REALN and BAM_RECAL before remerging the BAMs.
# defaults to true except for RNA-seq, amplicon-based capture and CHIP
# setting this to false means BAM_DUP_TYPE, BAM_REALN and BAM_RECAL are done on a genome-wide level
ifeq ($(CAPTURE_METHOD),RNA)
SPLIT_CHR ?= false
endif
ifeq ($(CAPTURE_METHOD),PCR)
SPLIT_CHR ?= false
endif
ifeq ($(CAPTURE_METHOD),CHIP)
SPLIT_CHR ?= false
endif
SPLIT_CHR ?= true


# This runs deduplication, GATK IndelRealigner and BaseRecalibrator
# deduplication: [none|markdup|rmdup]: markdup is Picard MarkDuplicates, rmdup is samtools rmdup
# defaults to markdup and true except for RNA-seq, amplicon-based capture and CHIP
ifeq ($(CAPTURE_METHOD),RNA)
ifeq ($(POST_PROCESS_RNA_BAM),true)
BAM_DUP_TYPE ?= markdup
endif
BAM_DUP_TYPE ?= none
BAM_REALN ?= false
BAM_RECAL ?= false
endif
ifeq ($(CAPTURE_METHOD),PCR)
BAM_DUP_TYPE ?= none
BAM_REALN ?= false
BAM_RECAL ?= false
endif
ifeq ($(CAPTURE_METHOD),CHIP)
BAM_DUP_TYPE ?= none
BAM_REALN ?= false
BAM_RECAL ?= false
endif
BAM_DUP_TYPE ?= markdup 
BAM_REALN ?= true 
BAM_RECAL ?= true 
CHOOSE_CHR_FOR_RECAL ?= true

#########################################################
############### Additional BAM processing ###############
#########################################################

# if you have a set of aligned BAM files you want to put
# through (any step) of the post-alignment sequence above
# set this to true, and check and set the post-alignment
# parameters in the previous section accordingly, then
# put the bam files in the unprocessed_bam direcory
BAM_REPROCESS ?= false

BAM_PHRED64 ?= false # obsolete parameter
MERGE_SPLIT_BAMS ?= false  # internal parameter

# if you need to downsample BAM files, for example to
# make a pooled normal BAM file, then you set this.
# see the samtools view -s option
SAMTOOLS_DOWNSAMPLE_FACTOR ?= 1001.026

# For QC, genotyping is usually done. For performance reasons,
# we usually use chromosome 2 for genotyping (a large chromosome
# that is not usually gained/lost in cancer. For very small panel,
# this is not good. This is also specified in some genome_inc panel
# config files. Values in genome_inc configs overide the values here.
# GENOTYPE_CHR can be "allchr" or the chromosome of your choice.
ifeq ($(CAPTURE_METHOD),PCR)
GENOTYPE_CHR ?= allchr
else
ifeq ($(CAPTURE_METHOD),CHIP)
GENOTYPE_CHR ?= allchr
else
ifeq ($(REF),b37)
GENOTYPE_CHR ?= 2
else
ifeq ($(REF),b37_GRCm38)
GENOTYPE_CHR ?= 2
else
ifeq ($(REF),b37_hbv_hcv)
GENOTYPE_CHR ?= 2
else
ifeq ($(REF),b37_hbv_hcv_GRCm38)
GENOTYPE_CHR ?= 2
else
ifeq ($(REF),hg38)
GENOTYPE_CHR ?= chr2
else
ifeq ($(REF),hg38_gallchrus6a)
GENOTYPE_CHR ?= chr2
else
ifeq ($(REF),hg38_GRCm38)
GENOTYPE_CHR ?= chr2
else
GENOTYPE_CHR ?= allchr
endif
endif
endif
endif
endif
endif
endif
endif
endif

##################################
########## Aligners ##############
##################################

HISAT2_OPTS = --dta 

################################################
###### Somatic variant callers #################
################################################

MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_GT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.5
MUTECT_OTHER_OPTS ?= 
MUTECT2_OTHER_OPTS ?= 

CAVEMAN_NUM_CORES ?= 12
CAVEMAN_NORMAL_CONTAMINATION ?= 0.1
CAVEMAN_USE_CN ?= true
CAVEMAN_CF_THRESHOLD ?= 0.3
CAVEMAN_TUM_CN_DEFAULT ?= 2
CAVEMAN_NORM_CN_DEFAULT ?= 2

USE_BAM_CLIPOVERLAP ?= true

MANTA_NUM_CORES ?= 12
MANTA_MEM ?= 24

SVABA_NUM_CORES ?= 8

DEEPSOMATIC_NUM_CORES ?= 12
ifeq ($(CAPTURE_METHOD),BAITS)
  DEEPSOMATIC_MODEL ?= $(if $(findstring true,$(TUMOR_ONLY)),WES_TUMOR_ONLY,WES)
else
  DEEPSOMATIC_MODEL ?= $(if $(findstring true,$(TUMOR_ONLY)),WGS_TUMOR_ONLY,WGS)
endif

################################################
###### Somatic variant filters #################
################################################

# after somatic variant callers, the calls go through
# a series of filters. If you don't want the filters
# set this to false
FILTER_VARIANTS ?= true

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
MIN_NORMAL_DEPTH ?= 5
MIN_TUMOR_DEPTH ?= 10
MIN_TUMOR_AD ?= 3
MIN_TN_AD_RATIO ?= 5.0
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
ifeq ($(TUMOR_ONLY),false)
MIN_NORMAL_DEPTH ?= 10
MIN_TUMOR_DEPTH ?= 10
MIN_TUMOR_AD ?= 8
MIN_TN_AD_RATIO ?= 10.0
MIN_AF ?= 0.01
endif
ifeq ($(TUMOR_ONLY),true)
MIN_TUMOR_DEPTH ?= 20
MIN_TUMOR_AD ?= 8
MIN_AF ?= 0.1
MAX_POP ?= 0.005
MAX_HOM ?= 4
endif
endif

PIPEIT_INTERMEDIATE_FILES ?= false

MIN_MQ ?= 10
MIN_AF_SNP ?= 0.01
MIN_AF_INDEL ?= 0.02

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
MUT_CALLER ?= tvc
endif

ifeq ($(findstring tvc,$(MUT_CALLER)),tvc)
USE_NFT ?= false
ANN_NFT ?= true
PON_VCF ?= tvc/pon.tvc.vcf.gz
endif

ifeq ($(findstring pipeit,$(MUT_CALLER)),pipeit)
USE_NFT ?= false
ANN_NFT ?= false
endif

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
USE_NFT ?= true
ANN_NFT ?= false
PON_VCF ?= mutect2/pon.mutect2.vcf.gz
endif

FILTER_OXOG ?= false
FILTER_FFPE ?= false

PON_MIN_SAMPLES ?= 2

VCF_PASS_MAX_FILTERS ?= 2

#################################################
###### Germline variant callers #################
#################################################

GATK_INTERVALS_COUNT = 25
GATK_HARD_FILTER_SNPS ?= false
GATK_HARD_FILTER_INDELS ?= false

VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL_SNPS = 99.5
VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL_INDELS = 99.0

## These are used if GATK_HARD_FILTER_SNPS ##
HAPCALL_SNP_MQ_THRESHOLD ?= 40.0
HAPCALL_SNP_QD_THRESHOLD ?= 2.0
HAPCALL_SNP_FS_THRESHOLD ?= 60.0
HAPCALL_SNP_HAP_SCORE_THRESHOLD ?= 13.0
HAPCALL_SNP_MQ_RANKSUM_THRESHOLD ?= -12.5
HAPCALL_SNP_READPOS_RANKSUM_THRESHOLD ?= -8.0

## These are used if GATK_HARD_FILTER_INDELS ##
HAPCALL_INDEL_INBREED_COEFF_THRESHOLD ?= -0.8
HAPCALL_INDEL_QD_THRESHOLD ?= 2.0
HAPCALL_INDEL_FS_THRESHOLD ?= 200.0
HAPCALL_INDEL_HAP_SCORE_THRESHOLD ?= 13.0
HAPCALL_INDEL_MQ_RANKSUM_THRESHOLD ?= -20.0
HAPCALL_INDEL_READPOS_RANKSUM_THRESHOLD ?= -8.0



##########################################
########## Variant Annotation ############
##########################################

ANNOTATE_VARIANTS ?= true
USE_SUFAM ?= true
ANN_FACETS ?= false

MUTATION_SUMMARY_FORMAT ?= EXCEL
INCLUDE_LNCRNA_IN_SUMMARY ?= false

SNP_EFF_FLAGS ?= -canon -noNextProt -no-intergenic -no-downstream -no-upstream -no-utr -noStats -noInteraction

EXAC_INFO_FIELDS ?= ExACnontcga_AC,ExACnontcga_AF
EXACNONPSYCH_INFO_FIELDS ?= ALL

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
# CALLER_PREFIX ?= mutect strelka_indels
CALLER_PREFIX ?= mutect2 strelka2_indels strelka2_snvs
endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
CALLER_PREFIX ?= tvc_snps tvc_indels
endif

SKIP_3D_HOTSPOTS ?= yes
##################################
########## Copy number ###########
##################################

FACETS_GATK_VARIANTS ?= false
FACETS_MINGC ?= 0
FACETS_MAXGC ?= 1

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
 ifeq ($(findstring NONE,$(PANEL)),NONE)
  FACETS_SNP_PILEUP_MIN_DEPTH ?= 5
  FACETS_SNP_PILEUP_MAX_DEPTH ?= 500
  FACETS_SNP_PILEUP_PSEUDO_SNPS ?= 1000
  FACETS_WINDOW_SIZE ?= 200
  FACETS_PRE_CVAL ?= 80
  FACETS_CVAL ?= 200
  FACETS_MAX_SEGS ?= 500
  FACETS_HETMARKFREQ_MINMARKS ?= 10000
 else
  FACETS_SNP_PILEUP_MIN_DEPTH ?= 25
  FACETS_SNP_PILEUP_MAX_DEPTH ?= 1000
  FACETS_SNP_PILEUP_PSEUDO_SNPS ?= 100
  FACETS_WINDOW_SIZE ?= 200
  FACETS_PRE_CVAL ?= 25
  FACETS_CVAL ?= 150
  FACETS_MAX_SEGS ?= 300
  FACETS_HETMARKFREQ_MINMARKS ?= 1000
 endif
endif


ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
FACETS_SNP_PILEUP_MIN_DEPTH ?= 50
FACETS_SNP_PILEUP_MAX_DEPTH ?= 5000
FACETS_SNP_PILEUP_PSEUDO_SNPS ?= 50
FACETS_WINDOW_SIZE ?= 10
FACETS_PRE_CVAL ?= 25
FACETS_CVAL ?= 75
FACETS_MAX_SEGS ?= 200
FACETS_HETMARKFREQ_MINMARKS ?= 500
endif

FACETS_SNP_PILEUP_MINMAPQ ?= 1
FACETS_SNP_PILEUP_MINBASEQ ?= 13


FACETS_MIN_NHET ?= 15
FACETS_SPECIAL_CASES ?=
FACETS_RUN_GENE_CN ?= false
FACETS_GENE_CN_OPTS ?=
FACETS_GENE_CN_PLOT_OPTS ?=
FACETS_HETMARKFREQ_THRESHOLD ?= 0.025

ifeq ($(CAPTURE_METHOD),PCR)
CBS_SEG_SD ?= 1.5
CBS_SEG_SMOOTH ?= 15
CBS_SEG_ALPHA ?= 0.0000000001
CBS_TRIM ?= 0.1
CBS_CLEN ?= 5
CBS_MIN_N_DEPTH ?= 100
endif

ifeq ($(CAPTURE_METHOD),RNA)
CBS_MIN_N_DEPTH ?= 30
CBS_MAX_N_DEPTH ?= 5000
CBS_MIN_T_DEPTH ?= 1
#CBS_OUTLIER_SD_SCALE ?= 1
#CBS_SEG_SD ?= 1
endif

CBS_SEG_SD ?= 2
CBS_SEG_SMOOTH ?= 10
CBS_SEG_ALPHA ?= 0.000001
CBS_TRIM ?= 0.025
CBS_CLEN ?= 10
CBS_OUTLIER_SD_SCALE ?= 2.5

CBS_EXCL_N_OUTLIER_PC ?= 0.05
CBS_MIN_N_DEPTH ?= 0

CBS_MULTIPARAM_SEGMENT ?= false
CBS_SEG_SDS ?= 1.5 2 2.5 3
CBS_SEG_SMOOTHS ?= 5 10
CBS_SEG_ALPHAS ?= 0.01 0.0001 0.000001 0.0000000001

VARSCAN_GENE_CN_OPTS = $(FACETS_GENE_CN_OPTS)

CNVKIT_CORRELATION ?= $(GROUP_DIR)/ref_nobackup/cnvkit_reference/cnvkit_correlation_tcga-lihc.tsv
CNVKIT_GENE_CN_OPTS ?= $(FACETS_GENE_CN_OPTS)


ifeq ($(CAPTURE_METHOD),PCR)
VARSCAN_CNV_MAX_SEG_SIZE = 1000
CBS_MPILEUP_OPTS = -q 1 -d 20000
endif

VARSCAN_CNV_MAX_SEG_SIZE ?= 100
CBS_MPILEUP_OPTS ?= -q 1

##################################
###### MUTATIONAL SIGNATURES #####
##################################
SIG_PROFILER_COSMIC_VERSION ?= 3.4

SIG_PROFILER_COSMIC_EXOME ?= 
ifeq ($(CAPTURE_METHOD),BAITS)
SIG_PROFILER_COSMIC_EXOME = True
endif

ifeq ($(REF),hg38)
 SIG_PROFILER_COSMIC_GENOME ?= GRCh38
else ifeq ($(REF),hg38_GRCm38)
 SIG_PROFILER_COSMIC_GENOME ?= GRCh38
else ifeq ($(REF),hg38_gallus6a)
 SIG_PROFILER_COSMIC_GENOME ?= GRCh38
else ifeq ($(REF),b37)
 SIG_PROFILER_COSMIC_GENOME ?= GRCh37
else ifeq ($(REF),b37_GRCm38)
 SIG_PROFILER_COSMIC_GENOME ?= GRCh37
else ifeq ($(REF),b37_hbv_hcv)
 SIG_PROFILER_COSMIC_GENOME ?= GRCh37
else ifeq ($(REF),b37_hbv_hcv_GRCm38)
 SIG_PROFILER_COSMIC_GENOME ?= GRCh37
else ifeq ($(REF),hg19_ionref)
 SIG_PROFILER_COSMIC_GENOME ?= GRCh37
else ifeq ($(REF),GRCm38)
 SIG_PROFILER_COSMIC_GENOME ?= mm10
endif

SIG_PROFILER_COSMIC_SIGNATURE_DB ?= 
SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS ?= 



##################################
############# GISTIC ############
##################################

GISTIC_THRESHOLD ?= 0.3
GISTIC_JS ?= 15

##################################
########## Gene expression #######
##################################

DESEQ_CONDITION ?= condition
DESEQ_REF_CONDITION ?= ref
DESEQ_ALT_CONDITION ?= alt
DESEQ_PHENO_FILE ?= pheno.txt


#################################
####### ChIP-seq ################
#################################

MOSAICS_FRAG_LEN ?= 200
MOSAICS_BIN_SIZE ?= 200
MOSAICS_NUM_CORES ?= 8
MOSAICS_PARALLEL ?= TRUE
MOSAICS_MAXGAP ?= 200
MOSAICS_MINSIZE ?= 50
MOSAICS_THRES ?= 10

MACS2_PVALUE ?= 0.01
MACS2_READLENGTH ?= 100
MACS2_NUM_CORES ?= 4

##################################
########## Fusions ###############
##################################

STAR_FUSION_MIN_JUNCTION_READS ?= 2
STAR_FUSION_MIN_SUM_FRAGS ?= 3
STAR_FUSION_MAX_PROMISCUITY ?= 3
STAR_FUSION_MIN_NOVEL_JUNCTION_SUPPORT ?= 3
STAR_FUSION_MIN_ALT_PCT_JUNC ?= 10

STAR_FUSION_MIN_JUNCTION_READS_NORMAL ?= 0
STAR_FUSION_MIN_SUM_FRAGS_NORMAL ?= 2
STAR_FUSION_MAX_PROMISCUITY_NORMAL ?= 5
STAR_FUSION_MIN_NOVEL_JUNCTION_SUPPORT_NORMAL ?= 2
STAR_FUSION_MIN_ALT_PCT_JUNC_NORMAL ?= 5

##################################
####### CLONALITY ################
##################################

ABSOLUTE_CALLER_PREFIX ?=

ABSOLUTE_STEP_3 ?= false
ABSOLUTE_MAKE_MUTS_OPT ?=
ABSOLUTE_MAKE_SEGS_OPT ?=

ABSOLUTE_SIGMA_P ?= 0
ABSOLUTE_MAX_SIGMA_H ?= 0.07
ABSOLUTE_MIN_PLOIDY ?= 0.9
ABSOLUTE_MAX_PLOIDY ?= 8
ABSOLUTE_DISEASE ?= NA
ABSOLUTE_PLATFORM ?= Illumina_WES
ABSOLUTE_MAX_SEG ?= 3500
ABSOLUTE_COPYNUMTYPE ?= total
ABSOLUTE_MAX_NEG_GENOME ?= 0.05
ABSOLUTE_MAX_NON_CLONAL ?= 0.2
ABSOLUTE_MIN_MUT_AF ?= 0

PYCLONE_PRIOR ?= parental_copy_number

################################
###### Miscellaneous ###########
###############################

DECONSTRUCTSIGS_NUMITER ?= 100
DECONSTRUCTSIGS_ASSOC_SIGS ?= 

PVACSEQ_ALGO ?= NetMHCpan
PVACSEQ_USE_LOCAL_IEDB ?= true

RSEM_NUM_CORES ?= 4
RSEM_CALC_PME ?=

ifeq ($(STRAND_SPECIFICITY),NONE)
RSEM_STRANDEDNESS ?= none
endif
ifeq ($(STRAND_SPECIFICITY),FIRST_READ_TRANSCRIPTION_STRAND)
RSEM_STRANDEDNESS ?= forward
endif
ifeq ($(STRAND_SPECIFICITY),SECOND_READ_TRANSCRIPTION_STRAND)
RSEM_STRANDEDNESS ?= reverse
endif

RSEM_OPTIONS ?= --alignments --no-bam-output -p $(RSEM_NUM_CORES) --strandedness $(RSEM_STRANDEDNESS) $(if $(RSEM_CALC_PME),--calc-pme)

STAR_OPTIONS ?= --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --alignSJDBoverhangMin 1 --outSAMprimaryFlag AllBestScore --outMultimapperOrder Random --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimSegmentReadGapMax 3 --chimOutType Junctions WithinBAM SoftClip --twopassMode Basic

endif
CONFIG_INC = true
