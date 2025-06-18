###############################################################
####### DEFINE PATHS/NAMES FOR SCRIPTS, BINARIES, MODULES #####
###############################################################

SHELL ?= /bin/bash
UNZIP ?= /usr/bin/unzip
GZIP ?= gzip -f
GUNZIP ?= gunzip -c

SINGULARITY_EXEC = singularity exec
R_LIBS ?= $(PWD)/R_lib/Rocker_4.2.2
R4_LIBS ?= $(R_LIBS)
R ?= $(SINGULARITY_EXEC) -B $(R_LIBS):/mnt/Rlib R_RStudio/Rocker_rstudio_4.2.2.sif R
R4 ?= $(R)
RSCRIPT ?= $(SINGULARITY_EXEC) -B $(R_LIBS):/mnt/Rlib R_RStudio/Rocker_rstudio_4.2.2.sif Rscript
R4SCRIPT ?= $(RSCRIPT)

FACETS_SNP_PILEUP = htstools/snp-pileup
FACETS = $(R4SCRIPT) facets_pipeline/R_scripts/facetsRun.R
FACETS_HETMARKFREQ = $(R4SCRIPT) facets_pipeline/R_scripts/facetsHetMarkFreq.R



##################################
########## facets ###########
##################################

FACETS_GATK_VARIANTS ?= false
FACETS_MINGC ?= 0
FACETS_MAXGC ?= 1


FACETS_SNP_PILEUP_MIN_DEPTH ?= 5
FACETS_SNP_PILEUP_MAX_DEPTH ?= 500
FACETS_SNP_PILEUP_PSEUDO_SNPS ?= 1000
FACETS_WINDOW_SIZE ?= 200
FACETS_PRE_CVAL ?= 80
FACETS_CVAL ?= 200
FACETS_MAX_SEGS ?= 500
FACETS_HETMARKFREQ_MINMARKS ?= 10000

FACETS_SNP_PILEUP_MINMAPQ ?= 1
FACETS_SNP_PILEUP_MINBASEQ ?= 13

FACETS_MIN_NHET ?= 15
FACETS_SPECIAL_CASES ?=
FACETS_RUN_GENE_CN ?= false
FACETS_GENE_CN_OPTS ?=
FACETS_GENE_CN_PLOT_OPTS ?=
FACETS_HETMARKFREQ_THRESHOLD ?= 0.025


###############################################################
############################## REF ############################
###############################################################

ANN_DIR = ref/annotation_data/hg38
DBSNP_TARGETS_INTERVALS ?= $(ANN_DIR)/dbsnp_151.hg38.vcf.gz
DBSNP_COMMON ?= $(ANN_DIR)/dbsnp_151.hg38.common.vcf.gz

###############################################################
#################### MISCELLANEOUS ############################
###############################################################

# define $(,) and $( ) for escaping commas and spaces
, := ,
# $( ) variable not working in make v4.3 (empty variable name error for '$(space) :=':
# space :=
# space +=
# $(space) :=
# $(space) +=
E :=
space := $E $E

NOW := $(shell date +"%F")
MKDIR = mkdir -p -m 775

<< = $(word 2,$^)
<<< = $(word 3,$^)
<<<< = $(word 4,$^)
<<<<< = $(word 5,$^)

define R_INIT
dir.create('$(@D)', showWarnings = F, recursive = T)
qw <- function(x) unlist(strsplit(x, "[[:space:]]+"))
endef