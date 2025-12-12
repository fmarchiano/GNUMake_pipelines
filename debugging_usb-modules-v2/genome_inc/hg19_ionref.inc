## This is the reference for the Ion Torrent-version of hg19
## Only use this with Ion Torrent data that have been aligned 
## by the Torrent Server

ifndef HG19_IONREF_INC
REF = hg19_ionref
include usb-modules-v2/genome_inc/common.inc

REF_FASTA ?= $(GENOME_DIR)/hg19_ionref.fasta
REF_DICT = $(REF_FASTA:.fasta=.dict)

DBSNP ?= $(ANN_DIR)/dbsnp_138.hg19.vcf

ifeq ($(findstring NONE,$(PANEL)),NONE)
else
include usb-modules-v2/genome_inc/$(REF)/$(PANEL).inc
endif

RESOURCE_REQ ?= medium
TARGETS_FILE_INTERVALS ?= $(REF_FASTA:.fasta=.bed)
DBSNP_TARGETS_INTERVALS ?= $(DBSNP)
TARGETS_FILE_INTERVALS_POOLS_CNA ?= $(TARGETS_FILE_INTERVALS_POOLS)

ifeq ($(findstring bed,$(TARGETS_FILE_INTERVALS)),bed)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE_INTERVALS) | cut -f1 | sort -V | uniq)
else
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE_INTERVALS) | cut -f1 | sort -V| uniq)
endif


### ANNOTATION SETS ####

ANNOVAR_HUMANDB ?= $(REF_DIR)/annovar_humandb.hg19.2018-04-16
CLINVAR ?= $(ANN_DIR)/clinvar.vcf
EXAC_NONTCGA ?= $(ANN_DIR)/ExAC_nonTCGA.r0.3.1.sites.vep.edited.vcf
COSMIC ?= $(ANN_DIR)/CosmicCodingMuts.sorted.vcf

#HAPLOTYPE_INSUF_BED = $(ANN_DIR)/haploinsufficient.genes.bed
#CANCER_GENE_CENSUS_BED = $(ANN_DIR)/cancer_gene_census.bed
#KANDOTH_BED = $(ANN_DIR)/kandoth127smg.bed
#LAWRENCE_BED = $(ANN_DIR)/lawrence_cancer5000s.bed
#DGD_BED = $(ANN_DIR)/dgd_Hsa_all_v71.bed
GENE_SETS_GMT = $(ANN_DIR)/gene_ann_20200201.gmt
GENE_SETS_LIST = $(ANN_DIR)/gene_ann_20200201.list

CANCER_HOTSPOT_VCF = $(ANN_DIR)/hotspots_changv2_gao_nc.vcf

CENTROMERE_TABLE = $(ANN_DIR)/centromeres_hg19.txt

### TOOL-SPECIFIC REF ###

SNP_EFF_GENOME = GRCh37.75
DB_NSFP ?= $(REF_DIR)/dbNSFP/GRCh37/dbNSFP2.9.txt.gz

endif
HG19_IONREF_INC = true

