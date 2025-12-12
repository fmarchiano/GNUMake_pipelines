ifndef GRCM38_INC
REF = GRCm38
include usb-modules-v2/genome_inc/common.inc


REF_FASTA ?= $(GENOME_DIR)/GRCm38_68.resorted.fa
REF_DICT = $(REF_FASTA:.fa=.dict)

DBSNP := $(ANN_DIR)/mgp.v5.merged.snps_all.indels.dbSNP142_noGT.vcf

ifeq ($(findstring NONE,$(PANEL)),NONE)
else
include usb-modules-v2/genome_inc/$(REF)/$(PANEL).inc
endif

RESOURCE_REQ ?= high
#TARGETS_FILE_INTERVALS ?= $(ANN_DIR)/human_g1k_v37_decoy.bed
DBSNP_TARGETS_INTERVALS ?= $(DBSNP)

ifeq ($(findstring bed,$(TARGETS_FILE_INTERVALS)),bed)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE_INTERVALS) | cut -f1 | sort -V | uniq)
else
CHROMOSOMES ?= $(shell seq 1 19) X $(if $(findstring true,$(INCLUDE_CHR_Y)),Y) MT
endif

RECAL_CHR = 2

### ANNOTATION SETS ####

KNOWN_INDELS := $(ANN_DIR)/mgp.v5.merged.indels.dbSNP142.normed.noGT.vcf

GENCODE_GTF = $(ANN_DIR)/gencode.vM11.annotation_nochr.gtf
GENCODE_GENE_GTF = $(ANN_DIR)/gencode.vM11.annotation_nochr.genes.gtf
GENCODE_CODING_GENE_GTF = $(ANN_DIR)/gencode.vM11.annotation_nochr.protein_coding_genes.gtf
RIBOSOMAL_INTERVALS = $(ANN_DIR)/ribosomalRNA.interval_list
GENE_REF_FLAT = $(ANN_DIR)/refFlat.txt

#### TOOL SPECIFIC OPTIONS ####

ifeq ($(findstring scicore_ubuntu,$(HPC)),scicore_ubuntu)
  SNP_EFF_GENOME = GRCm38.99
else
  SNP_EFF_GENOME = GRCm38.86
endif


endif
GRCM38_INC = true
