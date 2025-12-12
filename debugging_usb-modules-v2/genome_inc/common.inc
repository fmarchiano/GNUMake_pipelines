#######################################################
##################### SCICORE #########################
#######################################################

ifeq ($(findstring scicore,$(HPC)),scicore)
  REF_DIR = $(GROUP_DIR)/ref_nobackup

  GENOME_DIR = $(REF_DIR)/genomes
  ANN_DIR = $(REF_DIR)/annotation_data/$(REF)
  BED_DIR = $(REF_DIR)/target_beds/$(REF)

  ifeq ($(CAPTURE_METHOD),NONE)
    FACETS_TARGETS_INTERVALS ?= $(DBSNP_COMMON)
  else
    FACETS_TARGETS_INTERVALS ?= $(DBSNP_TARGETS_INTERVALS)
  endif

  STAR_GENOME_DIR ?= $(REF_DIR)/genome_index/star_$(REF)
  STAR_CTAT_DIR ?= $(ANN_DIR)/STAR-FusionFilter

  RSEM_INDEX ?= $(REF_DIR)/genome_index/rsem_$(REF)/rsem_index

  HISAT2_REF ?= $(REF_DIR)/genome_index/hisat2_$(REF)/genome_snp_tran

  KALLISTO_INDEX ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).transcripts.shortname.k31.idx
  KALLISTO_GTF ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).primary_assembly.annotation.gtf
  KALLISTO_CHR ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).primary_assembly.annotation.chromosomes.len

  PIZZLY_GTF_CACHE ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).primary_assembly.annotation.gtf.pizzly.cache
  PIZZLY_FASTA ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).transcripts.shortname.fa.gz
endif


#######################################################
#################### HUMANITAS ########################
#######################################################


ifeq ($(HPC),humanitas)
  REF_DIR = $(GROUP_DIR)/ref

  GENOME_DIR = $(REF_DIR)/genomes/$(REF)
  ANN_DIR = $(GENOME_DIR)/annotation
  BED_DIR = $(REF_DIR)/target_beds/$(REF)

  ifeq ($(CAPTURE_METHOD),NONE)
    FACETS_TARGETS_INTERVALS ?= $(DBSNP_COMMON)
  else
    FACETS_TARGETS_INTERVALS ?= $(DBSNP_TARGETS_INTERVALS)
  endif

  STAR_GENOME_DIR ?= $(GENOME_DIR)/index/star_$(REF)
  STAR_CTAT_DIR ?=

  RSEM_INDEX ?= $(GENOME_DIR)/index/rsem_$(REF)/rsem_index

  HISAT2_REF ?= $(GENOME_DIR)/index/hisat2_$(REF)/genome_snp_tran

  KALLISTO_INDEX ?= $(GENOME_DIR)/index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).transcripts.shortname.k31.idx
  KALLISTO_GTF ?= $(GENOME_DIR)/index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).primary_assembly.annotation.gtf
  KALLISTO_CHR ?= $(GENOME_DIR)/index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).primary_assembly.annotation.chromosomes.len

  PIZZLY_GTF_CACHE ?= $(GENOME_DIR)/index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).primary_assembly.annotation.gtf.pizzly.cache
  PIZZLY_FASTA ?= $(GENOME_DIR)/index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).transcripts.shortname.fa.gz
endif

#######################################################
##################### UBELIX ##########################
#######################################################

ifeq ($(HPC),ubelix)
  REF_DIR = $(GROUP_DIR)/ref

  GENOME_DIR = $(REF_DIR)/genomes/$(REF)
  ANN_DIR = $(REF_DIR)/anno/$(REF)
  BED_DIR = $(REF_DIR)/target_beds/$(REF)

  ifeq ($(CAPTURE_METHOD),NONE)
    FACETS_TARGETS_INTERVALS ?= $(DBSNP_COMMON)
  else
    FACETS_TARGETS_INTERVALS ?= $(DBSNP_TARGETS_INTERVALS)
  endif

  STAR_GENOME_DIR ?= $(REF_DIR)/genome_index/star_$(REF)/STAR_2.7.11b
  STAR_CTAT_DIR ?= $(ANN_DIR)/STAR-FusionFilter

  RSEM_INDEX ?= $(REF_DIR)/genome_index/rsem_$(REF)/rsem_index

  HISAT2_REF ?= $(REF_DIR)/genome_index/hisat2_$(REF)/genome_snp_tran

  KALLISTO_INDEX ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).transcripts.shortname.k31.idx
  KALLISTO_GTF ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).primary_assembly.annotation.gtf
  KALLISTO_CHR ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).primary_assembly.annotation.chromosomes.len

  PIZZLY_GTF_CACHE ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).primary_assembly.annotation.gtf.pizzly.cache
  PIZZLY_FASTA ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).transcripts.shortname.fa.gz
endif
