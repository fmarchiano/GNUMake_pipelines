ifndef HG38_GALLUS6A_INC
REF = hg38_gallus6a
include usb-modules-v2/genome_inc/common.inc

REF_FASTA ?= $(GENOME_DIR)/Homo_sapiens_assembly38__Gallus_gallus_000002315.5_GRCg6a.fasta
REF_DICT = $(REF_FASTA:.fasta=.dict)

DBSNP ?= $(ANN_DIR)/dbsnp_146.hg38.vcf.gz

ifeq ($(findstring NONE,$(PANEL)),NONE)
else
include usb-modules-v2/genome_inc/$(REF)/$(PANEL).inc
endif

RESOURCE_REQ ?= high
TARGETS_FILE_INTERVALS ?= $(REF_FASTA:.fasta=.bed)
DBSNP_TARGETS_INTERVALS ?= $(DBSNP)

ifeq ($(findstring bed,$(TARGETS_FILE_INTERVALS)),bed)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE_INTERVALS) | cut -f1 | sort -V | uniq)
else
CHROMOSOMES ?= $(foreach chr $(shell seq 1 22),chr$(chr)) chrX $(if $(findstring true,$(INCLUDE_CHR_Y)),chrY) chrM
endif

RECAL_CHR = chr2

# GENCODE_VER ?= v33
TARGETS_FILE_GENES ?= $(ANN_DIR)/Homo_sapiens_assembly38__Gallus_gallus_000002315.5_GRCg6a.gene.band.txt

### ANNOTATION SETS ####

CONTEST_POP_FILE ?= $(ANN_DIR)/hapmap_3.3_grch38_pop_stratified_af.vcf.gz

KNOWN_INDELS ?= $(ANN_DIR)/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
OMNI ?= $(ANN_DIR)/1000G_omni2.5.hg38.vcf.gz
HAPMAP ?= $(ANN_DIR)/hapmap_3.3.hg38.vcf.gz
R1000G ?= $(ANN_DIR)/1000G_phase1.snps.high_confidence.hg38.vcf.gz

CLINVAR = $(ANN_DIR)/clinvar_20190114.vcf.gz
EXAC_NONTCGA = $(ANN_DIR)/ExAC.r0.3.nontcga.sites.edited.vcf.gz
EXAC_NONPSYCH ?= $(ANN_DIR)/ExAC.r0.3.nonpsych.sites.edited.vcf.gz
#CADD_SNV ?= $(ANN_DIR)/CADD_whole_genome_SNVs_header.vcf.gz
#CADD_INDEL ?= $(ANN_DIR)/CADD_Indels_1000Gphase3_indels_header_sorted.vcf.gz
COSMIC ?= $(ANN_DIR)/CosmicCodingMuts_v87.vcf.gz
#COSMIC ?= $(ANN_DIR)/CosmicCodingMuts_v89_hg38.vcf

#HAPLOTYPE_INSUF_BED = $(ANN_DIR)/haploinsufficient.genes.bed
#CANCER_GENE_CENSUS_BED = $(ANN_DIR)/cancer_gene_census.bed
#KANDOTH_BED = $(ANN_DIR)/kandoth127smg.bed
#LAWRENCE_BED = $(ANN_DIR)/lawrence_cancer5000s.bed
#DGD_BED = $(ANN_DIR)/dgd_Hsa_all_v71.bed
CANCER_HOTSPOT_VCF = $(ANN_DIR)/hotspots_changv2_gao_nc.vcf.gz
GENE_SETS_GMT = $(ANN_DIR)/gene_ann_20200201.gmt
GENE_SETS_LIST = $(ANN_DIR)/gene_ann_20200201.list

CENTROMERE_TABLE = $(ANN_DIR)/centromeres_hg38.txt

DGV_FILE = $(ANN_DIR)/DGV_GRCh38_hg38_variants_2016-08-31.txt

GENCODE_GTF = $(ANN_DIR)/Homo_sapiens_assembly38__Gallus_gallus_000002315.5_GRCg6a.gtf
GENCODE_GENE_GTF = $(ANN_DIR)/Homo_sapiens_assembly38__Gallus_gallus_000002315.5_GRCg6a.gene.gtf
GENCODE_CODING_GENE_GTF = $(ANN_DIR)/Homo_sapiens_assembly38__Gallus_gallus_000002315.5_GRCg6a.gene.protein_coding.gtf
GENCODE_NONCODING_GENE_GTF = $(ANN_DIR)/Homo_sapiens_assembly38__Gallus_gallus_000002315.5_GRCg6a.gene.non_coding.gtf

RIBOSOMAL_INTERVALS ?= $(ANN_DIR)/Homo_sapiens_assembly38__Gallus_gallus_000002315.5_GRCg6a.annotation_ribosomal.interval_list
GENE_REF_FLAT = $(ANN_DIR)/Homo_sapiens_assembly38__Gallus_gallus_000002315.5_GRCg6a.refFlat

### TOOL-SPECIFIC REF ###


GISTIC_REF ?= $(REF_DIR)/gistic/hg38.UCSC.add_miR.160920.refgene.mat


ifeq ($(findstring scicore_ubuntu,$(HPC)),scicore_ubuntu)
  SNP_EFF_GENOME = GRCh38.99
else
  SNP_EFF_GENOME = GRCh38.86
endif

DB_NSFP ?= $(REF_DIR)/dbNSFP/GRCh38/dbNSFP3.2a.txt.gz

#MOSAICS_CHRFILE = $(GENOME_DIR)/human_g1k_v37_decoy.mosaics.chrfile

MUTSIGCV_REF_DIR ?= $(REF_DIR)/MutSigCV
MUTSIGCV_COVERAGE_REF ?= $(MUTSIGCV_REF_DIR)/coverage.exome_full192.txt
MUTSIGCV_COV_REF ?= $(MUTSIGCV_REF_DIR)/gene.covariates.txt
MUTSIGCV_DICT_REF ?= $(MUTSIGCV_REF_DIR)/mutation_type_dictionary_file.txt
MUTSIGCV_SEQ_REF_DIR ?= $(MUTSIGCV_REF_DIR)/chr_files_hg38/

#YOUN_AND_SIMON_SEQ_DATA ?= $(ANN_DIR)/youn_and_simon_sequence_data.RData

MSISENSOR_REF ?= $(REF_DIR)/msi_reference/hg38.txt
CNVKIT_GENE_RESOURCE ?= $(ANN_DIR)/cnvkit_gene_resource.txt

endif
HG38_GALLUS6A_INC = true

