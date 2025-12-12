ifndef B37_HBV_HCV_GRCM38_INC
REF = b37_hbv_hcv_GRCm38
include usb-modules-v2/genome_inc/common.inc

REF_FASTA ?= $(GENOME_DIR)/human_g1k_v37_decoy_HBV_NC_003977.2_HCV_GRCm38.fasta
REF_DICT = $(REF_FASTA:.fasta=.dict)

DBSNP ?= $(ANN_DIR)/dbsnp_138.b37.vcf

ifeq ($(findstring NONE,$(PANEL)),NONE)
else
include usb-modules-v2/genome_inc/$(REF)/$(PANEL).inc
endif

RESOURCE_REQ ?= high
TARGETS_FILE_INTERVALS ?= $(REF_FASTA:.fasta=.bed)
DBSNP_TARGETS_INTERVALS ?= $(ANN_DIR)/dbsnp_138.b37_snps_monoallelic.vcf
TARGETS_FILE_GENES ?= $(ANN_DIR)/GRCh37_protein_coding_gene_symbols.txt

ifeq ($(findstring bed,$(TARGETS_FILE_INTERVALS)),bed)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE_INTERVALS) | cut -f1 | sort -V | uniq)
else
CHROMOSOMES ?= $(shell seq 1 22) X $(if $(findstring true,$(INCLUDE_CHR_Y)),Y) MT
endif

RECAL_CHR = 2

### ANNOTATION SETS ####

CONTEST_POP_FILE ?= $(ANN_DIR)/hg19_population_stratified_af_hapmap_3.3.vcf

KNOWN_INDELS ?= $(ANN_DIR)/Mills_and_1000G_gold_standard.indels.b37.vcf
OMNI ?= $(ANN_DIR)/1000G_omni2.5.b37.vcf
HAPMAP ?= $(ANN_DIR)/hapmap_3.3.b37.vcf
R1000G ?= $(ANN_DIR)/1000G_phase1.snps.high_confidence.b37.vcf

CLINVAR = $(ANN_DIR)/clinvar.vcf
EXAC_NONTCGA = $(ANN_DIR)/ExAC_nonTCGA.r0.3.1.sites.vep.edited.vcf
COSMIC ?= $(ANN_DIR)/CosmicCodingMuts.sorted.vcf

#HAPLOTYPE_INSUF_BED = $(ANN_DIR)/haploinsufficient.genes.bed
#CANCER_GENE_CENSUS_BED = $(ANN_DIR)/cancer_gene_census.bed
#KANDOTH_BED = $(ANN_DIR)/kandoth127smg.bed
#LAWRENCE_BED = $(ANN_DIR)/lawrence_cancer5000s.bed
#DGD_BED = $(ANN_DIR)/dgd_Hsa_all_v71.bed
GENE_SETS_GMT = $(ANN_DIR)/gene_ann_20200201.gmt
GENE_SETS_LIST = $(ANN_DIR)/gene_ann_20200201.list

CANCER_HOTSPOT_VCF = $(ANN_DIR)/hotspots_changv2_gao_nc.vcf

CENTROMERE_TABLE = $(ANN_DIR)/centromeres_b37.txt

DGV_FILE = $(ANN_DIR)/DGV_GRCh37_hg19_variants_2016-05-15.txt



GENCODE_GTF = $(ANN_DIR)/Homo_sapiens.GRCh37.82.chr_gsa_faidx_HBV_HCV_GRCm38.gtf
GENCODE_GENE_GTF = $(ANN_DIR)/Homo_sapiens.GRCh37.82.chr_gsa_faidx_HBV_HCV_GRCm38.genes.gtf
GENCODE_CODING_GENE_GTF = $(ANN_DIR)/Homo_sapiens.GRCh37.82.chr_gsa_faidx_HBV_HCV_GRCm38.genes.protein_coding.gtf
GENCODE_NONCODING_GENE_GTF = $(ANN_DIR)/Homo_sapiens.GRCh37.82.chr_gsa_faidx_HBV_HCV_GRCm38.genes.non_coding.gtf

RIBOSOMAL_INTERVALS ?= $(ANN_DIR)/ribosomalRNA.interval_list
GENE_REF_FLAT = $(ANN_DIR)/refFlat.txt

### TOOL-SPECIFIC REF ###
GISTIC_REF ?= $(REF_DIR)/gistic/hg19.UCSC.add_miR.140312.refgene.mat

SNP_EFF_GENOME = GRCh37.75
DB_NSFP ?= $(REF_DIR)/dbNSFP/GRCh37/dbNSFP2.9.txt.gz

MUTSIGCV_REF_DIR ?= $(REF_DIR)/MutSigCV
MUTSIGCV_COVERAGE_REF ?= $(MUTSIGCV_REF_DIR)/coverage.exome_full192.txt
MUTSIGCV_COV_REF ?= $(MUTSIGCV_REF_DIR)/gene.covariates.txt
MUTSIGCV_DICT_REF ?= $(MUTSIGCV_REF_DIR)/mutation_type_dictionary_file.txt
MUTSIGCV_SEQ_REF_DIR ?= $(MUTSIGCV_REF_DIR)/chr_files_hg19/

YOUN_AND_SIMON_SEQ_DATA ?= $(ANN_DIR)/youn_and_simon_sequence_data.RData

MSISENSOR_REF ?= $(REF_DIR)/msi_reference/b37.txt
CNVKIT_GENE_RESOURCE ?= $(ANN_DIR)/cnvkit_gene_resource_GRCh37.tsv

endif
B37_HBV_HCV_GRCM38_INC = true

