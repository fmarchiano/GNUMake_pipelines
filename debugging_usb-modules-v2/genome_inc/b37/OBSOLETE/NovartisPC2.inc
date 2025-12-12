TARGETS_FILE_COVERED_INTERVALS ?= $(BED_DIR)/pc2_novartis_Covered_noHeader_noAnn.bed
TARGETS_FILE_INTERVALS ?= $(BED_DIR)/pc2_novartis_Padded_noHeader_noAnn.bed
TARGETS_FILE_MUTATIONS ?= $(BED_DIR)/pc2_novartis_Padded_noHeader_noAnn.bed
DBSNP_TARGETS_INTERVALS ?= $(BED_DIR)/pc2_novartis.dbsnp_138.vcf
RESOURCE_REQ ?= low
SNP_EFF_FLAGS ?= -canon -no NEXT_PROT -no-downstream -no-intergenic -no-upstream -no-utr

