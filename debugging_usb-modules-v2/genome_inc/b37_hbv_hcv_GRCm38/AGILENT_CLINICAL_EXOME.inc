# Bed files and annotation files
TARGETS_FILE_COVERED_INTERVALS ?= $(BED_DIR)/S06588914_Agilent_Clinical_Exome_Covered_noHeader_noAnn.bed
TARGETS_FILE_INTERVALS ?= $(BED_DIR)/S06588914_Agilent_Clinical_Exome_Padded_noHeader_noAnn.bed
TARGETS_FILE_MUTATIONS ?= $(BED_DIR)/S06588914_Agilent_Clinical_Exome_Padded_noHeader_noAnn.bed
TARGETS_FILE_GENES ?= $(BED_DIR)/S06588914_Agilent_Clinical_Exome_TargetGeneSymbols.txt
DBSNP_TARGETS_INTERVALS ?= $(BED_DIR)/S06588914_Agilent_Clinical_Exome.dbsnp_138.vcf
DBSNP_TARGETS_INTERVALS_COMMON ?= $(BED_DIR)/S06588914_Agilent_Clinical_Exome.dbsnp_138.common.vcf

# Resources
RESOURCE_REQ ?= medium

# Annotation and reporting parameters
# SNP_EFF_FLAGS ?= -canon -no NEXT_PROT -no-downstream -no-intergenic -no-upstream -no-utr

