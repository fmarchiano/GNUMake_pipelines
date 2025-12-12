# Bed files and annotation files
TARGETS_FILE_COVERED_INTERVALS ?= $(BED_DIR)/Agilent_CRE_V2_S30409818_Covered.bed
TARGETS_FILE_INTERVALS ?= $(BED_DIR)/Agilent_CRE_V2_S30409818_Padded.bed
TARGETS_FILE_MUTATIONS ?= $(BED_DIR)/Agilent_CRE_V2_S30409818_Padded.bed
TARGETS_FILE_GENES ?= $(BED_DIR)/S06588914_Agilent_Clinical_Exome_TargetGeneSymbols.txt
DBSNP_TARGETS_INTERVALS ?= $(BED_DIR)/Agilent_CRE_V2_S30409818_Padded.dbsnp_138.vcf
DBSNP_TARGETS_INTERVALS_COMMON ?= $(BED_DIR)/Agilent_CRE_V2_S30409818_Padded.dbsnp_138.common.monoallelic.vcf

# Resources
RESOURCE_REQ ?= medium

# Annotation and reporting parameters
# SNP_EFF_FLAGS ?= -canon -no NEXT_PROT -no-downstream -no-intergenic -no-upstream -no-utr
