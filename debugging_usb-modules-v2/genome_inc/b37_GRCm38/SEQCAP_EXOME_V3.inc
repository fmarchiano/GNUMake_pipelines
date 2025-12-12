# Bed files and annotation files
TARGETS_FILE_COVERED_INTERVALS ?= $(BED_DIR)/SeqCap_EZ_Exome_v3_hg19_primary_targets.bed
TARGETS_FILE_INTERVALS ?= $(TARGETS_FILE_COVERED_INTERVALS)
TARGETS_FILE_INTERVALS ?= $(TARGETS_FILE_COVERED_INTERVALS)
TARGETS_FILE_MUTATIONS ?= $(TARGETS_FILE_COVERED_INTERVALS)
DBSNP_TARGETS_INTERVALS ?= $(BED_DIR)/SeqCap_EZ_Exome_v3_hg19_primary_targets.dbsnp_138.vcf
DBSNP_TARGETS_INTERVALS_COMMON ?= $(BED_DIR)SeqCap_EZ_Exome_v3_hg19_primary_targets.dbsnp_138.common.vcf

# Resources
RESOURCE_REQ ?= high

# Annotation and reporting parameters
# SNP_EFF_FLAGS ?= -canon -no NEXT_PROT -no-downstream -no-intergenic -no-upstream -no-utr

