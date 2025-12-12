# Bed files and annotation files
TARGETS_FILE_COVERED_INTERVALS ?= $(BED_DIR)/Agilent_All_Exon_V6+COSMIC_r2_S07604715_Covered.bed
TARGETS_FILE_INTERVALS ?= $(BED_DIR)/Agilent_All_Exon_V6+COSMIC_r2_S07604715_Padded.bed
TARGETS_FILE_MUTATIONS ?= $(BED_DIR)/Agilent_All_Exon_V6+COSMIC_r2_S07604715_Padded.bed
# TARGETS_FILE_GENES ?= $(BED_DIR)/wxs.GeneSymbols.txt
DBSNP_TARGETS_INTERVALS ?= $(BED_DIR)/Agilent_All_Exon_V6+COSMIC_r2_S07604715_dbsnp_146.vcf
DBSNP_TARGETS_INTERVALS_COMMON ?= $(BED_DIR)/Agilent_All_Exon_V6+COSMIC_r2_S07604715_dbsnp_146.common.monoallelic.vcf
# Resources
RESOURCE_REQ ?= medium

# Annotation and reporting parameters
# SNP_EFF_FLAGS ?= -canon -no NEXT_PROT -no-downstream -no-intergenic -no-upstream -no-utr
