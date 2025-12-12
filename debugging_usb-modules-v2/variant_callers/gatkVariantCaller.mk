# GATK variant caller module
# Author: Raymond Lim <raylim@mm.st> & Fong Chun Chan <fongchunchan@gmail.com>
# 

include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/variantCaller.inc

LOGDIR ?= log/gatk.$(NOW)

VPATH ?= bam

VARIANT_TYPES ?= gatk_snps gatk_indels
PHONY += gatk


#gatk : gatk/gatk_indels.filtered.vcf.gz gatk/gatk_snps.filtered.vcf.gz
gatk : $(foreach type,$(VARIANT_TYPES),vcf/all$(PROJECT_PREFIX).$(type).$(VCF_ANNS_SUFFIX).vcf) \
$(foreach type,$(VARIANT_TYPES),tables/all$(PROJECT_PREFIX).$(type).$(VCF_ANNS_SUFFIX).tab.txt)
# vcf/all$(PROJECT_PREFIX).gatk_indels.$(VCF_ANNS_SUFFIX).vcf.gz vcf/all$(PROJECT_PREFIX).gatk_snps.$(VCF_ANNS_SUFFIX).vcf.gz
#gatk : gatk_vcfs $(if $(findstring NONE,$(PANEL)),gatk_vcf_stats,gatk_tables)
#gatk_vcfs : $(foreach type,$(VARIANT_TYPES),$(call MAKE_VCF_FILE_LIST,$(type)) $(addsuffix .idx,$(call MAKE_VCF_FILE_LIST,$(type))))
#gatk_tables : $(foreach type,$(VARIANT_TYPES),$(call MAKE_TABLE_FILE_LIST,$(type)))
#gatk_reports : $(foreach type,$(VARIANT_TYPES),reports/$(type).dp_ft.grp)
#gatk_vcf_stats : $(foreach type,$(VARIANT_TYPES),$(call VCF_STATS,$(type)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

include usb-modules-v2/vcf_tools/vcftools.mk
include usb-modules-v2/variant_callers/gatk.mk

