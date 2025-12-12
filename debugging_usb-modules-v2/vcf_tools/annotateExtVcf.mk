include usb-modules-v2/Makefile.inc

CALLER_PREFIX ?=
FILTER_VARIANTS = false
ANNOTATE_VARIANTS ?= true

ifeq ($(findstring SOMATIC,$(ANALYSIS_TYPE)),SOMATIC)
USE_SUFAM = false
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc
else
include usb-modules-v2/variant_callers/variantCaller.inc
endif

# CALLER_PREFIX should be mutect or anything that goes before the filtering suffix
# if filtering is required then set FILTER_VARIANTS to true
# or if custom filtering is required then specify VCF_FILTER_SUFFIX
# if annotation is required then set ANNOTATE_VARIANTS to true
# or for custom annotation, specificy VCF_ANNS_SUFFIX

LOGDIR ?= log/ann_ext_vcf.$(NOW)

ann_ext_vcf : ext_vcfs ext_tables

ext_vcfs : $(foreach ext_name,$(CALLER_PREFIX),$(call MAKE_VCF_FILE_LIST,$(ext_name)))
ext_tables : $(foreach ext_name,$(CALLER_PREFIX),$(call MAKE_TABLE_FILE_LIST,$(ext_name)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : ann_ext_vcf

include usb-modules-v2/vcf_tools/vcftools.mk
