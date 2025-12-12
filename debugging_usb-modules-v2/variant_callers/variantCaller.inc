ifndef VARIANT_CALLER_INC

### Filters

ifeq ($(findstring true,$(FILTER_VARIANTS)),true)
VCF_FILTERS += $(if $(TARGETS_FILE_INTERVALS),target_ft) dp_ft
VCF_FILTERS += $(if $(or $(findstring GRCm38,$(REF)),$(findstring SOMATIC,$(ANALYSIS_TYPE))),hotspot,) pass
else
VCF_FILTERS =
endif
VCF_FILTER_SUFFIX = $(subst $(space),.,$(strip $(call VCF_FILTERS)))

### Annotation
ifeq ($(findstring true,$(ANNOTATE_VARIANTS)),true)

VCF_ANNS = norm left_align post_bcftools $(if $(findstring GRCm38,$(REF)),eff dbsnp, \
	eff dbsnp $(if $(findstring SOMATIC,$(ANALYSIS_TYPE)),cosmic exac_nontcga gene_ann,) \
	$(if $(findstring GERMLINE,$(ANALYSIS_TYPE)),exac_nonpsych cadd,) clinvar nsfp)
else 
VCF_ANNS =
endif
VCF_ANNS_SUFFIX = $(subst $(space),.,$(strip $(VCF_ANNS)))

####

VCF_SUFFIX = $(subst $(space),.,$(strip $1 $(VCF_FILTER_SUFFIX) $(VCF_ANNS_SUFFIX)))

VCF_SUFFIXES = $(foreach type,$1,$(call VCF_SUFFIX,$(type)))

MAKE_VCF_FILE_LIST = $(foreach sample,$(SAMPLES),\
	$(foreach suff,$(call VCF_SUFFIXES,$1), \
	vcf/$(sample).$(suff).vcf))

ifeq ($(findstring false,$(INCLUDE_LNCRNA_IN_SUMMARY)),false)
TABLE_EFF_TYPES ?= nonsynonymous_synonymous_hotspot
else
TABLE_EFF_TYPES ?= nonsynonymous_synonymous_hotspot_lincRNA_upstream
endif

TABLE_EXTENSIONS = tab $(foreach eff,$(TABLE_EFF_TYPES),tab.$(eff))
TABLE_SUFFIXES = $(foreach suff,$(call VCF_SUFFIXES,$1),\
	 $(foreach ext,$(TABLE_EXTENSIONS),$(suff).$(ext)))

DOWNSTREAM_EFF_TYPES ?=	$(word 1,$(TABLE_EFF_TYPES))
DOWMSTREAM_VCF_TABLE_SUFFIX ?= $(call VCF_SUFFIXES,$1).tab.$(DOWNSTREAM_EFF_TYPES)

MAKE_TABLE_FILE_LIST = $(foreach sample,$(SAMPLES),\
	$(foreach suff,$(call TABLE_SUFFIXES,$1),\
	tables/$(sample).$(suff).txt))
MAKE_TABLE_FILE_LIST += $(foreach suff,$(call TABLE_SUFFIXES,$1),alltables/all$(PROJECT_PREFIX).$(suff).txt)

VCF_STATS = $(foreach sample,$(SAMPLES),\
	$(foreach suff,$(call VCF_SUFFIXES,$1), \
	vcf/$(sample).$(suff).vcf.stats))
VCF_STATS += $(foreach suff,$(call VCF_SUFFIXES,$1),alltables/all$(PROJECT_PREFIX).$(suff).vcf.stats)

endif
VARIANT_CALLER_INC = true

