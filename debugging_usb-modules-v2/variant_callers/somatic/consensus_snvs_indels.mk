include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

# Usage: make TYPE=snvs or make TYPE=indels
TYPE ?= snvs
CALLERS ?= strelka2_$(TYPE) mutect2
CALLERS := $(strip $(CALLERS))
CALLER_STRING := $(subst $(space),.,$(CALLERS))
MIN_CALLERS ?= 2

LOGDIR ?= log/consensus_$(TYPE).$(NOW)

PHONY += all

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

#clean intermediate files 
all: consensus_tables consensus_vcfs
	rm -f vcf/*.norm2.vcf.gz \
		vcf/*.norm2.vcf.gz.csi \
		vcf/*consensus.$(CALLER_STRING).vcf \
		vcf/*fail_indels.$(CALLER_STRING).vcf \
		vcf/*union*

ifeq ($(TYPE),indels)
consensus_vcfs : $(call MAKE_VCF_FILE_LIST,consensus.$(CALLER_STRING)) $(call MAKE_VCF_FILE_LIST,fail_indels.$(CALLER_STRING))
consensus_tables : $(call MAKE_TABLE_FILE_LIST,consensus.$(CALLER_STRING)) $(call MAKE_TABLE_FILE_LIST,fail_indels.$(CALLER_STRING))
else
consensus_vcfs : $(call MAKE_VCF_FILE_LIST,consensus.$(CALLER_STRING))
consensus_tables : $(call MAKE_TABLE_FILE_LIST,consensus.$(CALLER_STRING))
endif


# Preprocess step to normalize and annotate each caller's VCF
define preprocess
vcf/$1_$2.$3.$(VCF_FILTER_SUFFIX).norm2.vcf.gz : vcf/$1_$2.$3.$(VCF_FILTER_SUFFIX).vcf
vcf/$1_$2.$3.$(VCF_FILTER_SUFFIX).norm2.vcf.gz.csi : vcf/$1_$2.$3.$(VCF_FILTER_SUFFIX).norm2.vcf.gz
endef
# Loop through each sample pair and caller to preprocess each VCF
$(foreach pair,$(SAMPLE_PAIRS), \
	$(foreach caller,$(CALLERS), \
		$(eval $(call preprocess,$(tumor.$(pair)),$(normal.$(pair)),$(caller))) \
	) \
)

# Generate consensus VCF for each sample pair and all callers with bcftools isec
define consensus_snvs
vcf/$1_$2.consensus.$(CALLER_STRING).vcf : $(foreach caller,$(CALLERS),vcf/$1_$2.$(caller).$(VCF_FILTER_SUFFIX).norm2.vcf.gz.csi)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE),"\
	$$(BCFTOOLS) isec -n+$$(MIN_CALLERS) -w1 $$(basename $$^) -o $$@")

vcf/$1_$2.consensus.$(CALLER_STRING).$(VCF_FILTER_SUFFIX).$(VCF_ANNS_SUFFIX).vcf : vcf/$1_$2.consensus.$(CALLER_STRING).vcf
endef

# Generate consensus VCF for each sample pair and all callers with bcftools isec, in case of indels also keeping the excluded entries
#PS: exluded files for indels are pretty big, in the future we could consider using SURVIVOR and allow a certain slope for indel coordinates differences instead
define consensus_indels
vcf/$1_$2.consensus.$(CALLER_STRING).vcf : $(foreach caller,$(CALLERS),vcf/$1_$2.$(caller).$(VCF_FILTER_SUFFIX).norm2.vcf.gz.csi)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE),"\
	$$(BCFTOOLS) isec -n+$$(MIN_CALLERS) -w1 $$(basename $$^) -o $$@")

vcf/$1_$2.consensus.$(CALLER_STRING).norm2.vcf.gz : vcf/$1_$2.consensus.$(CALLER_STRING).vcf
vcf/$1_$2.consensus.$(CALLER_STRING).norm2.vcf.gz.csi : vcf/$1_$2.consensus.$(CALLER_STRING).norm2.vcf.gz

vcf/$1_$2.union.$(CALLER_STRING).vcf : $(foreach caller,$(CALLERS),vcf/$1_$2.$(caller).$(VCF_FILTER_SUFFIX).norm2.vcf.gz.csi)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE),"\
	$$(BCFTOOLS) isec -n+1 -w1 $$(basename $$^) -o $$@")

vcf/$1_$2.union.$(CALLER_STRING).norm2.vcf.gz : vcf/$1_$2.union.$(CALLER_STRING).vcf
vcf/$1_$2.union.$(CALLER_STRING).norm2.vcf.gz.csi : vcf/$1_$2.union.$(CALLER_STRING).norm2.vcf.gz

vcf/$1_$2.fail_indels.$(CALLER_STRING).vcf : vcf/$1_$2.union.$(CALLER_STRING).norm2.vcf.gz.csi vcf/$1_$2.consensus.$(CALLER_STRING).norm2.vcf.gz.csi
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE),"\
	$$(BCFTOOLS) isec -n=1 -w1 $$(basename $$<) $$(basename $$<<) -o $$@")

vcf/$1_$2.fail_indels.$(CALLER_STRING).$(VCF_FILTER_SUFFIX).$(VCF_ANNS_SUFFIX).vcf : vcf/$1_$2.fail_indels.$(CALLER_STRING).vcf
vcf/$1_$2.consensus.$(CALLER_STRING).$(VCF_FILTER_SUFFIX).$(VCF_ANNS_SUFFIX).vcf : vcf/$1_$2.consensus.$(CALLER_STRING).vcf
endef

# Loop through each sample pair and generate consensus VCFs for snvs and indels accordingly
ifeq ($(TYPE),snvs)
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call consensus_snvs,$(tumor.$(pair)),$(normal.$(pair)))))
endif

ifeq ($(TYPE),indels)
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call consensus_indels,$(tumor.$(pair)),$(normal.$(pair)))))
endif

include usb-modules-v2/vcf_tools/vcftools.mk
