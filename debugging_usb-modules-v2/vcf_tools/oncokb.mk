include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

# Usage: consensus or single caller

TUMOR_TYPE ?= default
TYPE ?= snvs
TARGET_TYPE ?= gene
CONSENSUS ?= consensus
CALLER ?= strelka2_$(TYPE) mutect2
CALLER := $(strip $(CALLER))
CALLER_STRING := $(subst $(space),.,$(CALLER))

LOGDIR ?= log/oncokb.$(NOW)

PHONY += all

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

all: results

ifeq ($(CONSENSUS),consensus)
  CALLER_SUFFIX := consensus.$(CALLER_STRING)
else
  CALLER_SUFFIX := $(CALLER_STRING)
endif

ifeq ($(filter $(TYPE),snvs indels),$(TYPE))
results : oncokb/snvs_indels/allTN.$(CALLER_SUFFIX).$(TUMOR_TYPE).oncokb.maf
else ifeq ($(TYPE),cna)
results : oncokb/cna/all_thresholded.by_genes.$(TUMOR_TYPE).oncokb.txt
else ifeq ($(TYPE),fusion)
results : oncokb/fusion/SV.consensus.fusion.$(TUMOR_TYPE).oncokb.tsv
else ifeq ($(TYPE),sv)
results : oncokb/sv/SV.consensus.distrupted_$(TARGET_TYPE).$(TUMOR_TYPE).oncokb.tsv
endif

# convert VCF to MAF
define convert_vcf_maf
maf/$1_$2.$(CALLER_SUFFIX).$(call VCF_SUFFIX,).maf : vcf/$1_$2.$(CALLER_SUFFIX).$(call VCF_SUFFIX,).vcf
endef

# oncokb MAF annotator
define oncokb_snvs_indels_annotator
oncokb/snvs_indels/$1_$2.$(CALLER_SUFFIX).$(TUMOR_TYPE).oncokb.maf : maf/$1_$2.$(CALLER_SUFFIX).$(call VCF_SUFFIX,).maf
	$$(MKDIR) $$(@D)
	$$(call RUN,4,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(VEP_ONCOKB_MODULE)," 	sleep 5 && python $$(ONCOKB_MAF_ANNO) -i $$< -o $$@ $$(if $$(findstring hg38,$$(REF)), -r GRCh38) -t $$(TUMOR_TYPE) -b $$(ONCOKB_TOKEN) -d")
endef

# oncokb CNA annotator
define oncokb_cna_annotator
oncokb/cna/all_thresholded.by_genes.$(TUMOR_TYPE).oncokb.txt : gistic/gistic_cnv300000/all_thresholded.by_genes.txt
	$$(MKDIR) $$(@D)
	$$(call RUN,4,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(VEP_ONCOKB_MODULE),"\
	sleep 5 && python $$(ONCOKB_CNA_ANNO) -i $$< -o $$@ -t $$(TUMOR_TYPE) -b $$(ONCOKB_TOKEN) -z -d")
endef

#need to format fusion in the right format for oncokb
# oncokb FUSION annotator
define oncokb_fusion_annotator
oncokb/fusion/SV.consensus.fusion.$(TUMOR_TYPE).oncokb.tsv : consSV/oncokb_input/SV.consensus.fusion.tsv
	$$(MKDIR) $$(@D)
	$$(call RUN,4,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(VEP_ONCOKB_MODULE),"\
	sleep 5 && python $$(ONCOKB_FUSION_ANNO) -i $$< -o $$@ -t $$(TUMOR_TYPE) -b $$(ONCOKB_TOKEN) -d")
endef
#need to format gene/promoters in the right format for oncokb
# oncokb SV annotator
define oncokb_sv_annotator
oncokb/sv/SV.consensus.distrupted_$(TARGET_TYPE).$(TUMOR_TYPE).oncokb.tsv : consSV/oncokb_input/SV.consensus.distrupted_$(TARGET_TYPE).tsv
	$$(MKDIR) $$(@D)
	$$(call RUN,4,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(VEP_ONCOKB_MODULE),"\
	sleep 5 && python $$(ONCOKB_SV_ANNO) -i $$< -o $$@ -t $$(TUMOR_TYPE) -b $$(ONCOKB_TOKEN) -d")
endef

#execute rules based on TYPE
ifeq ($(filter $(TYPE),snvs indels),$(TYPE))
$(foreach pair,$(SAMPLE_PAIRS), \
		$(eval $(call convert_vcf_maf,$(tumor.$(pair)),$(normal.$(pair)))) \
)

$(foreach pair,$(SAMPLE_PAIRS), \
		$(eval $(call oncokb_snvs_indels_annotator,$(tumor.$(pair)),$(normal.$(pair)))) \
)

oncokb/snvs_indels/allTN.$(CALLER_SUFFIX).$(TUMOR_TYPE).oncokb.maf  : $(foreach pair,$(SAMPLE_PAIRS),oncokb/snvs_indels/$(pair).$(CALLER_SUFFIX).$(TUMOR_TYPE).oncokb.maf)
	$(INIT) \
	{ \
	echo -e "TUMOR_NORMAL\t$(shell grep -v '^#' $< | sed -n 1p)"; \
	for x in $^; do \
		pair=$$(basename $$x | cut -d. -f1); \
		grep -v '^#' $$x | sed 1d | awk -v p=$$pair '{print p "\t" $$0}'; \
	done; \
	} > $@

else ifeq ($(TYPE),cna)
$(eval $(oncokb_cna_annotator))
else ifeq ($(TYPE),fusion)
$(eval $(oncokb_fusion_annotator))
else ifeq ($(TYPE),sv)
$(eval $(oncokb_sv_annotator))
endif

include usb-modules-v2/vcf_tools/vcftools.mk
