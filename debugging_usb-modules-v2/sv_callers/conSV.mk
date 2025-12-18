.DEFAULT_GOAL := conSV

# run brass
MUT_CALLER = conSV

include usb-modules-v2/Makefile.inc

LOGDIR = log/conSV.$(NOW)
PHONY += conSV

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

SV_CALLERS ?= brass delly gridss manta svaba
SV_CALLERS := $(strip $(SV_CALLERS))
SV_CALLER_STRING := $(subst $(space),.,$(SV_CALLERS))
MIN_CALLERS ?= 2
MIN_DISTANCE_BP ?= 200
SINGULARITY_BIND = $(IMG_DIR)/oncoliner/

$(info sample pair: $(SAMPLE_PAIRS))
$(info callers: $(SV_CALLERS))

conSV : $(foreach pair,$(SAMPLE_PAIRS), $(foreach caller,$(SV_CALLERS), conSV/$(caller)/$(tumor.$(pair)).vcf ))

# create folder and organise data
#BRASS no filter implemented so taken as it is
define preprocess_brass
conSV/brass/$1.vcf : brass/$1_$2/$1_vs_$2.annot.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	cp $$< $$@")
endef

#the other callers have PASS filter implemented
define preprocess_delly
conSV/delly/$1.vcf : delly/$1.somatic.SVpass.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	cp $$< $$@")
endef

define preprocess_gridss
conSV/gridss/$1.vcf : gridss/$1.somatic.SVpass.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	cp $$< $$@")
endef

define preprocess_manta
conSV/manta/$1.vcf : manta/$1/results/variants/somaticSV.SVpass.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	cp $$< $$@")
endef

define preprocess_svaba
conSV/svaba/$1.vcf : svaba/$1.svaba.somatic.sv.SVpass.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	cp $$< $$@")
endef

# Preprocess step for SV callers to populate the fodler and rename files appropriately
define preprocess_dispatch
$(info Preprocessing SV caller '$3' for sample pair: $1 (tumor), $2 (normal))
$(if $(filter brass,$3),$(call preprocess_brass,$1,$2))
$(if $(filter delly,$3),$(call preprocess_delly,$1,$2))
$(if $(filter gridss,$3),$(call preprocess_gridss,$1,$2))
$(if $(filter manta,$3),$(call preprocess_manta,$1,$2))
$(if $(filter svaba,$3),$(call preprocess_svaba,$1,$2))
$(if $(filter-out brass delly gridss manta svaba,$3),$(error Unknown SV caller '$3'))
endef
$(foreach pair,$(SAMPLE_PAIRS), \
  $(foreach caller,$(SV_CALLERS), \
    $(eval $(call preprocess_dispatch,$(tumor.$(pair)),$(normal.$(pair)),$(caller))) \
  ) \
)

# Preprocess step to normalize and annotate each SV caller's filtered VCF
define variant_preprocess
conSV/variantExtr/$1_$2.$3.varExtr.vcf : conSV/$3/$1_$2.%.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	$$(VARIANT_EXTRACTOR) python3 VariantExtractor.py $$< $$@")
endef
# Loop through each sample pair and caller to preprocess each VCF
#$(foreach pair,$(SAMPLE_PAIRS), \
	$(foreach caller,$(SV_CALLERS), \
		$(eval $(call variant_preprocess,$(tumor.$(pair)),$(normal.$(pair)),$(caller))) \
	) \
)