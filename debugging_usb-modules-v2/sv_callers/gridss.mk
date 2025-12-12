#template for gridss
MUT_CALLER = gridss

# Run gridss on tumour-normal matched pairs

include usb-modules-v2/Makefile.inc

LOGDIR ?= log/gridss.$(NOW)
PHONY += gridss

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

gridss : gridss_vcfs

gridss_vcfs : $(foreach pair,$(SAMPLE_PAIRS),gridss/$(tumor.$(pair)).somatic.vcf)

# Somatic analysis depends on PoN files
define gridss-tumor-normal

# Main call for tumor-normal pairs
gridss/$1.vcf : bam/$2.bam bam/$1.bam
	mkdir -p gridss/gridss_tmp && \
	$$(call RUN,8,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(GRIDSS) gridss \
	$$(if $$(findstring hg38,$$(REF)),-b $$(BED_DIR)/ENCFF356LFX.bed) \
	$$(if $$(findstring b37,$$(REF)),-b $$(BED_DIR)/ENCFF001TDO.bed) \
	-r $$(REF_FASTA) \
	-s preprocess$$(,)assemble$$(,)call \
	--skipsoftcliprealignment \
	--workingdir gridss/gridss_tmp \
	-o $$@ \
	$$^")

# Somatic filtering
gridss/$1.somatic.vcf : gridss/$1.vcf $(SV_PON_DIR)
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(GRIDSS) Rscript $$(IMG_DIR)/GRIDSS/gridss-2.13.2/scripts/gridss_somatic_filter \
	--pondir $$(<<) \
	--input $$< \
	--output $$@ \
	--fulloutput gridss/$1.full.vcf.gz \
	--scriptdir $$(IMG_DIR)/GRIDSS/gridss-2.13.2/scripts/")

endef

$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call gridss-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))