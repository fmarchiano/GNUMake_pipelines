MUT_CALLER = manta

# Run manta on tumour-normal matched pairs

include usb-modules-v2/Makefile.inc

LOGDIR ?= log/manta.$(NOW)
PHONY += manta

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

manta : manta_vcfs
manta_vcfs : $(foreach pair,$(SAMPLE_PAIRS),manta/$(tumor.$(pair))/results/variants/somaticSV.vcf.gz)

define manta-tumor-normal
manta/$1/runWorkflow.py : bam/$1.bam bam/$2.bam
	$$(call RUN,1,1G,$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	$$(MANTA) configManta.py \
	--tumorBam $$< \
	--normalBam $$(<<) \
	--referenceFasta $$(REF_FASTA) \
	$$(if $$(findstring BAITS,$$(CAPTURE_METHOD)),--exome,) \
	$$(if $$(findstring hg38,$$(REF)),--callRegions $$(BED_DIR)/hg38_main_chr.bed.gz,) \
	--runDir $$(@D)")

# manta uses little RAM, 2G per cpu should be enough
manta/$1/results/variants/somaticSV.vcf.gz : manta/$1/runWorkflow.py
	$$(call RUN,$$(MANTA_NUM_CORES),$$(MANTA_MEM)G,$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(MANTA) $$< -j $$(MANTA_NUM_CORES) -g $$(MANTA_MEM)")

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call manta-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/vcf_tools/vcftools.mk
