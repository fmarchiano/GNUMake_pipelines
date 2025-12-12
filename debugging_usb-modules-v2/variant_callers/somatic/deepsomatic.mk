MUT_CALLER = deepsomatic

# Run deepsomatic on tumour-normal matched pairs

include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/deepsomatic.$(NOW)
PHONY += deepsomatic

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

deepsomatic : deepsomatic_vcfs deepsomatic_tables

deepsomatic_vcfs : $(call MAKE_VCF_FILE_LIST,deepsomatic) $(addsuffix .idx,$(call MAKE_VCF_FILE_LIST,deepsomatic))
deepsomatic_tables : $(call MAKE_TABLE_FILE_LIST,deepsomatic)

# TUMOR/NORMAL
ifeq ($(TUMOR_ONLY),false)
define deepsomatic-tumor-normal
deepsomatic/$1_$2.vcf.gz : bam/$1.bam bam/$2.bam
	$$(call RUN,$$(DEEPSOMATIC_NUM_CORES),$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_LONG),$$(SINGULARITY_MODULE),"\
	$$(DEEPSOMATIC) run_deepsomatic \
	--num_shards $$(DEEPSOMATIC_NUM_CORES) \
	--ref=$$(REF_FASTA) \
	--reads_tumor=$$< \
	--reads_normal=$$(<<) \
	--output_vcf=$$@ \
	--model_type=$$(DEEPSOMATIC_MODEL) \
	--logging_dir=$$(basename $$(basename $$@)).log --runtime_report")

vcf/$1_$2.deepsomatic.vcf : deepsomatic/$1_$2.deepsomatic.vcf.gz
	$$(INIT) zcat $$< > $$@

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call deepsomatic-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

# TUMOR_ONLY
ifeq ($(TUMOR_ONLY),true)
define deepsomatic-tumor-only
deepsomatic/$1.vcf.gz : bam/$1.bam
	$$(call RUN,$$(DEEPSOMATIC_NUM_CORES),$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_LONG),$$(SINGULARITY_MODULE),"\
	$$(DEEPSOMATIC) run_deepsomatic \
	--use_default_pon_filtering=true \
	--num_shards $$(DEEPSOMATIC_NUM_CORES) \
	--ref=$$(REF_FASTA) \
	--reads_tumor=$$< \
	--output_vcf=$$@ \
	--model_type=$$(DEEPSOMATIC_MODEL) \
	--logging_dir=$$(basename $$(basename $$@)).log --runtime_report")

# pipeline expects FA, but we have VAF
vcf/$1.deepsomatic.vcf : deepsomatic/$1.vcf.gz
	$$(call RUN,1,1G,$$(RESOURCE_REQ_VSHORT),,"\
	zcat $$< | sed 's/VAF/FA/' > $$@")

endef
$(foreach sample,$(SAMPLES),$(eval $(call deepsomatic-tumor-only,$(sample))))
endif



include usb-modules-v2/vcf_tools/vcftools.mk
