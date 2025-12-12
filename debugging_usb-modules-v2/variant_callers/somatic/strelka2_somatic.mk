MUT_CALLER = strelka2_somatic

# Run strelka2 on tumour-normal matched pairs

include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

##### DEFAULTS ######
STRELKA2_NUM_CORES ?= 8

LOGDIR ?= log/strelka2_somatic.$(NOW)
PHONY += strelka2 strelka2_vcfs strelka2_tables

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

strelka2 : strelka2_vcfs strelka2_tables

VARIANT_TYPES := strelka2_snvs strelka2_indels

strelka2_vcfs : $(foreach type,$(VARIANT_TYPES),$(call MAKE_VCF_FILE_LIST,$(type)))
strelka2_tables : $(foreach type,$(VARIANT_TYPES),$(call MAKE_TABLE_FILE_LIST,$(type)))

define strelka2-tumor-normal
# If USE_BAM_CLIPOVERLAP then look for BAMs in bam_clipoverlap folder.
strelka2/$1_$2/runWorkflow.py : $(if $(findstring true,$(USE_BAM_CLIPOVERLAP)),bam_clipoverlap,bam)/$1.bam $(if $(findstring true,$(USE_BAM_CLIPOVERLAP)),bam_clipoverlap,bam)/$2.bam
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(PERL_MODULE) $$(SINGULARITY_MODULE),"\
	rm -rf $$(@D) && \
	$$(CONFIGURE_STRELKA2_SOMATIC) $$(if $$(findstring NONE,$$(PANEL)),,--exome) --tumorBam $$< --normalBam $$(<<) \
	--referenceFasta $$(REF_FASTA) --config $$(STRELKA2_SOMATIC_CONFIG) --runDir $$(@D)")

# Because we use 8 cores, mem per cpu request becomes unreasonably high
# (and job won't run if RESOURCE_REQ_LOW_MEM is set to a very high value).
# Use a fixed mem value instead (strelka2 uses very little RAM). For most WGS, even 2GB is enough
strelka2/$1_$2/results/variants/somatic.indels.vcf.gz : strelka2/$1_$2/runWorkflow.py
	$$(call RUN,$$(STRELKA2_NUM_CORES),4GB,$$(RESOURCE_REQ_SHORT),$$(SINGULARITY_MODULE),"\
	$$(STRELKA2) $$(<D)/runWorkflow.py -j $$(STRELKA2_NUM_CORES) -m local")

strelka2/$1_$2/results/variants/somatic.snvs.vcf.gz : strelka2/$1_$2/results/variants/somatic.indels.vcf.gz


vcf/$1_$2.%.vcf : strelka2/vcf/$1_$2.%.vcf
	$$(INIT) perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' $$< > $$@ && $$(RM) $$<

strelka2/vcf/$1_$2.strelka2_snvs.vcf : strelka2/$1_$2/results/variants/somatic.snvs.vcf.gz
	$$(INIT) $$(FIX_STRELKA2_VCF_SNVS) $$< > $$@

strelka2/vcf/$1_$2.strelka2_indels.vcf : strelka2/$1_$2/results/variants/somatic.indels.vcf.gz
	$$(INIT) $$(FIX_STRELKA2_VCF_INDELS) $$< > $$@

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call strelka2-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/vcf_tools/vcftools.mk

