MUT_CALLER := tvc

include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

TVC_NUM_CORES ?= 4

LOGDIR ?= log/tvc_somatic.$(NOW)


PHONY += tvc_somatic tvc_somatic_vcfs tvc_somatic_tables 
tvc_somatic : tvc_somatic_vcfs tvc_somatic_tables 

VARIANT_TYPES ?= tvc_snps tvc_indels
tvc_somatic_vcfs : $(foreach type,$(VARIANT_TYPES),$(call MAKE_VCF_FILE_LIST,$(type)) $(addsuffix .idx,$(call MAKE_VCF_FILE_LIST,$(type))))
tvc_somatic_tables : $(foreach type,$(VARIANT_TYPES),$(call MAKE_TABLE_FILE_LIST,$(type)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)


define tvc-somatic-vcf
# TVC used to produce TSVC_variants.vcf and TSVC_variants.vcf.gz, but in the current TVC version we only get TSVC_variants.vcf that is actually bgzipped and with a .tbi index.
# Rename it to avoid downstream errors
tvc/vcf/$1_$2/TSVC_variants.vcf.gz : bam/$1.bam bam/$1.bam.bai bam/$2.bam bam/$2.bam.bai
	$$(call RUN,$$(TVC_NUM_CORES),$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_LONG),$(OPENBLAS_MODULE) $(PYTHON_MODULE),"\
	$$(TVC) -i $$< -n $$(word 3,$$^) -r $$(REF_FASTA) -o $$(@D) -N $$(TVC_NUM_CORES) \
	$$(if $$(TARGETS_FILE_INTERVALS),-b $$(TARGETS_FILE_INTERVALS)) -p $$(TVC_SOMATIC_JSON) -m $$(TVC_MOTIF) \
	-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) -g $$(basename $$(notdir $$<)) &&\
	mv $$(basename $$@) $$@ &&mv $$(basename $$@).tbi $$@.tbi")

tvc/vcf/$1_$2/isec/0000.vcf : tvc/vcf/$1_$2/TSVC_variants.multiallelic_ft.norm.left_align.vcf tvc/vcf/$1_$2/TSVC_variants.biallelic_ft.vcf tvc/vcf/$1_$2/TSVC_variants.multiallelic_ft.norm.left_align.vcf.gz tvc/vcf/$1_$2/TSVC_variants.biallelic_ft.vcf.gz tvc/vcf/$1_$2/TSVC_variants.multiallelic_ft.norm.left_align.vcf.gz.tbi tvc/vcf/$1_$2/TSVC_variants.biallelic_ft.vcf.gz.tbi
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE),"\
	$$(BCFTOOLS) isec -O v -p $$(dir $$@) $$(word 3,$$^) $$(word 4,$$^)")

tvc/vcf/$1_$2/isec/0001.vcf : tvc/vcf/$1_$2/isec/0000.vcf
	
tvc/vcf/$1_$2/isec/0003.vcf : tvc/vcf/$1_$2/isec/0000.vcf
	

tvc/vcf/$1_$2/TSVC_variants_final.vcf : tvc/vcf/$1_$2/isec/0000.post_bcftools.vcf tvc/vcf/$1_$2/isec/0001.post_bcftools.vcf tvc/vcf/$1_$2/isec/0003.post_bcftools.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),$$(JAVA8_MODULE),"\
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) \
	$$(foreach vcf,$$^,-V $$(vcf)) -o $$@ --assumeIdenticalSamples")

endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call tvc-somatic-vcf,$(tumor.$(pair)),$(normal.$(pair)))))

#include usb-modules-v2/vcf_tools/vcftools.mk
include usb-modules-v2/variant_callers/TVC.mk
#include usb-modules-v2/variant_callers/somatic/pon.mk
