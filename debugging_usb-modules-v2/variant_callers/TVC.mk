include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/variantCaller.inc

MUT_CALLER = tvc
TVC_NUM_CORES ?= 4

LOGDIR ?= log/tvc.$(NOW)

VPATH ?= bam
VARIANT_TYPES ?= tvc_snps tvc_indels

PHONY += tvc tvc_vcfs tvc_tables
tvc : tvc_vcfs tvc_tables
tvc_vcfs : $(foreach type,$(VARIANT_TYPES),$(call MAKE_VCF_FILE_LIST,$(type)) $(addsuffix .idx,$(call MAKE_VCF_FILE_LIST,$(type))))
tvc_tables : $(foreach type,$(VARIANT_TYPES),$(call MAKE_TABLE_FILE_LIST,$(type)))


# TVC used to produce TSVC_variants.vcf and TSVC_variants.vcf.gz, but in the current TVC version we only get TSVC_variants.vcf that is actually bgzipped and with a .tbi index.
# Rename it to avoid downstream errors
tvc/dbsnp/%/TSVC_variants.vcf.gz : bam/%.bam
	$(call RUN,$(TVC_NUM_CORES),$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_LONG),,"\
	$(TVC) -s $(DBSNP_TARGETS_INTERVALS) -i $< -r $(REF_FASTA) -o $(@D) \
	-N $(TVC_NUM_CORES) -p $(TVC_SOMATIC_JSON) -m $(TVC_MOTIF) \
	$(if $(TARGETS_FILE_INTERVALS),-b $(TARGETS_FILE_INTERVALS)) \
	-t $(TVC_ROOT_DIR) --primer-trim-bed $(PRIMER_TRIM_BED) &&\
	mv $(basename $@) $@ &&mv $(basename $@).tbi $@.tbi")


define tvc-vcf
tvc/vcf/$1/TSVC_variants.vcf.gz : bam/$1.bam bam/$1.bai
	$$(call RUN,$$(TVC_NUM_CORES),$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),,"\
	$$(TVC) -i $$< -r $$(REF_FASTA) -o $$(@D) -N $(TVC_NUM_CORES) \
	$$(if $(TARGETS_FILE_INTERVALS),-b $$(TARGETS_FILE_INTERVALS)) \
	-p $$(TVC_SOMATIC_JSON) -m $$(TVC_MOTIF) \
	-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) &&\
	mv $$(basename $$@) $$@ &&mv $$(basename $$@).tbi $$@.tbi")


tvc/vcf/$1/isec/0000.vcf : tvc/vcf/$1/TSVC_variants.multiallelic_ft.norm.left_align.vcf tvc/vcf/$1/TSVC_variants.biallelic_ft.vcf tvc/vcf/$1/TSVC_variants.multiallelic_ft.norm.left_align.vcf.gz tvc/vcf/$1/TSVC_variants.biallelic_ft.vcf.gz tvc/vcf/$1/TSVC_variants.multiallelic_ft.norm.left_align.vcf.gz.tbi tvc/vcf/$1/TSVC_variants.biallelic_ft.vcf.gz.tbi
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE),"\
	$$(BCFTOOLS) isec -O v -p $$(dir $$@) $$(word 3,$$^) $$(word 4,$$^)")

tvc/vcf/$1/isec/0001.vcf : tvc/vcf/$1/isec/0000.vcf
	
tvc/vcf/$1/isec/0003.vcf : tvc/vcf/$1/isec/0000.vcf
	

tvc/vcf/$1/TSVC_variants_final.vcf : tvc/vcf/$1/isec/0000.post_bcftools.vcf tvc/vcf/$1/isec/0001.post_bcftools.vcf tvc/vcf/$1/isec/0003.post_bcftools.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),$$(JAVA8_MODULE),"\
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) -R $$(REF_FASTA) \
	$$(foreach vcf,$$^,-V $$(vcf)) -o $$@ --assumeIdenticalSamples")

endef
$(foreach sample,$(SAMPLES), \
	$(eval $(call tvc-vcf,$(sample))))

tvc/vcf/%/TSVC_variants.snps.vcf : tvc/vcf/%/TSVC_variants_final.vcf
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(VCFTOOLS_MODULE),"\
	$(FIX_TVC_VCF) < $< | $(VCFTOOLS) --vcf - --remove-indels --recode --recode-INFO-all --out $@ && \
	mv $@.recode.vcf $@")

tvc/vcf/%/TSVC_variants.indels.vcf : tvc/vcf/%/TSVC_variants_final.vcf
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(VCFTOOLS_MODULE),"\
	$(FIX_TVC_VCF) < $< | $(VCFTOOLS) --vcf -  --keep-only-indels --recode --recode-INFO-all --out $@ && \
	mv $@.recode.vcf $@")

vcf/%.tvc_snps.vcf : tvc/vcf/%/TSVC_variants.snps.vcf
	$(INIT) $(VCF_SORT) $(REF_DICT) $< > $@

vcf/%.tvc_indels.vcf : tvc/vcf/%/TSVC_variants.indels.vcf
	$(INIT) $(VCF_SORT) $(REF_DICT) $< > $@

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

include usb-modules-v2/vcf_tools/vcftools.mk
