MUT_CALLER = svaba

# Run SvABA on tumour-normal matched pairs

include usb-modules-v2/Makefile.inc

LOGDIR ?= log/svaba.$(NOW)
PHONY += svaba

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

svaba : svaba_vcfs
svaba_vcfs : $(foreach pair,$(SAMPLE_PAIRS),svaba/$(tumor.$(pair)).svaba.somatic.sv.vcf)

define svaba-tumor-normal
svaba/$1.contigs.bam : bam/$1.bam bam/$2.bam
	$$(call RUN,$$(SVABA_NUM_CORES),$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(SVABA) svaba run \
	-t $$< \
	-n $$(<<) \
	-p $$(SVABA_NUM_CORES) \
	-a svaba/$$(basename $$(notdir $$<)) \
	-G $$(REF_FASTA) \
	$$(if $$(findstring BAITS,$$(CAPTURE_METHOD)),--region $$(TARGETS_FILE_COVERED_INTERVALS),) \
	-D $$(DBSNP_TARGETS_INTERVALS) \
	--blacklist $$(if $$(findstring hg38,$$(REF)),$$(BED_DIR)/human.hg38.excl.bed,\
	$$(if $$(findstring b37,$$(REF)),$$(BED_DIR)/human.hg19.excl.bed,\
	$$(if $$(findstring GRCm38,$$(REF)),$$(BED_DIR)/mouse.mm10.excl.bed,)))")

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call svaba-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

%.svaba.somatic.sv.vcf : %.contigs.bam
	

include usb-modules-v2/vcf_tools/vcftools.mk
