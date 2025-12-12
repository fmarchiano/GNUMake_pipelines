

include usb-modules-v2/Makefile.inc
include usb-modules-v2/aligners/align.inc

ALIGNER := bwaaln
BWA_ALN_OPTS ?=
BWA_NUM_CORES ?= 8

LOGDIR ?= log/bwaaln.$(NOW)

VPATH ?= unprocessed_bam

..DUMMY := $(shell mkdir -p version; $(BWA) &> version/bwaaln.txt; echo "options: $(BWA_ALN_OPTS)" >> version/bwaaln.txt )
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: bwaaln

ifeq ($(strip $(PRIMARY_ALIGNER)),bwaaln)
BWA_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
else
BWA_BAMS = $(foreach sample,$(SAMPLES),bwaaln/bam/$(sample).bwaaln.bam)
endif

bwaaln : $(BWA_BAMS) $(addsuffix .bai,$(BWA_BAMS))

bam/%.bam : bwaaln/bam/%.bwaaln.$(BAM_SUFFIX)
	$(INIT) ln -f $< $@

bwaaln/sai/%.sai : fastq/%.fastq.gz
	$(call RUN,$(BWA_NUM_CORES),$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_MEDIUM),$(BWA_MODULE),"\
	$(BWA_ALN) $(BWA_ALN_OPTS) -t $(BWA_NUM_CORES) $(REF_FASTA) $(<) > $(@) ")

ifeq ($(MERGE_SPLIT_FASTQ),true)
bwaaln/bam/%.bwaaln.bam : bwaaln/sai/%.1.sai $(if $(findstring true,$(PAIRED_END)),bwaaln/sai/%.2.sai) \
fastq/%.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),fastq/%.2.fastq.gz)
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9]\+//'`; \
	$(call RUN,$(BWA_NUM_CORES),$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_MEDIUM),$(BWA_MODULE) $(SAMTOOLS_MODULE),"\
	$(if $(findstring true,$(PAIRED_END)),$(BWA_SAMPE) -P,$(BWA_SAMSE)) \
	-r \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(REF_FASTA) $^ | \
	$(SAMTOOLS) view -uhS - | $(SAMTOOLS) sort -T $(TMPDIR)/$(@F).tmp -m 20G -o $@")
else
define align-split-fastq
bwaaln/bam/$1.bwaaln.bam : bwaaln/sai/$1.1.sai $(if $(findstring true,$(PAIRED_END)),bwaaln/sai/$1.2.sai) \
fastq/$1.1.fastq.gz $$(if $$(findstring true,$$(PAIRED_END)),fastq/$1.2.fastq.gz)
	LBID=`echo "$1" | sed 's/_[A-Za-z0-9]\+//'`; \
	$$(call RUN,$$(BWA_NUM_CORES),$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(BWA_MODULE) $$(SAMTOOLS_MODULE),"\
	$$(if $$(findstring true,$$(PAIRED_END)),$$(BWA_SAMPE) -P,$$(BWA_SAMSE)) \
	-r \"@RG\tID:$1\tLB:$$$${LBID}\tPL:$${SEQ_PLATFORM}\tSM:$$$${LBID}\" $$(REF_FASTA) $$^ | \
	$$(SAMTOOLS) view -uhS - | $$(SAMTOOLS) sort -T $$(TMPDIR)/$$(@F).tmp -m 20G -o $$@")
endef
$(foreach sample,$(SPLIT_SAMPLES),$(foreach split,$(split.$(sample)),$(eval $(call align-split-fastq,$(split)))))

define merged-bwaaln-bam
bwaaln/bam/$1.bwaaln.header.sam : $$(foreach split,$2,bwaaln/bam/$$(split).bwaaln.bam)
	$$(INIT) module load $$(SAMTOOLS_MODULE); $$(SAMTOOLS) view -H $$< | grep -v '^@RG' > $$@.tmp; \
	for bam in $$^; do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $$(RM) $$@.tmp

bwaaln/bam/$1.bwaaln.bam : bwaaln/bam/$1.bwaaln.header.sam $$(foreach split,$2,bwaaln/bam/$$(split).bwaaln.bam) 
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SAMTOOLS_MODULE),"\
	$$(SAMTOOLS) merge -f -h $$< $$(@) $$(filter %.bam,$$^) && $$(RM) $$<")
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-bwaaln-bam,$(sample),$(split.$(sample)))))


endif

include usb-modules-v2/fastq_tools/fastq.mk
include usb-modules-v2/bam_tools/processBam.mk
include usb-modules-v2/aligners/align.mk
include usb-modules-v2/variant_callers/gatk.mk
