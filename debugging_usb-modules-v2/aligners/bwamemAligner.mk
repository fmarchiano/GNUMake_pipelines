
include usb-modules-v2/Makefile.inc
include usb-modules-v2/aligners/align.inc

ALIGNER := bwamem
BWA_MEM_OPTS ?= -M
BWA_NUM_CORES ?= 8

LOGDIR ?= log/bwamem.$(NOW)

VPATH ?= unprocessed_bam

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: bwamem

ifeq ($(strip $(PRIMARY_ALIGNER)),bwamem)
BWA_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
else
BWA_BAMS = $(foreach sample,$(SAMPLES),bwamem/bam/$(sample).bwamem.bam)
endif

bwamem : $(BWA_BAMS) $(addsuffix .bai,$(BWA_BAMS))

bam/%.bam : bwamem/bam/%.bwamem.$(BAM_SUFFIX)
	$(INIT) ln -f $< $@

ifeq ($(MERGE_SPLIT_FASTQ),true)
bwamem/bam/%.bwamem.bam : fastq/%.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),fastq/%.2.fastq.gz)
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9\-]\+//'`; \
	$(call RUN,$(BWA_NUM_CORES),$(RESOURCE_REQ_VHIGH_MEM),$(RESOURCE_REQ_MEDIUM),$(BWA_MODULE) $(SAMTOOLS_MODULE),"\
	$(BWA_MEM) -t $(BWA_NUM_CORES) $(BWA_MEM_OPTS) \
	-R \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(REF_FASTA) $^ | \
	$(SAMTOOLS) view -bhS - > $@")
else
define align-split-fastq
bwamem/bam/$1.bwamem.bam : fastq/$1.1.fastq.gz $$(if $$(findstring true,$$(PAIRED_END)),fastq/$1.2.fastq.gz)
	LBID=`echo "$1" | sed 's/_[A-Za-z0-9\-]\+//'`; \
	$$(call RUN,$$(BWA_NUM_CORES),$$(RESOURCE_REQ_VHIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(BWA_MODULE) $$(SAMTOOLS_MODULE),"\
	$$(BWA_MEM) -t $$(BWA_NUM_CORES) $$(BWA_MEM_OPTS) \
	-R \"@RG\tID:$1\tLB:$$$${LBID}\tPL:$${SEQ_PLATFORM}\tSM:$$$${LBID}\" $$(REF_FASTA) $$^ | \
	$$(SAMTOOLS) view -bhS - > $$@")
endef
$(foreach sample,$(SPLIT_SAMPLES),$(foreach split,$(split.$(sample)),$(eval $(call align-split-fastq,$(split)))))

define merged-bwamem-bam
bwamem/bam/$1.bwamem.header.sam : $$(foreach split,$2,bwamem/bam/$$(split).bwamem.bam)
	$$(INIT) module load $$(SAMTOOLS_MODULE); $$(SAMTOOLS) view -H $$< | grep -v '^@RG' > $$@.tmp; \
	for bam in $$^; do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $$(RM) $$@.tmp

bwamem/bam/$1.bwamem.bam : bwamem/bam/$1.bwamem.header.sam $$(foreach split,$2,bwamem/bam/$$(split).bwamem.bam) 
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SAMTOOLS_MODULE),"\
	$$(SAMTOOLS) merge -f -h $$< $$(@) $$(filter %.bam,$$^) && $$(RM) $$<")
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-bwamem-bam,$(sample),$(split.$(sample)))))


endif

include usb-modules-v2/fastq_tools/fastq.mk
include usb-modules-v2/bam_tools/processBam.mk
include usb-modules-v2/aligners/align.mk
include usb-modules-v2/variant_callers/gatk.mk
