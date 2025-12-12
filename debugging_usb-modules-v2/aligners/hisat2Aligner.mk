# This module is for the hisat2 aligner
include usb-modules-v2/Makefile.inc
include usb-modules-v2/aligners/align.inc

ALIGNER := hisat2

HISAT2_NUM_CORES ?= 8

ifeq ($(STRAND_SPECIFICITY),FIRST_READ_TRANSCRIPTION_STRAND)
HISAT2_OPTS += --rna-strandness FR
endif
ifeq ($(STRAND_SPECIFICITY),SECOND_READ_TRANSCRIPTION_STRAND)
HISAT2_OPTS += --rna-strandness RF
endif


LOGDIR ?= log/hisat2.$(NOW)

..DUMMY := $(shell mkdir -p version; $(HISAT2) --version &> version/hisat2.txt; echo "options: $(HISAT2_OPTS)" >> version/hisat2.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : hisat2_bams
	
ifeq ($(strip $(PRIMARY_ALIGNER)),hisat2)
HISAT2_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
else
HISAT2_BAMS = $(foreach sample,$(SAMPLES),hisat2/bam/$(sample).hisat2.bam)
endif

hisat2_bams : $(HISAT2_BAMS) $(addsuffix .bai,$(HISAT2_BAMS))

bam/%.bam : hisat2/bam/%.hisat2.$(BAM_SUFFIX)
	$(INIT) ln -f $(<) $(@) 

hisat2/bam/%.raw.bam hisat2/unmapped_bam/%.raw_unmapped.bam : fastq/%.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),fastq/%.2.fastq.gz)
	$(call RUN,$(HISAT2_NUM_CORES),$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(HISAT2_MODULE) $(SAMTOOLS_MODULE),"\
	$(MKDIR) hisat2/bam hisat2/unmapped_bam; \
	LBID=`echo \"$*\" | sed 's/_.\+//'`; \
	$(HISAT2) $(HISAT2_OPTS) -x $(HISAT2_REF) -p $(HISAT2_NUM_CORES) \
	--rg-id $* --rg \"SM:\$${LBID}\" \
	--rg PL:${SEQ_PLATFORM} --rg \"LB:\$${LBID}\" \
	-S hisat2/bam/$*.raw.sam --un hisat2/unmapped_bam/$*.raw_unmapped.sam \
	-1 $(<) $(if $(findstring true,$(PAIRED_END)),-2 $(<<)) && \
	$(SAMTOOLS) view -Sbh hisat2/bam/$*.raw.sam > hisat2/bam/$*.raw.bam && \
	$(SAMTOOLS) view -Sbh hisat2/unmapped_bam/$*.raw_unmapped.sam > hisat2/unmapped_bam/$*.raw_unmapped.bam && \
	$(RM) hisat2/bam/$*.raw.sam hisat2/unmapped_bam/$*.raw_unmapped.sam")

hisat2/bam/%.hisat2.bam : hisat2/bam/%.raw.reordered.sorted.bam
	$(INIT) mv $< $@
	
hisat2/unmapped_bam/%.hisat2_unmapped.bam : hisat2/unmapped_bam/%.raw_unmapped.reordered.sorted.bam
	$(INIT) mv $< $@
	
include usb-modules-v2/fastq_tools/fastq.mk
include usb-modules-v2/bam_tools/processBam.mk
include usb-modules-v2/aligners/align.mk
