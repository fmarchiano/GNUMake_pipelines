# Merges overlapping paired-end reads using ClipOverlap from bamUtil (https://genome.sph.umich.edu/wiki/BamUtil:_clipOverlap)
# When two mates overlap, this tool will clip the record's whose clipped region would have the lowest average quality.
# This step is needed if running Strelka/Strelka2 to avoid double-counting overlapping pairs.

include usb-modules-v2/Makefile.inc
LOGDIR ?= log/bam_clipoverlap.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : bam_clipoverlap

ifeq ($(PAIRED_END),true)
bam_clipoverlap : $(foreach sample,$(SAMPLES),bam_clipoverlap/$(sample).bam)

bam_clipoverlap/%.bam : bam/%.bam
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_SHORT),$(SAMTOOLS_MODULE),"\
	$(BAMUTIL) clipOverlap --poolSize 100000000 --in $< --out $@ &&samtools index $@")
else
bam_clipoverlap:
	$(info *******************************************************************)
	$(info *********** Did nothing because PAIRED_END was not true ***********)
	$(info *******************************************************************)
endif
