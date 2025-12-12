include usb-modules-v2/Makefile.inc

LOGDIR ?= log/downsampleBam.$(NOW)

BAM = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
downsample_bam : $(BAM) $(addsuffix .bai,$(BAMS))

bam/%.bam : unprocessed_bam/%.downsampled.bam
	$(INIT) ln -f $(<) $(@)

include usb-modules-v2/bam_tools/processBam.mk
