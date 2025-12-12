# various bam processing steps
##### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc

LOGDIR ?= log/fixRG.$(NOW)

BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
fixed_bams : $(BAMS) $(addsuffix .bai,$(BAMS))

bam/%.bam : unprocessed_bam/%.rg.bam
	$(INIT) ln -f $(<) $(@)


include usb-modules-v2/bam_tools/processBam.mk
