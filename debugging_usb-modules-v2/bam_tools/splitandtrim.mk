# various bam processing steps
##### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc

LOGDIR ?= log/splitandtrim.$(NOW)

BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
split_and_trim : $(BAMS) $(addsuffix .bai,$(BAMS))

bam/%.bam : unprocessed_bam/%.splitntrim.bam
	$(INIT) ln -f $(<) $(@)


include usb-modules-v2/bam_tools/processBam.mk
