include usb-modules-v2/Makefile.inc
include usb-modules-v2/config.inc

LOGDIR ?= log/poolednorm_bam.$(NOW)

poolednorm_bam : bam/poolednorm.bam $(addsuffix .bai,bam/poolednorm.bam)

unprocessed_bam/poolednorm.bam : $(foreach normal,$(POOLED_NORM_SAMPLES),bam/$(normal).downsampled.bam) $(addsuffix .bai,$(foreach normal,$(POOLED_NORM_SAMPLES),bam/$(normal).downsampled.bam))
	$(MKDIR) unprocessed_bam; \
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_MEDIUM),$(SAMTOOLS_MODULE),"\
	$(SAMTOOLS) merge -f $@ $(filter %.bam,$^) && $(RM) $^")

include usb-modules-v2/bam_tools/processBam.mk
include usb-modules-v2/bam_tools/fixRG.mk



