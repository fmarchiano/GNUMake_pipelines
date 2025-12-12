include usb-modules-v2/Makefile.inc

STRINGTIE_NUM_CORES ?= 4


LOGDIR ?= log/ballgown.$(NOW)

.PHONY: ballgown
.DELETE_ON_ERROR:

ballgown : $(foreach sample,$(SAMPLES),stringtie/$(sample).ballgown.gtf)

stringtie/%.stringtie.gtf : hisat2/bam/%.hisat2.bam
	$(call RUN,$(STRINGTIE_NUM_CORES),$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(STRINGTIE_MODULE),"\
	$(STRINGTIE) -p $(STRINGTIE_NUM_CORES) -G $(GENCODE_GTF) -o $@ -l $* $<")

stringtie/all.stringtie_merged.gtf : $(foreach sample,$(SAMPLES),stringtie/$(sample).stringtie.gtf)
	ls $^ > stringtie/all.stringtie_merge_list.txt; \
	$(call RUN,$(STRINGTIE_NUM_CORES),$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(STRINGTIE_MODULE),"\
	$(STRINGTIE) --merge -p $(STRINGTIE_NUM_CORES) -G $(GENCODE_GTF) -o $@ stringtie/all.stringtie_merge_list.txt && \
	$(RM) stringtie/all.stringtie_merge_list.txt")

stringtie/%.ballgown.gtf : stringtie/all.stringtie_merged.gtf hisat2/bam/%.hisat2.bam
	$(call RUN,$(STRINGTIE_NUM_CORES),$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(STRINGTIE_MODULE),"\
	$(STRINGTIE) -p $(STRINGTIE_NUM_CORES) -B -G $< $(<<) > $@")

include usb-modules-v2/aligners/hisat2Aligner.mk
include usb-modules-v2/bam_tools/processBam.mk
