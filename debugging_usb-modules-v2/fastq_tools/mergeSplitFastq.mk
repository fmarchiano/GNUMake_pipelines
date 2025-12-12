# merge split fastqs for workflows like defuse

include usb-modules-v2/Makefile.inc

LOGDIR ?= log/merge_split_fastq.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : fastq

fastq: $(foreach sample,$(SAMPLES),fastq/$(sample).1.fastq.gz fastq/$(sample).2.fastq.gz)

define merged-fastq
fastq/$1.%.fastq.gz : $$(foreach split,$2,fastq/$$(split).%.fastq.gz)
	$$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),,"cat $$(^) > $$(@)")
endef
$(foreach sample,$(SAMPLES),$(eval $(call merged-fastq,$(sample),$(split.$(sample)))))

define merged-fastq2
fastq/$1.1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),,"cat $$(^) > $$(@)")
fastq/$1.2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),,"cat $$(^) > $$(@)")
endef
$(foreach sample,$(SAMPLES),$(eval $(call merged-fastq2,$(sample),$(split.$(sample)))))

