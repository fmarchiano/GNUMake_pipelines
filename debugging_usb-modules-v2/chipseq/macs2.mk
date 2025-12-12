include usb-modules-v2/Makefile.inc

LOGDIR ?= log/macs2.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : macs2

# if there is an input sample, list it as "normal" in sample_sets.txt
ifdef SAMPLE_PAIRS
macs2 : $(foreach pair,$(SAMPLE_PAIRS),macs2/$(pair)_summits.bed)

define macs2-callpeaks
macs2/$1_$2_summits.bed : bam/$1.bam bam/$2.bam
	$$(call RUN,$$(MACS2_NUM_CORES),$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_SHORT),$$(MACS2_MODULE),"\
	$$(MACS2) callpeak -t $$(<) -c $$(<<) \
	-f $(if $$(findstring true,$$(PAIRED_END)),BAMPE,BAM) \
	-g hs -p $$(MACS2_PVALUE) -n $1_$2 --outdir $$(@D)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call macs2-callpeaks,$(tumor.$(pair)),$(normal.$(pair)))))
endif

# if there is no input sample, process each sample individually
ifndef SAMPLE_PAIRS
macs2 : $(foreach sample,$(SAMPLES),macs2/$(sample)_summits.bed)

define macs2-callpeaks
macs2/$1_summits.bed : bam/$1.bam
	$$(call RUN,$$(MACS2_NUM_CORES),$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_SHORT),$$(MACS2_MODULE),"\
	$$(MACS2) callpeak -t $$(<) \
	-f $(if $$(findstring true,$$(PAIRED_END)),BAMPE,BAM) \
	-g hs -p $$(MACS2_PVALUE) -n $1 --outdir $$(@D)")
endef
$(foreach sample,$(SAMPLES),$(eval $(call macs2-callpeaks,$(sample))))
endif


macs2/%_peaks.narrowPeak : macs2/%_summits.bed
	

%.homer_ann.txt : %
	cut -f1-5 $< | sed 's/^/chr/' - > %.bed && |
