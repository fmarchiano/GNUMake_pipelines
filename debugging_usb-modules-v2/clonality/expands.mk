# run expands for determining tumor ploidy

include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/expands.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: expands

expands : $(foreach pair,$(SAMPLE_PAIRS),expands/output/$(pair).expands.RData)
#expands : $(foreach pair,$(SAMPLE_PAIRS),expands/input/$(pair).mutations.txt expands/input/$(pair).segments.txt)

expands/output/%.expands.RData : expands/input/%.mutations.txt expands/input/%.segments.txt
	$(MKDIR) expands/output; \
	$(call RUN,$(EXPANDS_NUM_CORES),$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(EXPANDS) --mutations $(<) --segs $(<<) --numcores $(EXPANDS_NUM_CORES) --outPrefix $(subst .RData,,$@)")

define expands_make_input
expands/input/$1_$2.mutations.txt : $$(foreach prefix,$$(CALLER_PREFIX),tables/$1_$2.$$(call DOWMSTREAM_VCF_TABLE_SUFFIX,$$(prefix)).txt)
	$$(MKDIR) expands/input; \
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(RBIND) --tumorNormal $$^ > $$@.tmp1 && \
	$$(MUTATION_SUMMARY_RSCRIPT) --outFile $$@.tmp2 --outputFormat TXT $$@.tmp1 && \
	$$(EXPANDS_MAKE_INPUT) --outFile $$@ --type mutations $$@.tmp2 && \
	$$(RM) $$@.tmp1 $$@.tmp2")

expands/input/$1_$2.segments.txt : absolute/step3/reviewed/SEG_MAF/$1_$2.segtab.txt
	$$(MKDIR) expands/input; \
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(EXPANDS_MAKE_INPUT) --outFile $$@ --type cna $$<")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call expands_make_input,$(tumor.$(pair)),$(normal.$(pair)))))

#include usb-modules-v2/clonality/absoluteSeq.mk