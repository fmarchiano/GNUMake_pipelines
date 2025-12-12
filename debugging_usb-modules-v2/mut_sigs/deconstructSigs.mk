include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/deconstructSigs.$(NOW)

ifneq ($(words $(CALLER_PREFIX)),1)
  $(info CALLER_PREFIX contains more than one variant caller)
  $(info Choose only one by executing: make deconstruct_sigs CALLER_PREFIX=<variant caller>)
  $(info  )
  exit:
	val=1 && exit $${val}
endif

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: deconstructSigs

deconstructSigs : deconstructSigs/all$(PROJECT_PREFIX).deconstructSigs.txt

# if bootstrap needs to be done, then do each sample separately
# otherwise just run everything together
ifdef DECONSTRUCTSIGS_NUMITER
deconstructSigs/all$(PROJECT_PREFIX).deconstructSigs.txt : $(foreach pair,$(SAMPLE_PAIRS),deconstructSigs/$(pair).deconstructSigs.txt)
	$(INIT) head -1 $< > $@; \
	for ds in $^; do tail -n +2 $$ds >> $@; done
else
deconstructSigs/all$(PROJECT_PREFIX).deconstructSigs.txt : deconstructSigs/all$(PROJECT_PREFIX).deconstructSigs.RData
	
endif

deconstructSigs/all$(PROJECT_PREFIX).deconstructSigs.RData : summary/mutation_summary$(PROJECT_PREFIX).$(DOWNSTREAM_EFF_TYPES).detected.txt
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
		$(DECONSTRUCTSIGS) --outPrefix=$(subst .RData,,$@) $<")

define deconstruct-sigs
deconstructSigs/$1_$2.deconstructSigs.RData : $$(foreach prefix,$$(CALLER_PREFIX),tables/$1_$2.$$(call DOWMSTREAM_VCF_TABLE_SUFFIX,$$(prefix)).txt)
	$$(call RUN,$$(DECONSTRUCTSIGS_NUMCORES),$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(RBIND) --tumorNormal $$^ | grep -v interrogation_Absent > $$@.tmp && \
	$$(DECONSTRUCTSIGS) --outPrefix $$(subst .RData,,$$@) \
	--num_iter $$(DECONSTRUCTSIGS_NUMITER) --num_cores $$(DECONSTRUCTSIGS_NUMCORES) \
	$$(if $$(DECONSTRUCTSIGS_ASSOC_SIGS),--associated $$(DECONSTRUCTSIGS_ASSOC_SIGS)) \
	$$@.tmp && \
	$$(RM) $$@.tmp")

deconstructSigs/$1_$2.deconstructSigs.txt : deconstructSigs/$1_$2.deconstructSigs.RData
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call deconstruct-sigs,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/summary/mutationSummary.mk
