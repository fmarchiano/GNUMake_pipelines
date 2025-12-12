include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/mutational_patterns.$(NOW)

ifneq ($(words $(CALLER_PREFIX)),1)
  $(info CALLER_PREFIX contains more than one variant caller)
  $(info Choose only one by executing: make mutational_patterns CALLER_PREFIX=<variant caller>)
  $(info  )
  exit:
	val=1 && exit $${val}
endif

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: mutational_patterns

mutational_patterns : mutational_patterns/all$(PROJECT_PREFIX).$(MUT_SIG_COSMIC).RData

mutational_patterns/all$(PROJECT_PREFIX).$(MUT_SIG_COSMIC).RData : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(CALLER_PREFIX).*.hotspot.pass.vcf)
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_SHORT),$(R_MODULE),"\
	$(MUTATIONALPATTERNS) --outPrefix mutational_patterns/all$(PROJECT_PREFIX) $^")

