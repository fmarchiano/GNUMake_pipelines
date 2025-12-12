include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/sig_profiler_assignment.$(NOW)

ifneq ($(words $(CALLER_PREFIX)),1)
  $(info CALLER_PREFIX contains more than one variant caller)
  $(info Choose only one by executing: make sig_profiler_assignment CALLER_PREFIX=<variant caller>)
  $(info  )
  exit:
	val=1 && exit $${val}
endif

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: sig_profiler_assignment

sig_profiler_assignment : SigProfilerAssignment/JOB_METADATA_SPA.txt

#SigProfilerAssignment/input/$(foreach pair,$(SAMPLE_PAIRS),$(pair).$(CALLER_PREFIX).*.hotspot.pass.vcf) : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(CALLER_PREFIX).*.hotspot.pass.vcf)
#	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
#	$(RM) SigProfilerAssignment/input/* &&\
#	ln $^ SigProfilerAssignment/input")

SigProfilerAssignment/JOB_METADATA_SPA.txt : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(CALLER_PREFIX).*.hotspot.pass.vcf)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_MEDIUM),$(SIG_PROFILER_MODULE),"\
	$(RM) -r SigProfilerAssignment/* && mkdir SigProfilerAssignment/input && \
	ln $^ SigProfilerAssignment/input && \
	$(SIG_PROFILER_ASSIGNMENT) --samples SigProfilerAssignment/input --output SigProfilerAssignment \
	--cosmic_version $(SIG_PROFILER_COSMIC_VERSION) \
	$(if $(findstring True,$(SIG_PROFILER_COSMIC_EXOME)),--exome) \
	--genome_build $(SIG_PROFILER_COSMIC_GENOME) \
	$(if $(SIG_PROFILER_COSMIC_SIGNATURE_DB),--signature_database $(SIG_PROFILER_COSMIC_SIGNATURE_DB)) \
	$(if $(SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS),--signature_database $(SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS))")
