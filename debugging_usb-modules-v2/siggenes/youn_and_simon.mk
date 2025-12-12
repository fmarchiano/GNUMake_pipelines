include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/youn_and_simon.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : youn_and_simon

youn_and_simon : siggenes/all$(PROJECT_PREFIX).youn_and_simon.sig_genes.txt

siggenes/all$(PROJECT_PREFIX).youn_and_simon.sig_genes.txt : siggenes/all$(PROJECT_PREFIX).youn_and_simon_input.nonsilent.maf siggenes/all$(PROJECT_PREFIX).youn_and_simon_input.silent.maf
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(R_MODULE),"\
 	$(YOUN_AND_SIMON) --outFile $@ --nonsilentFile $(<) --silentFile $(<<) \
	--sequenceDataFile $(YOUN_AND_SIMON_SEQ_DATA) --numCases $(words $(SAMPLE_PAIRS))")

siggenes/all$(PROJECT_PREFIX).youn_and_simon_input.silent.maf : $(if $(findstring EXCEL,$(MUTATION_SUMMARY_FORMAT)),summary/mutation_summary$(PROJECT_PREFIX).$(DOWNSTREAM_EFF_TYPES).xlsx,summary/mutation_summary$(PROJECT_PREFIX).$(DOWNSTREAM_EFF_TYPES).txt)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(YOUN_AND_SIMON_MAKE_INPUT) --outPrefix $(@D)/all$(PROJECT_PREFIX).youn_and_simon_input $<")

siggenes/all$(PROJECT_PREFIX).youn_and_simon_input.nonsilent.maf : siggenes/all$(PROJECT_PREFIX).youn_and_simon_input.silent.maf
	
