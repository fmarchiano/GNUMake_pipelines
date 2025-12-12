include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/mutsigcv.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : mutsigcv

mutsigcv : siggenes/all$(PROJECT_PREFIX).mutsigcv.sig_genes.txt

siggenes/all$(PROJECT_PREFIX).mutsigcv.sig_genes.txt : siggenes/all$(PROJECT_PREFIX).mutsigcv_input.maf
	$(call RUN,1,$(RESOURCE_REQ_VHIGH_MEM),$(RESOURCE_REQ_SHORT),,"$(MUTSIGCV) $(MCR) $^ \
	$(MUTSIGCV_COVERAGE_REF) $(MUTSIGCV_COV_REF) all$(PROJECT_PREFIX).mutsigcv $(MUTSIGCV_DICT_REF) $(MUTSIGCV_SEQ_REF_DIR) && \
	mv all$(PROJECT_PREFIX).mutsigcv.* siggenes")

siggenes/all$(PROJECT_PREFIX).mutsigcv_input.maf : $(if $(findstring EXCEL,$(MUTATION_SUMMARY_FORMAT)),summary/mutation_summary$(PROJECT_PREFIX).$(DOWNSTREAM_EFF_TYPES).xlsx,summary/mutation_summary$(PROJECT_PREFIX).$(DOWNSTREAM_EFF_TYPES).txt)
	$(call RUN,1,$(RESOURCE_REQ_VHIGH_MEM),$(RESOURCE_REQ_SHORT),$(R_MODULE),"\
	$(MUTSIGCV_MAKE_INPUT) --outFile $@ $<")
