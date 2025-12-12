include usb-modules-v2/Makefile.inc

LOGDIR ?= log/lst.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: lst

lst : lst/all$(PROJECT_PREFIX).lst.txt

lst/all$(PROJECT_PREFIX).lst.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(LST) --centromere_file $(CENTROMERE_TABLE) --outFile $@ $^")

#lst/all$(PROJECT_PREFIX).lst_class.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.out) lst/all$(PROJECT_PREFIX).lst.txt
#	$(INIT) \
	
