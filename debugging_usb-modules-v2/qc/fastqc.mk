# vim: set ft=make :
# Run Fastqc on bam files

include usb-modules-v2/Makefile.inc

LOGDIR ?= log/fastqc.$(NOW)

.PHONY: fastqc
.SECONDARY: 

#fastqc : $(foreach sample,$(SAMPLES),fastqc/$(sample)_fastqc/summary.txt) fastqc/all_summary.txt
fastqc : $(foreach sample,$(SAMPLES),fastqc/$(sample)_fastqc.zip)
fastqc/%_fastqc.zip : bam/%.bam
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(FASTQC_MODULE),"\
	$(FASTQC) -o fastqc $^")

fastqc/%_fastqc/summary.txt : fastqc/%_fastqc.zip
	$(INIT) $(UNZIP) -o -d fastqc $< &> $(LOG) && touch $@

fastqc/all_summary.txt : $(foreach sample,$(SAMPLES),fastqc/$(sample)_fastqc/summary.txt)
	$(INIT) module load $(R_MODULE); $(FASTQC_SUMMARY_PLOT) --outPrefix fastqc/all_summary $^ &> $(LOG)
