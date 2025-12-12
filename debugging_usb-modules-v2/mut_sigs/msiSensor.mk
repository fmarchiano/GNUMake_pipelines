include usb-modules-v2/Makefile.inc

LOGDIR ?= log/msisensor.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: msisensor

msisensor : msisensor/all$(PROJECT_PREFIX).msisensor.txt

msisensor/all$(PROJECT_PREFIX).msisensor.txt : $(foreach pair,$(SAMPLE_PAIRS),msisensor/$(pair).msisensor)
	$(INIT) head -1 $< | sed 's/^/TUMOR_NORMAL\t/' > $@; \
	for msi in $^; do \
		samplename=`basename $$msi | sed 's/\.msisensor//;'`; \
		sed '/^Total_Number/d;' $$msi | sed "s/^/$$samplename\t/" >> $@; \
	done

define msisensor-msi
msisensor/$1_$2.msisensor : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call RUN,4,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_MEDIUM),,"\
		$$(MSISENSOR) msi -d $$(MSISENSOR_REF) -b 4 -t $$(<) -n $$(<<) -o $$(subst .complete,,$$@) \
		$$(if $$(TARGETS_FILE_INTERVALS),-e $$(TARGETS_FILE_INTERVALS)) || exit 0; \
		$$(RM) $$(@)_dis $$(@)_germline $$(@)_somatic; \
		if grep \"Total time consumed\" $$@; then exit 0; else exit 1; fi")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call msisensor-msi,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/bam_tools/processBam.mk
