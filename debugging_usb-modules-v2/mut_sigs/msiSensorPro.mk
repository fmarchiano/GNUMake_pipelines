include usb-modules-v2/Makefile.inc

LOGDIR ?= log/msisensorpro.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: msisensorpro

msisensorpro : msisensorpro/all$(PROJECT_PREFIX).msisensorpro.txt

msisensorpro/all$(PROJECT_PREFIX).msisensorpro.txt : $(foreach pair,$(SAMPLE_PAIRS),msisensorpro/$(pair).msisensorpro)
	$(INIT) head -1 $< | sed 's/^/TUMOR_NORMAL\t/' > $@; \
	for msi in $^; do \
		samplename=`basename $$msi | sed 's/\.msisensorpro//;'`; \
		sed '/^Total_Number/d;' $$msi | sed "s/^/$$samplename\t/" >> $@; \
	done

define msisensorpro-msi
msisensorpro/$1_$2.msisensorpro : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call RUN,4,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_MEDIUM),,"\
		export LD_LIBRARY_PATH=/home/ng_piscuoglio/pipeline/usr/local/htslib-1.19.1:$LD_LIBRARY_PATH; \
		$$(MSISENSORPRO) msi -d $$(MSISENSORPRO_REF) -b 4 -t $$(<) -n $$(<<) -o $$(subst .complete,,$$@) \
		$$(if $$(TARGETS_FILE_INTERVALS),-e $$(TARGETS_FILE_INTERVALS)) || exit 0; \
		$$(RM) $$(@)_dis $$(@)_unstable $$(@)_all;")
#		if grep \"Total time consumed\" $$@; then exit 0; else exit 1; fi")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call msisensorpro-msi,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/bam_tools/processBam.mk
