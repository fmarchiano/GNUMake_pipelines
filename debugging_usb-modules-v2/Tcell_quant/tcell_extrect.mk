include usb-modules-v2/Makefile.inc

LOGDIR ?= log/tcell_extrect.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : tcell_extrect

tcell_extrect : tcell_extrect/all.resTcellExTRECT.txt

define extrect
tcell_extrect/$1.resTcellExTRECT.txt : bam/$1.bam facets/cncf/$1_$2.out facets/cncf/$1_$2.cncf.txt
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(R4_MODULE) $$(SAMTOOLS_MODULE),"\
	$$(TCELLEXTRECT) --sample_id $1 \
	--genome_build $$(REF) \
	--target_bed $$(TARGETS_FILE_INTERVALS) \
	--purity `grep Purity $$(filter %.out,$$^) | cut -f2 -d'=' | tr -d ' ' | sed 's/NA/0.1/'` \
	--cncf $$(filter %.cncf.txt,$$^) \
	--outdir tcell_extrect \
	$$(filter %.bam,$$^)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call extrect,$(tumor.$(pair)),$(normal.$(pair)))))

tcell_extrect/all.resTcellExTRECT.txt : $(foreach pair,$(SAMPLE_PAIRS),tcell_extrect/$(tumor.$(pair)).resTcellExTRECT.txt) 
	$(INIT) \
	cat $^ | sed '1p;/sample/d' > $@
