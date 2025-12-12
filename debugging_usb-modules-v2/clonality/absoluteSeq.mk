include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/absoluteSeq.$(NOW)

ifndef ABSOLUTE_CALLER_PREFIX
  $(info ABSOLUTE_CALLER_PREFIX is not set)
  $(info  )
  exit:
	val=1 && exit $${val}
endif

ifneq ($(words $(ABSOLUTE_CALLER_PREFIX)),1)
  $(info warning: ABSOLUTE_CALLER_PREFIX contains more than one variant caller)
  $(info  )
endif


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: absolute
#absolute_rdata absolute_reviewed absolute_tables

whoami = $(shell whoami)

absolute : $(if $(findstring true,$(ABSOLUTE_STEP_3)),absolute/step3/reviewed/all$(PROJECT_PREFIX).$(whoami).ABSOLUTE.table.txt,absolute/step2/all$(PROJECT_PREFIX).PP-calls_tab.review.$(whoami).txt)

define absolute_make_mutations
absolute/mutations/$1_$2.mutations.txt : $$(foreach prefix,$$(ABSOLUTE_CALLER_PREFIX),tables/$1_$2.$$(call DOWMSTREAM_VCF_TABLE_SUFFIX,$$(prefix)).txt)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE) ,"\
	$$(SINGULARITY_EXEC) $$(ABSOLUTE_IMG) Rscript $(MODULE_SCRIPTS_DIR)/rbind.R --tumorNormal $$^ > $$@.tmp1 && \
	$$(SINGULARITY_EXEC) $$(ABSOLUTE_IMG) Rscript usb-modules-v2/summary/mutation_summary_excel.v2.R --outFile $$@.tmp2 --outputFormat TXT $$@.tmp1 && \
	$$(SINGULARITY_EXEC) $$(ABSOLUTE_IMG) $$(ABSOLUTE_MAKE_MUTS) --sample `basename $$@ .mutations.txt` --outFile $$@ $$@.tmp2 $$(ABSOLUTE_MAKE_MUTS_OPT) && \
	$$(RM) $$@.tmp1 $$@.tmp2")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call absolute_make_mutations,$(tumor.$(pair)),$(normal.$(pair)))))

absolute/segments/%.seg.txt : facets/cncf/%.cncf.txt
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(SINGULARITY_MODULE) ,"\
	$(SINGULARITY_EXEC) $(ABSOLUTE_IMG) $(ABSOLUTE_MAKE_SEGS) --sample `basename $@ .seg.txt` --outFile $@ $(ABSOLUTE_MAKE_SEGS_OPT) $<")

absolute/absolute_parameters_$(PROJECT_PREFIX).txt : $(foreach pair,$(SAMPLE_PAIRS),absolute/mutations/$(pair).mutations.txt absolute/segments/$(pair).seg.txt)
	$(INIT) echo "include patient sample seg.dat.fn sigma.p max.sigma.h min.ploidy max.ploidy \
	primary.disease platform sample.name results.dir max.as.seg.count copy_num_type \
	max.neg.genome max.non.clonal maf.fn min.mut.af output.fn.base" | perl -p -e "s/\s+/\t/g;" > $@ && \
	echo "" >> $@ && \
	for sample_pair in $(SAMPLE_PAIRS); do \
		echo "Y $$sample_pair $$sample_pair absolute/segments/$$sample_pair.seg.txt \
		$(ABSOLUTE_SIGMA_P) $(ABSOLUTE_MAX_SIGMA_H) $(ABSOLUTE_MIN_PLOIDY) $(ABSOLUTE_MAX_PLOIDY) \
		$(ABSOLUTE_DISEASE) $(ABSOLUTE_PLATFORM) $$sample_pair absolute/step1/$$sample_pair \
		$(ABSOLUTE_MAX_SEG) $(ABSOLUTE_COPYNUMTYPE) $(ABSOLUTE_MAX_NEG_GENOME) \
		$(ABSOLUTE_MAX_NON_CLONAL) absolute/mutations/$$sample_pair.mutations.txt $(ABSOLUTE_MIN_MUT_AF) $$sample_pair" | \
		perl -p -e "s/\s+/\t/g;" >> $@ && echo "" >> $@; \
	done

define absolute_step1
absolute/step1/$1_$2/$1_$2.ABSOLUTE.RData : absolute/absolute_parameters_$(PROJECT_PREFIX).txt absolute/mutations/$1_$2.mutations.txt absolute/segments/$1_$2.seg.txt 
	$$(call RUN,$$(ABSOLUTE_NUM_CORE),$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(SINGULARITY_MODULE),"\
	$$(SINGULARITY_EXEC) $$(ABSOLUTE_IMG) $$(ABSOLUTE_STEP1) --params $$< --numCores $$(ABSOLUTE_NUM_CORE) --sample `basename $$@ .ABSOLUTE.RData`")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call absolute_step1,$(tumor.$(pair)),$(normal.$(pair)))))

absolute/step2/all$(PROJECT_PREFIX).PP-modes.data.RData : $(foreach pair,$(SAMPLE_PAIRS),absolute/step1/$(pair)/$(pair).ABSOLUTE.RData)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(SINGULARITY_MODULE),"\
	$(SINGULARITY_EXEC) $(ABSOLUTE_IMG) $(ABSOLUTE_STEP2) --outdir $(@D) --obj.name all$(PROJECT_PREFIX) $^")

absolute/step2/all$(PROJECT_PREFIX).PP-calls_tab.txt : absolute/step2/all$(PROJECT_PREFIX).PP-modes.data.RData
	

absolute/step2/all$(PROJECT_PREFIX).PP-calls_tab.review.$(whoami).txt : absolute/step2/all$(PROJECT_PREFIX).PP-calls_tab.txt
	sed -e "s/^/\t/" -e "s/^\tarray/override\tarray/" $^ > $@

ifeq ($(ABSOLUTE_STEP_3),true)
absolute/step3/reviewed/all$(PROJECT_PREFIX).$(whoami).ABSOLUTE.table.txt : absolute/step2/all$(PROJECT_PREFIX).PP-modes.data.RData absolute/step2/all$(PROJECT_PREFIX).PP-calls_tab.review.$(whoami).txt
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(SINGULARITY_MODULE),"\
	$(SINGULARITY_EXEC) $(ABSOLUTE_IMG) $(ABSOLUTE_STEP3) --obj.name all$(PROJECT_PREFIX) --analyst $(whoami) \
	--outdir absolute/step3 --modes.fn $(<) --pp.calls $(<<)")

define absolute_step3
absolute/step3/reviewed/SEG_MAF/$1_$2.segtab.txt : absolute/step3/reviewed/all$(PROJECT_PREFIX).$(whoami).ABSOLUTE.table.txt
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call absolute_step3,$(tumor.$(pair)),$(normal.$(pair)))))
endif
