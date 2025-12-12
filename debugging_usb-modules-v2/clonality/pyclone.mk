include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/pyclone.$(NOW)

ifneq ($(words $(CALLER_PREFIX)),1)
  $(info CALLER_PREFIX contains more than one variant caller)
  $(info Choose only one by executing: make pyclone CALLER_PREFIX=<variant caller>)
  $(info  )
  exit:
	val=1 && exit $${val}
endif

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : pyclone

pyclone : $(foreach normal_sample,$(NORMAL_SAMPLES),pyclone/tables/$(normal_sample).run2.clusters.txt \
pyclone/tables/$(normal_sample).run1.clusters.signatures.RData pyclone/tables/$(normal_sample).run2.clusters.signatures.RData)

# There are edge cases with no mutations...
CHECK_PYCLONE_CONFIG = if [ `grep tumour_content $1 | wc -l` -eq 0 ] ; then touch $2; else $3; fi

# not used yet, and probably does not take care of empty files
pyclone/all$(PROJECT_PREFIX).num_clust.txt : $(foreach normal_sample,$(NORMAL_SAMPLES),pyclone/tables/$(normal_sample).cluster.txt)
	$(INIT) \
	for cl in $^; do awk '$3>1&&$4>0.05' $cl | cut -f1 | sort | uniq -c >>$@;
	done
	
pyclone/tables/%.clusters.txt : pyclone/tables/%.loci.txt
	$(INIT) \
	if [ `wc -l $< | cut -f1 -d' '` -gt 0 ]; then \
		ml $(R_MODULE);\
		$(PYCLONE_MAKE_CLUSTER_TABLE) --outFile $@ $^; \
	else \
		touch $@; \
	fi

pyclone/tables/%.clusters.signatures.RData : pyclone/tables/%.loci.txt pyclone/tables/%.clusters.txt
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(PYCLONE_DECONSTRUCTSIGS) --outPrefix $(subst .RData,,$@) \
	--num_iter $(DECONSTRUCTSIGS_NUMITER) --num_cores $(DECONSTRUCTSIGS_NUMCORES) $<")

pyclone/tables/%.loci.txt : pyclone/configs/%.yaml pyclone/runs/%/alpha.tsv.bz2
	$(call CHECK_PYCLONE_CONFIG,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_LONG),$(PYCLONE_MODULE),"\
	$(PYCLONE) build_table --config_file $< --table_type loci \
	--out_file $(subst cluster,loci,$@) --burnin $(PYCLONE_BURNIN)"))

pyclone/runs/%/alpha.tsv.bz2 : pyclone/configs/%.yaml
	$(MKDIR) $(@D); \
	$(call CHECK_PYCLONE_CONFIG,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_LONG),$(PYCLONE_MODULE),"\
	$(PYCLONE) run_analysis --config_file $< --seed $(PYCLONE_SEED)"))

define pyclone_make_config
pyclone/configs/$1.run$3.yaml : $$(foreach tumor,$2,facets/cncf/$$(tumor)_$1.out pyclone/mutations/$$(tumor)_$1.run$3.mutations.yaml)
	$$(INIT) $$(MKDIR) pyclone/configs; \
	echo -n "num_iters: " > $$@; \
	echo $$(PYCLONE_ITER) >> $$@; \
	echo "base_measure_params:" >> $$@; \
	echo "  alpha: 1" >> $$@; \
	echo "  beta: 1" >> $$@; \
	echo "concentration:" >> $$@; \
	echo "  value: 1.0" >> $$@; \
	echo "  prior:" >> $$@; \
	echo "    shape: 1.0" >> $$@; \
	echo "    rate: 0.001" >> $$@; \
	echo "density: pyclone_beta_binomial" >> $$@; \
	echo "beta_binomial_precision_params:" >> $$@; \
	echo "  value: 1000" >> $$@; \
	echo "  prior:" >> $$@; \
	echo "    shape: 1.0" >> $$@; \
	echo "    rate: 0.0001" >> $$@; \
	echo "  proposal:" >> $$@; \
	echo "    precision: 0.01" >> $$@; \
	echo -n "working_dir: " >> $$@; echo `pwd` >> $$@; \
	echo -n "trace_dir: " >> $$@; echo "pyclone/runs/$1.run$3" >> $$@; \
	echo "samples:" >> $$@; \
	for cncf in $$(filter %.out,$$^); do \
		samplename=`basename $$$$cncf | cut -f1 -d'_'`; \
		tnname=`basename $$$$cncf | cut -f1 -d'.'`; \
		if [ `wc -l pyclone/mutations/$$$${tnname}.run$3.mutations.yaml | cut -f1 -d' '` -gt 1 ]; then \
			echo "  $$$$samplename:" >> $$@; \
			echo "    tumour_content: " >> $$@; \
			echo -n "      value: " >> $$@; \
			echo `grep Purity $$$$cncf | cut -f2 -d'=' | sed 's/NA/0.1/;'` >> $$@; \
			echo -n "    mutations_file: pyclone/mutations/" >> $$@; \
			echo "$$$${tnname}.run$3.mutations.yaml" >> $$@; \
			echo "    error_rate: 0.001"  >> $$@; \
		fi; \
	done
endef
$(foreach set,$(SAMPLE_SETS),$(foreach run,1 2,\
	$(eval $(call pyclone_make_config,$(lastword $(subst _,$(space),$(set))),\
	$(wordlist 1,$(shell expr $(words $(subst _,$(space),$(set))) - 1),$(subst _,$(space),$(set))),$(run)))))

pyclone/mutations/%.mutations.yaml : pyclone/mutations/%.mutations.txt
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(PYCLONE_MODULE),"\
	$(PYCLONE) build_mutations_file --prior $(PYCLONE_PRIOR) --in_file $< --out_file $@")

#define pyclone_make_mutations
#pyclone/mutations/$1_$2.mutations.txt : $$(foreach prefix,$$(CALLER_PREFIX),tables/$1_$2.$$(call DOWMSTREAM_VCF_TABLE_SUFFIX,$$(prefix)).txt)
#	$$(MKDIR) pyclone/mutations; \
#	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
#	$$(RBIND) --tumorNormal $$^ > $$@.tmp1 && \
#	$$(MUTATION_SUMMARY_RSCRIPT) --outFile $$@.tmp2 --outputFormat TXT $$@.tmp1 && \
#	$$(PYCLONE_MAKE_MUT_TXT) --outFile $$@ $$@.tmp2 && \
#	$$(RM) $$@.tmp1 $$@.tmp2")
#endef
#$(foreach pair,$(SAMPLE_PAIRS),\
#	$(eval $(call pyclone_make_mutations,$(tumor.$(pair)),$(normal.$(pair)))))
	
pyclone/mutations/%.run1.mutations.txt : $(foreach prefix,$(CALLER_PREFIX),tables/%.$(call DOWMSTREAM_VCF_TABLE_SUFFIX,$(prefix)).txt)
	$(MKDIR) pyclone/mutations; \
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(RBIND) --tumorNormal $^ > $@.tmp1 && \
	$(MUTATION_SUMMARY_RSCRIPT) --outFile $@.tmp2 --outputFormat TXT $@.tmp1 && \
	$(PYCLONE_MAKE_MUT_TXT) --outFile $@ $@.tmp2 && \
	$(RM) $@.tmp1 $@.tmp2")
	
define pyclone_make_mutations_run2	
pyclone/mutations/$1_$2.run2.mutations.txt : pyclone/mutations/$1_$2.run1.mutations.txt pyclone/tables/$2.run1.loci.txt
	$$(MKDIR) pyclone/mutations; \
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(PYCLONE_POSTPY) --maxSD $$(PYCLONE_POSTPY_MAXSD) --outFile $$@ $$^")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call pyclone_make_mutations_run2,$(tumor.$(pair)),$(normal.$(pair)))))

	
	
	
