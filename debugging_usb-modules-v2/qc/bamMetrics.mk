# generate bam interval metrics per sample

include usb-modules-v2/Makefile.inc

LOGDIR ?= log/metrics.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: bam_metrics #hs_metrics amplicon_metrics wgs_metrics rna_metrics #interval_report #non_ref_metrics

ifeq ($(CAPTURE_METHOD),NONE)
bam_metrics : wgs_metrics oxog_wgs flagstats alignment_summary_metrics dup artifacts_wgs
endif
ifeq ($(CAPTURE_METHOD),BAITS)
bam_metrics : hs_metrics oxog flagstats alignment_summary_metrics artifacts
endif
ifeq ($(CAPTURE_METHOD),PCR)
bam_metrics : amplicon_metrics flagstats alignment_summary_metrics per_base_depth
endif
ifeq ($(CAPTURE_METHOD),RNA)
bam_metrics : rna_metrics flagstats alignment_summary_metrics
#oxog_wgs flagstats alignment_summary_metrics dup
endif
ifeq ($(CAPTURE_METHOD),CHIP)
#bam_metrics : flagstats alignment_summary
bam_metrics : flagstats
endif

hs_metrics : metrics/all$(PROJECT_PREFIX).hs_metrics.txt metrics/all$(PROJECT_PREFIX).interval_hs_metrics.txt.gz
amplicon_metrics : metrics/all$(PROJECT_PREFIX).amplicon_metrics.txt metrics/all$(PROJECT_PREFIX).interval_amplicon_metrics.txt.gz
per_base_depth : $(foreach sample,$(SAMPLES),metrics/$(sample).per_base_depth.txt.gz)
wgs_metrics : metrics/all$(PROJECT_PREFIX).wgs_metrics.txt
rna_metrics : metrics/all$(PROJECT_PREFIX).rnaseq_metrics.txt metrics/all$(PROJECT_PREFIX).normalized_coverage.rnaseq_metrics.txt
#metrics/all.rnaseq_report/index.html
flagstats : metrics/all$(PROJECT_PREFIX).flagstats.txt
flagstatsQ30 : metrics/all$(PROJECT_PREFIX).flagstatsQ30.txt
alignment_summary_metrics : metrics/all$(PROJECT_PREFIX).alignment_summary_metrics.txt
#gc : $(foreach sample,$(SAMPLES),metrics/$(sample).gc_bias_metrics.txt)
artifacts : $(foreach sample,$(SAMPLES),metrics/$(sample).artifact_metrics.error_summary_metrics)
artifacts_wgs : $(foreach sample,$(SAMPLES),metrics/$(sample).wgs.artifact_metrics.error_summary_metrics)
oxog : metrics/all$(PROJECT_PREFIX).oxog_metrics.txt
oxog_wgs : $(foreach sample,$(SAMPLES),metrics/$(sample).wgs.oxog_metrics.txt)
dup : metrics/all$(PROJECT_PREFIX).dup_metrics.txt

#interval_report : metrics/interval_report/index.html
#non_ref_metrics : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_nonref_freq.txt)

# interval metrics per sample
metrics/%.hs_metrics.txt metrics/%.interval_hs_metrics.txt.gz : bam/%.bam bam/%.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE) $(SAMTOOLS_MODULE),"\
	TMP=`mktemp`.intervals; TMPCOVERED=`mktemp`.covered_intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP &&  grep -P \"\t\" $(TARGETS_FILE_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMPCOVERED &&  grep -P \"\t\" $(TARGETS_FILE_COVERED_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMPCOVERED; \
	$(call PICARD,CollectHsMetrics,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) INPUT=$< OUTPUT=metrics/$*.hs_metrics.txt \
	PER_TARGET_COVERAGE=metrics/$*.interval_hs_metrics.txt TARGET_INTERVALS=\$$TMPCOVERED REFERENCE_SEQUENCE=$(REF_FASTA) \
	BAIT_SET_NAME=hs BAIT_INTERVALS=\$$TMP &&\
	gzip -f metrics/$*.interval_hs_metrics.txt")

metrics/%.amplicon_metrics.txt metrics/%.interval_amplicon_metrics.txt.gz : bam/%.bam bam/%.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE) $(SAMTOOLS_MODULE),"\
	TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP && grep -P \"\t\" $(TARGETS_FILE_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(call PICARD,CollectTargetedPcrMetrics,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) INPUT=$< OUTPUT=$@ REFERENCE_SEQUENCE=$(REF_FASTA)\
	AMPLICON_INTERVALS=\$$TMP TARGET_INTERVALS=\$$TMP \
	PER_TARGET_COVERAGE=metrics/$*.interval_amplicon_metrics.txt COVERAGE_CAP=500000; \
	gzip metrics/$*.interval_amplicon_metrics.txt")

# To avoid crazy RAM usage, sort the target bed to match the chr order in the BAM, and use the '-sorted' option
metrics/%.per_base_depth.txt.gz : bam/%.bam bam/%.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE) $(SAMTOOLS_MODULE),"\
	TMP=`mktemp`.chromosomes; \
	$(SAMTOOLS) idxstats $< | cut -f 1-2 > \$$TMP; \
	$(PURGE_AND_LOAD) $(BEDTOOLS_MODULE); \
	$(BEDTOOLS) sort -faidx \$$TMP -i $(TARGETS_FILE_INTERVALS_MERGED) > \$$TMP.intervals; \
	$(BEDTOOLS) coverage -g \$$TMP -sorted -d -b $< -a \$$TMP.intervals | gzip > $@")
 
define amplicon-metrics-pools
POOLNAME=$$(shell basename $2)
metrics/$1.amplicon_metrics_$$(POOLNAME).txt : bam/$1.bam bam/$1.bam.bai $2
	$$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE) $$(SAMTOOLS_MODULE),"\
	TMP=`mktemp`.intervals; \
	$$(SAMTOOLS) view -H $$< | grep '^@SQ' > \$$$$TMP && grep -P \"\t\" $2 | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$$$1$$(,)\$$$$2+1$$(,)\$$$$3$$(,)\"+\"$$(,)NR }' >> \$$$$TMP; \
	$$(call PICARD,CollectTargetedPcrMetrics,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) INPUT=$$< OUTPUT=$$@ REFERENCE_SEQUENCE=$(REF_FASTA)\
	AMPLICON_INTERVALS=\$$$$TMP TARGET_INTERVALS=\$$$$TMP \
	COVERAGE_CAP=500000")
endef
$(if $(TARGETS_FILE_INTERVALS_POOLS),\
	$(foreach pool,$(TARGETS_FILE_INTERVALS_POOLS),\
		$(foreach sample,$(SAMPLES),\
			$(eval $(call amplicon-metrics-pools,$(sample),$(pool))))))

metrics/%.wgs_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
		$(call PICARD,CollectWgsMetrics,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
		INPUT=$< OUTPUT=$@ COUNT_UNPAIRED=true MINIMUM_MAPPING_QUALITY=30 REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.rnaseq_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
		$(call PICARD,CollectRnaSeqMetrics,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
		REF_FLAT=$(GENE_REF_FLAT) RIBOSOMAL_INTERVALS=$(RIBOSOMAL_INTERVALS) \
		STRAND_SPECIFICITY=$(STRAND_SPECIFICITY) \
		INPUT=$< OUTPUT=$@ VERBOSITY=ERROR REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.alignment_summary_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call PICARD,CollectAlignmentSummaryMetrics,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) INPUT=$< OUTPUT=$@ REFERENCE_SEQUENCE=$(REF_FASTA)")

# does not work, there's a conflict of java versions
#metrics/%.gc_bias_metrics.txt metrics/%.gc_bias_metrics_summary.txt : bam/%.bam bam/%.bam.bai
#	$(call LSCRIPT_MEM,12G,02:59:59,"$(LOAD_JAVA6_MODULE); $(LOAD_R_MODULE); \
#		$(call COLLECT_GCBIAS_METRICS,11G) \
#		INPUT=$< OUTPUT=$@ CHART_OUTPUT=$(addsuffix .pdf,$@) \
#		SUMMARY_OUTPUT=metrics/$*.gc_bias_metrics_summary.txt")

metrics/%.artifact_metrics.error_summary_metrics : bam/%.bam bam/%.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(SAMTOOLS_MODULE) $(JAVA8_MODULE),"\
	TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP &&  grep -P \"\t\" $(TARGETS_FILE_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(call PICARD,CollectSequencingArtifactMetrics,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) INPUT=$< OUTPUT=$(basename $@) \
	DB_SNP=$(DBSNP) INTERVALS=\$$TMP" REFERENCE_SEQUENCE=$(REF_FASTA))

metrics/%.wgs.artifact_metrics.error_summary_metrics : bam/%.bam bam/%.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call PICARD,CollectSequencingArtifactMetrics,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	INPUT=$< OUTPUT=$(basename $@) DB_SNP=$(DBSNP) REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.artifact_metrics.pre_adapter_detail_metrics :	metrics/%.artifact_metrics.error_summary_metrics
	

metrics/%.oxog_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(SAMTOOLS_MODULE) $(JAVA8_MODULE),"\
	TMP=`mktemp`.intervals; \
	$(SAMTOOLS) view -H $< | grep '^@SQ' > \$$TMP &&  grep -P \"\t\" $(TARGETS_FILE_INTERVALS) | \
	awk 'BEGIN {OFS = \"\t\"} { print \$$1$(,)\$$2+1$(,)\$$3$(,)\"+\"$(,)NR }' >> \$$TMP; \
	$(call PICARD,CollectOxoGMetrics,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) INPUT=$< OUTPUT=$@ \
	DB_SNP=$(DBSNP) INTERVALS=\$$TMP REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.wgs.oxog_metrics.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call PICARD,CollectOxoGMetrics,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) INPUT=$< OUTPUT=$@ DB_SNP=$(DBSNP) REFERENCE_SEQUENCE=$(REF_FASTA)")

metrics/%.flagstats.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(SAMTOOLS_MODULE),"\
	$(SAMTOOLS) flagstat $< > $@")

metrics/%.flagstatsQ30.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(SAMTOOLS_MODULE),"\
	$(SAMTOOLS) view -bh -q 30 $< | $(SAMTOOLS) flagstat - > $@")

metrics/%.idxstats.txt : bam/%.bam bam/%.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(SAMTOOLS_MODULE),"\
	$(SAMTOOLS) idxstats $< > $@")

# summarize metrics into one file
metrics/all$(PROJECT_PREFIX).hs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).hs_metrics.txt)
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/SAMPLE.*//; s/BAIT_SET/SAMPLE/; s/\s$$//' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.hs_metrics.txt}); \
		sed "/^#/d; /^BAIT/d; /^\$$/d; s/^hs/$$samplename/; s/\t\+$$//" $$metrics | grep "^$$samplename"; \
	done; \
	} > $@

# summarize interval metrics into one file
# (basename within make will strip extension; system basename will strip dir path, and optional extension)
metrics/all$(PROJECT_PREFIX).interval_hs_metrics.txt.gz : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_hs_metrics.txt.gz)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	zcat $< | sed '/^#/d; /^$$/d' | cut -f 1-6 > $(basename $@).tmp; \
	for metrics in $^; do \
		samplename=\$$(basename \$${metrics} .interval_hs_metrics.txt.gz); \
		zcat \$$metrics | sed '/^#/d; /^$$/d' | cut -f 7-8 | \
		sed \"s/mean_coverage/\$${samplename}_mean_coverage/; s/normalized_coverage/\$${samplename}_normalized_coverage/\" | \
		paste $(basename $@).tmp - > $(basename $@); \
		cp $(basename $@) $(basename $@).tmp; \
	done && gzip $(basename $@) && rm -f $(basename $@).tmp")

# summarize metrics into one file
metrics/all$(PROJECT_PREFIX).amplicon_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).amplicon_metrics.txt) $(if $(TARGETS_FILE_INTERVALS_POOLS),$(foreach pool,$(TARGETS_FILE_INTERVALS_POOLS),$(foreach sample,$(SAMPLES),metrics/$(sample).amplicon_metrics_$(shell basename $(pool)).txt)))
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/CUSTOM_AMPLICON_SET/SAMPLE/; s/\s$$//' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.amplicon_metrics.txt}); \
		head -8 $$metrics | sed "/^#/d; /^CUSTOM_AMPLICON_SET/d; /^\$$/d; s/^tmp/$$samplename/; s/\t\+$$//" ; \
	done; \
	} > $@

metrics/all$(PROJECT_PREFIX).wgs_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).wgs_metrics.txt)
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/^/SAMPLE\t/; s/\s$$//' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.wgs_metrics.txt}); \
		head -8 $$metrics | sed "/^#/d; /^GENOME_TERRITORY/d; /^\$$/d; s/^/$$samplename\t/; s/\t\+$$//" ; \
	done; \
	} > $@


# summarize interval metrics into one file
# (basename within make will strip extension; system basename will strip dir path, and optional extension)
metrics/all$(PROJECT_PREFIX).interval_amplicon_metrics.txt.gz : $(foreach sample,$(SAMPLES),metrics/$(sample).interval_amplicon_metrics.txt.gz)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	zcat $< | sed '/^#/d; /^$$/d' | cut -f 1-6 > $(basename $@).tmp; \
	for metrics in $^; do \
		samplename=\$$(basename \$${metrics} .interval_amplicon_metrics.txt.gz); \
		zcat \$$metrics | sed '/^#/d; /^$$/d' | cut -f 7-8 | \
		sed \"s/mean_coverage/\$${samplename}_mean_coverage/; s/normalized_coverage/\$${samplename}_normalized_coverage/\" | \
		paste $(basename $@).tmp - > $(basename $@); \
		cp $(basename $@) $(basename $@).tmp; \
	done && gzip $(basename $@) && rm -f $(basename $@).tmp")

$(info SAMPLE $(SAMPLES))
metrics/all$(PROJECT_PREFIX).rnaseq_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).rnaseq_metrics.txt)
	$(INIT) \
	grep '^PF_BASES' $< > $@ && \
	samp_col=`grep '^PF_BASES' $< | tr '\t' '\n' | grep -n SAMPLE | cut -f1 -d':'` && \
	for i in $^; do sample=`basename $$i .rnaseq_metrics.txt`; \
		grep -A1 '^PF_BASES' $$i | tail -1 | awk -v sample=$$sample -v col=$$samp_col 'BEGIN { OFS = "\t" } { $$col=sample; print $$0 }'  >> $@; \
	done; 

metrics/all$(PROJECT_PREFIX).normalized_coverage.rnaseq_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).rnaseq_metrics.txt)
	$(INIT) \
	grep -A101 '^normalized_position' $< | cut -f1 > $@ && \
	for i in $^; do sample=`basename $$i .rnaseq_metrics.txt`; \
		grep -A101 '^normalized_position' $$i | cut -f2 | sed "s/All_Reads/$$sample/" | paste $@ - > $@.tmp && mv $@.tmp $@; \
	done;

metrics/all$(PROJECT_PREFIX).rnaseq_report/index.html : metrics/all.rnaseq_metrics.txt metrics/all.normalized_coverage.rnaseq_metrics.txt
	$(INIT) module load $(R_MODULE); $(PLOT_RNASEQ_METRICS) --outDir $(@D) $^

metrics/all$(PROJECT_PREFIX).alignment_summary_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).alignment_summary_metrics.txt)
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/SAMPLE.*//; s/\s$$//; s/^/SAMPLE\t/;' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.alignment_summary_metrics.txt}); \
		sed "/^#/d; /^CATEGORY/d; /^READ_LENGTH/d; /^[0-9]/d; /^\$$/d; s/^/$$samplename\t/; s/\t\+$$//" $$metrics | grep "^$$samplename"; \
	done; \
	} >$@

metrics/all$(PROJECT_PREFIX).dup_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).dup_metrics.txt)
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/\s$$//; ' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.dup_metrics.txt}); \
		head -8 $$metrics | sed "/^#/d; /^LIBRARY/d; /^\$$/d; s/\t\+$$//"; \
	done; \
	} >$@

metrics/all$(PROJECT_PREFIX).oxog_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).oxog_metrics.txt)
	$(INIT) \
	{ \
	sed '/^$$/d; /^#/d; s/\s$$//;' $< | head -1; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.oxog_metrics.txt}); \
		sed "/^#/d; /^SAMPLE_ALIAS/d; /^\$$/d; s/^[^\t]\+\t/$$samplename\t/;" $$metrics | \
		grep "^$$samplename" | sort -r -k 11 | head -1; \
	done; \
	} >$@

metrics/all$(PROJECT_PREFIX).flagstats.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).flagstats.txt)
	$(INIT) \
	{ \
	echo -ne "category\t"; sed 's/^[0-9]\+ + [0-9]\+ //;' $< | tr '\n' '\t';  echo "";\
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.flagstats.txt}); echo -ne "$$samplename\t"; \
		cut -f1 -d' ' $$metrics | tr '\n' '\t'; echo ""; \
	done; \
	} >$@

metrics/all$(PROJECT_PREFIX).flagstatsQ30.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).flagstatsQ30.txt)
	$(INIT) \
	{ \
	echo -ne "category\t"; sed 's/^[0-9]\+ + [0-9]\+ //;' $< | tr '\n' '\t';  echo "";\
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.flagstats.txt}); echo -ne "$$samplename\t"; \
		cut -f1 -d' ' $$metrics | tr '\n' '\t'; echo ""; \
	done; \
	} >$@
	
include usb-modules-v2/bam_tools/processBam.mk
