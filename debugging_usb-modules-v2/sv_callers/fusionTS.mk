# Run coverageAnalysis from Torrent Suite
##### DEFAULTS ######

LOGDIR = log/fusion_ts.$(NOW)

##### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all

all : fusion_ts/Fusions.txt

define coverage-analysis
fusion_ts/coverageAnalysis/$1/$1.amplicon.cov.xls : bam/$1.bam
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(SINGULARITY_EXEC) -B $$(dir $$(REF_FASTA)) -B $$(dir $$(TARGETS_FILE_INTERVALS)) \
	$$(TS_IMG) /TS/plugin/coverageAnalysis/run_coverage_analysis.sh -ag -D ./fusion_ts/coverageAnalysis/$1 -B $$(TARGETS_FILE_INTERVALS) $$(REF_FASTA) $$<")

fusion_ts/filtering/$1.amplicon.cov.xls.filtered : fusion_ts/coverageAnalysis/$1/$1.amplicon.cov.xls
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(FUSION_TS_FILTER) \
	--min_total_reads $$(FUSION_TS_FILTER_MIN_TOTAL_READS) \
	--min_overall_e2e_percent $$(FUSION_TS_FILTER_MIN_OVERALL_E2E_PERCENT) \
	--min_breakpoint_reads $$(FUSION_TS_FILTER_MIN_BREAKPOINT_READS) \
	--min_breakpoint_e2e_percent $$(FUSION_TS_FILTER_MIN_BREAKPOINT_E2E_PERCENT) \
	--max_breakpoint_e2e_strandbias $$(FUSION_TS_FILTER_MAX_BREAKPOINT_E2E_STRANDBIAS) \
	--min_partner_genes_exprs $$(FUSION_TS_FILTER_MIN_PARTNER_GENES_EXPRS) \
	--out $$@ $$< ")

endef
$(foreach sample,$(SAMPLES),$(eval $(call coverage-analysis,$(sample))))

fusion_ts/Fusions.txt : $(foreach sample,$(SAMPLES),fusion_ts/filtering/$(sample).amplicon.cov.xls.filtered)
	$(INIT) echo -e 'sample\tcontig_id\tcontig_srt\tcontig_end\tregion_id\tattributes\tgc_count\toverlaps\tfwd_e2e\trev_e2e\ttotal_reads\tfwd_reads\trev_reads\tcov20x\tcov100x\tcov500x' > $@ ; \
	for fname in $^ ; do \
		sed '/contig_id/d' $$fname | sed "s/^/$$(basename $$fname .amplicon.cov.xls.out)\t/g" ; \
	done >> $@


