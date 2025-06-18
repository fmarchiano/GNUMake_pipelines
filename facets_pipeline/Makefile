# --- Include configuration ---
include facets_pipeline/Makefile.inc


# --- Filenames ---
SNPPILEUP_SUFFIX = q$(FACETS_SNP_PILEUP_MINMAPQ)_Q$(FACETS_SNP_PILEUP_MINBASEQ)_d$(FACETS_SNP_PILEUP_MAX_DEPTH)_r$(FACETS_SNP_PILEUP_MIN_DEPTH)_P$(FACETS_SNP_PILEUP_PSEUDO_SNPS)
FACETS_SUFFIX = $(SNPPILEUP_SUFFIX)_bin$(FACETS_WINDOW_SIZE)_mingc$(FACETS_MINGC)_maxgc$(FACETS_MAXGC)_nhet$(FACETS_MIN_NHET)_cval$(FACETS_CVAL)
FACETS_TARGETS_INTERVALS = ref/annotation_data/hg38/dbsnp_151.hg38.common.vcf.gz
# --- Tools ---
SNP_PILEUP := singularity exec -B $(PWD)/bam:/bam -B $(PWD)/ref/annotation_data/hg38/dbsnp_151.hg38.common.vcf.gz snp-pileup/snp-pileup.sif /usr/local/bin/htstools/snp-pileup

# --- Targets ---
.PHONY: snp-pileup

snp-pileup: facets/snp_pileup/$(TUMOR_BAM)_$(NORMAL_BAM)_$(SNPPILEUP_SUFFIX).bc.gz

# --- Rule Template ---
define snp-pileup

facets/snp_pileup/$1_$2_$$(SNPPILEUP_SUFFIX).bc.gz : bam/$1.bam bam/$2.bam $$(FACETS_TARGETS_INTERVALS)
	@echo "Running snp-pileup for TUMOR=$1, NORMAL=$2"
	mkdir -p facets/snp_pileup
	$$(SNP_PILEUP) \
		-A -P $$(FACETS_SNP_PILEUP_PSEUDO_SNPS) -d $$(FACETS_SNP_PILEUP_MAX_DEPTH) -g \
		-q $$(FACETS_SNP_PILEUP_MINMAPQ) -Q $$(FACETS_SNP_PILEUP_MINBASEQ) \
		-r $$(FACETS_SNP_PILEUP_MIN_DEPTH),0 \
		$(PWD)/$$(word 3,$$^) $$@ /bam/$$(notdir $$(word 2,$$^)) /bam/$$(notdir $$(word 1,$$^))
endef

# --- Evaluate Rule ---
$(eval $(call snp-pileup,$(TUMOR_BAM),$(NORMAL_BAM)))


# --- Targets ---
.PHONY: facets

facets: facets/cncfTN/$(TUMOR_BAM)_$(NORMAL_BAM)_$(FACETS_SUFFIX).done facets/cncf/$(TUMOR_BAM)_$(NORMAL_BAM).Rdata \
facets/cncf/$(TUMOR_BAM)_$(NORMAL_BAM).out facets/cncf/$(TUMOR_BAM)_$(NORMAL_BAM).cncf.txt facets/cncf/$(TUMOR_BAM)_$(NORMAL_BAM).cncf.pdf \
facets/cncf/$(TUMOR_BAM)_$(NORMAL_BAM).cncf.png facets/cncf/$(TUMOR_BAM)_$(NORMAL_BAM).logR.pdf facets/cncf/$(TUMOR_BAM)_$(NORMAL_BAM).HetMarkFreq.pdf 

ifdef FACETS_SPECIAL_CASES
FACETS_SUFFIX_FINAL = $(shell grep $(TUMOR_BAM)_$(NORMAL_BAM) $(FACETS_SPECIAL_CASES) | cut -f2)
endif
FACETS_SUFFIX_FINAL ?= $(FACETS_SUFFIX)

#use singularity on host HPC but R sif and library here present in the container
define facets
facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).done : facets/snp_pileup/$1_$2_$$(SNPPILEUP_SUFFIX).bc.gz
	@echo "Running facets for TUMOR=$1, NORMAL=$2"
	mkdir -p facets/cncfTN
	$$(FACETS) --pre_cval $$(FACETS_PRE_CVAL) \
	--minNDepth $$(FACETS_SNP_PILEUP_MIN_DEPTH) \
	--maxNDepth $$(FACETS_SNP_PILEUP_MAX_DEPTH) \
	--snp_nbhd $$(FACETS_WINDOW_SIZE) \
	--minGC $$(FACETS_MINGC) --maxGC $$(FACETS_MAXGC) \
	--cval $$(FACETS_CVAL) --genome $$(REF) --min_nhet $$(FACETS_MIN_NHET) \
	--max_segs $$(FACETS_MAX_SEGS) \
	--tumorName $1 --normalName $2 \
	--outPrefix facets/cncfTN/$1_$2_$$(FACETS_SUFFIX) $$<

facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).Rdata facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).cncf.txt facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).cncf.pdf facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).cncf.png facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).logR.pdf facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).out: facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).done

# Final Rdata is in facets/cncfTN, but other files can be in facets/cncfTN/rerun.
facets/cncf/$1_$2.Rdata : facets/cncfTN/$1_$2_$$(FACETS_SUFFIX_FINAL).Rdata
	mkdir -p facets/cncf
	ln -f $$(<) $$(@)

facets/cncf/$1_$2.out : facets/cncfTN/$1_$2_$$(FACETS_SUFFIX_FINAL).out
	{ \
	if test -f facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.out; then \
		mkdir -p facets/cncf/rerun && \
		ln -f facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.out facets/cncf/rerun/$1_$2_$$$$(echo facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.out | sed 's/.*rerun_cval/rerun_cval/')  && \
		ln -f -s rerun/$$$$(basename facets/cncf/rerun/$1_$2_rerun*.out) $$(@); \
	else mkdir -p facets/cncf && ln -f $$(<) $$(@); fi; \
	}

facets/cncf/$1_$2.cncf.txt : facets/cncfTN/$1_$2_$$(FACETS_SUFFIX_FINAL).cncf.txt
	{ \
	if test -f facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.cncf.txt; then \
		mkdir -p facets/cncf/rerun && \
		ln -f facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.cncf.txt facets/cncf/rerun/$1_$2_$$$$(echo facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.cncf.txt | sed 's/.*rerun_cval/rerun_cval/')  && \
		ln -f -s rerun/$$$$(basename facets/cncf/rerun/$1_$2_rerun*.cncf.txt) $$(@); \
	else mkdir -p facets/cncf && ln -f $$(<) $$(@); fi; \
}

facets/cncf/$1_$2.cncf.pdf : facets/cncfTN/$1_$2_$$(FACETS_SUFFIX_FINAL).cncf.pdf
	{ \
	if test -f facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.cncf.pdf; then \
		mkdir -p facets/cncf/rerun && \
		ln -f facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.cncf.pdf facets/cncf/rerun/$1_$2_$$$$(echo facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.cncf.pdf | sed 's/.*rerun_cval/rerun_cval/')  && \
		ln -f -s rerun/$$$$(basename facets/cncf/rerun/$1_$2_rerun*.cncf.pdf) $$(@); \
	else mkdir -p facets/cncf && ln -f $$(<) $$(@); fi; \
}

facets/cncf/$1_$2.cncf.png : facets/cncfTN/$1_$2_$$(FACETS_SUFFIX_FINAL).cncf.png
	{ \
	if test -f facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.cncf.png; then \
		mkdir -p facets/cncf/rerun && \
		ln -f facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.cncf.png facets/cncf/rerun/$1_$2_$$$$(echo facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.cncf.png | sed 's/.*rerun_cval/rerun_cval/')  && \
		ln -f -s rerun/$$$$(basename facets/cncf/rerun/$1_$2_rerun*.cncf.png) $$(@); \
	else mkdir -p facets/cncf && ln -f $$(<) $$(@); fi; \
}

facets/cncf/$1_$2.logR.pdf : facets/cncfTN/$1_$2_$$(FACETS_SUFFIX_FINAL).logR.pdf
	{ \
	if test -f facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.logR.pdf; then \
		mkdir -p facets/cncf/rerun && \
		ln -f facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.logR.pdf facets/cncf/rerun/$1_$2_$$$$(echo facets/cncfTN/rerun/$1_$2_$$(FACETS_SUFFIX).rerun*.logR.pdf | sed 's/.*rerun_cval/rerun_cval/')  && \
		ln -f -s rerun/$$$$(basename facets/cncf/rerun/$1_$2_rerun*.logR.pdf) $$(@); \
	else mkdir -p facets/cncf && ln -f $$(<) $$(@); fi; \
}

facets/cncf/$1_$2.HetMarkFreq.pdf : facets/cncf/$1_$2.HetMarkFreq.txt

facets/cncf/$1_$2.HetMarkFreq.txt : facets/cncf/$1_$2.cncf.txt facets/cncf
	@echo "Running facets hetmarkfreq for TUMOR=$1, NORMAL=$2"
	$$(FACETS_HETMARKFREQ) --minMarks $$(FACETS_HETMARKFREQ_MINMARKS) \
	--threshold $$(FACETS_HETMARKFREQ_THRESHOLD) \
	$$<

endef
$(eval $(call facets,$(TUMOR_BAM),$(NORMAL_BAM)))

# --- Targets ---
.PHONY: summary

summary: facets/cncf/all$(PROJECT_PREFIX).summary.txt facets/cncf/all$(PROJECT_PREFIX).cncf.txt facets/cncf/all$(PROJECT_PREFIX).cncf.pdf.tar.gz facets/cncf/all$(PROJECT_PREFIX).cncf.png.tar.gz facets/cncf/all$(PROJECT_PREFIX).HetMarkFreq.txt

SAMPLE_PAIRS = $(shell ls facets/cncf/*.out | cut -d/ -f3  | cut -d. -f1 | uniq)
comma = ,

facets/cncf/all$(PROJECT_PREFIX).summary.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).out)
	{ \
	cut -f2 -d' ' $< | tr '\n' '\t'; echo ""; \
	for metrics in $^; do \
		cut -f4 -d' ' $$metrics | tr '\n' '\t'; echo ""; \
	done; \
	} > $@

facets/cncf/all$(PROJECT_PREFIX).cncf.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt)
	head -1 $< | sed 's/^/TUMOR_NORMAL\t/' > $@; \
	for cncf in $^; do \
		samplename=`basename $$cncf | sed 's/\.cncf\.txt//;'`; \
		sed '/^chrom/d; s/^23/X/;' $$cncf | sed "s/^/$$samplename\t/" >> $@; \
	done

facets/cncf/all$(PROJECT_PREFIX).cncf.pdf.tar.gz : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.pdf) $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).logR.pdf)
	tar -czf $@ $^

facets/cncf/all$(PROJECT_PREFIX).cncf.png.tar.gz : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.png) 
	tar -czf $@ $^

facets/cncf/all$(PROJECT_PREFIX).HetMarkFreq.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).HetMarkFreq.txt)
	{ \
	grep -r MB $^ | sed -E -e 's/.+\/(.+).HetMarkFreq.txt.([0-9\.]+) MB/\1\t\2/g' | sort -k2$(comma)2nr;\
	} >$@