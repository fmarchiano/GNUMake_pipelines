include usb-modules-v2/Makefile.inc
include usb-modules-v2/config.inc
#include usb-modules-v2/variant_callers/gatk.inc

LOGDIR = log/deseq.$(NOW)

DESEQ_RNW = usb-modules-v2/rnaseq/deseq.Rnw

# pheno file: sample\tpheno with header

.DELETE_ON_ERROR: 
.SECONDARY: 

.PHONY : all

#deseq_results.txt : sumreads/geneCounts.txt
deseq_results.txt : star/all.ReadsPerGene.out.tab
	mkdir -p graphics; $(SWEAVE) $(DESEQ_RNW) \
	--condition $(DESEQ_CONDITION) \
	--refCondition $(DESEQ_REF_CONDITION) \
	--altCondition $(DESEQ_ALT_CONDITION) \
	--outFile $@ $< $(DESEQ_PHENO_FILE)


