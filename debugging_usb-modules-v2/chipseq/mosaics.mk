include usb-modules-v2/Makefile.inc

LOGDIR ?= log/mosaics.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : mosaics

ifeq ($(findstring true,$(PAIRED_END)),true)
MOSAICS_SUFFIX ?= bin$(MOSAICS_BIN_SIZE)
else
MOSAICS_SUFFIX ?= fragL$(MOSAICS_FRAG_LEN)_bin$(MOSAICS_BIN_SIZE)
endif

mosaics : $(foreach pair,$(SAMPLE_PAIRS),mosaics/rdata_$(MOSAICS_SUFFIX)/$(pair).rdata mosaics/peaks_$(MOSAICS_SUFFIX)/$(pair).peakTFBS.annotated.bed)

define mosaics-peaks
mosaics/rdata_$$(MOSAICS_SUFFIX)/$1_$2.rdata : mosaics/bin/$1.bam_$$(MOSAICS_SUFFIX).txt mosaics/bin/$2.bam_$$(MOSAICS_SUFFIX).txt
	$$(call RUN,$$(MOSAICS_NUM_CORES),$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(R_MODULE) $$(BEDTOOLS_MODULE),"\
	$$(MKDIR) mosaics mosaics/peaks_$$(MOSAICS_SUFFIX) mosaics/rdata_$$(MOSAICS_SUFFIX) mosaics/plots_$$(MOSAICS_SUFFIX); \
	$$(MOSAICS_RUN) --parallel $$(MOSAICS_PARALLEL) --num_cores $$(MOSAICS_NUM_CORES) \
	--plotFileLoc mosaics/plots_$$(MOSAICS_SUFFIX) \
	--peakFileLoc mosaics/peaks_$$(MOSAICS_SUFFIX) \
	--rdataFileLoc mosaics/rdata_$$(MOSAICS_SUFFIX) \
	--maxgap $$(MOSAICS_MAXGAP) --minsize $$(MOSAICS_MINSIZE) --thres $$(MOSAICS_THRES) $$^ && \
	$$(BEDTOOLS) sort -i mosaics/peaks_$$(MOSAICS_SUFFIX)/$$(subst .rdata,,$$(@F)).peakTFBS.bed -faidx $$(REF_FASTA).fai > tmp && \
	mv tmp mosaics/peaks_$$(MOSAICS_SUFFIX)/$$(subst .rdata,,$$(@F)).peakTFBS.bed && \
	head -1 mosaics/peaks_$$(MOSAICS_SUFFIX)/$$(subst .rdata,,$$(@F)).peakTFBS.txt > tmp && \
	tail -n +2 mosaics/peaks_$$(MOSAICS_SUFFIX)/$$(subst .rdata,,$$(@F)).peakTFBS.txt | \
	$$(BEDTOOLS) sort -faidx $$(REF_FASTA).fai >> tmp && mv tmp mosaics/peaks_$$(MOSAICS_SUFFIX)/$$(subst .rdata,,$$(@F)).peakTFBS.txt") 

mosaics/peaks_$$(MOSAICS_SUFFIX)/$1_$2.peakTFBS.bed : mosaics/rdata_$$(MOSAICS_SUFFIX)/$1_$2.rdata
	

mosaics/peaks_$$(MOSAICS_SUFFIX)/$1_$2.peakTFBS.txt : mosaics/rdata_$$(MOSAICS_SUFFIX)/$1_$2.rdata
	
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call mosaics-peaks,$(tumor.$(pair)),$(normal.$(pair)))))

mosaics/%.coding_genes.bed : mosaics/%.bed
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BEDTOOLS_MODULE),"\
	$(BEDTOOLS) closest -g $(REF_FASTA).fai -k 1 -t all -D b -a $< -b $(GENCODE_CODING_GENE_GTF) > $@")

mosaics/%.noncoding_genes.bed : mosaics/%.bed
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BEDTOOLS_MODULE),"\
	$(BEDTOOLS) closest -g $(REF_FASTA).fai -k 1 -t all -D b -a $< -b $(GENCODE_NONCODING_GENE_GTF) > $@")

mosaics/%.closest_genes.bed : mosaics/%.bed
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BEDTOOLS_MODULE),"\
	$(BEDTOOLS) closest -g $(REF_FASTA).fai -k 1 -t all -D b -a $< -b $(GENCODE_GENE_GTF) > $@")

mosaics/bin/%.bam_$(MOSAICS_SUFFIX).txt : bam/%.bam
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(R_MODULE) $(PERL_MODULE),"\
	$(MKDIR) mosaics mosaics/bin mosaics/wig; \
	$(MOSAICS_CONSTRUCTBINS) $(if $(findstring false,$(PAIRED_END)),--fragLen $(MOSAICS_FRAG_LEN),) \
	--binSize $(MOSAICS_BIN_SIZE) --pet $(PAIRED_END) \
	--outfileLoc mosaics/bin $(if $(MOSAICS_CHRFILE),--chrfile $(MOSAICS_CHRFILE)) $<")
