include usb-modules-v2/Makefile.inc
include usb-modules-v2/aligners/align.inc

LOGDIR ?= log/star.$(NOW)

.PHONY: star
.DELETE_ON_ERROR:

ALIGNER := star
override BAM_SORT := false
override BAM_FIX_RG := true

ifeq ($(strip $(PRIMARY_ALIGNER)),star)
STAR_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam) $(foreach sample,$(SAMPLES),star/$(sample).Aligned.sortedByCoord.out.bam)
else
STAR_BAMS = $(foreach sample,$(SAMPLES),star/$(sample).Aligned.sortedByCoord.out.bam)
endif

star : $(STAR_BAMS) $(addsuffix .bai,$(STAR_BAMS)) \
star/all$(PROJECT_PREFIX).ReadsPerGene.out.tab star/all$(PROJECT_PREFIX).ReadsPerGene.out.tab.coding \
star/all$(PROJECT_PREFIX).alignment_stats.txt \
$(foreach sample,$(SAMPLES),$(if $(findstring Junctions,$(STAR_OPTIONS)),star/$(sample).Chimeric.out.junction.gz,) $(if $(findstring SeparateSAMold,$(STAR_OPTIONS)),star/$(sample).Chimeric.out.sam.gz,) \
star/$(sample).Unmapped.out.mate1.gz $(if $(findstring true,$(PAIRED_END)),star/$(sample).Unmapped.out.mate2.gz) \
star/$(sample).Log.out.gz star/$(sample).Log.progress.out.gz)

star/%.Aligned.sortedByCoord.out.bam : fastq/%.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),fastq/%.2.fastq.gz)
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9\-]\+//'`; \
	$(call RUN,$(STAR_CPU),$(if $(findstring _,$(REF)),100G,50G),$(RESOURCE_REQ_MEDIUM),$(STAR_MODULE),"\
	$(STAR) --runMode alignReads --runThreadN $(STAR_CPU) \
	--genomeDir $(STAR_GENOME_DIR) \
	--readFilesIn $< $(if $(findstring true,$(PAIRED_END)),$(word 2,$^)) \
	--readFilesCommand gunzip -c \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMattrRGline ID:$* LB:$${LBID} PL:$(SEQ_PLATFORM) SM:$${LBID} \
	--outReadsUnmapped Fastx \
	--outFileNamePrefix $(@D)/$*. \
	--quantMode GeneCounts TranscriptomeSAM \
	$(STAR_OPTIONS) && \
	$(RMR) $(@D)/$*._STARgenome $(@D)/$*._STARpass1")

star/%.Aligned.toTranscriptome.out.bam : star/%.Aligned.sortedByCoord.out.bam
	

star/%.Chimeric.out.junction.gz : star/%.Aligned.sortedByCoord.out.bam
	@if [ -f '$(basename $@)' ]; then \
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"$(GZIP) $(basename $@)"); fi

ifeq ($(findstring SeparateSAMold,$(STAR_OPTIONS)),SeparateSAMold)
star/%.Chimeric.out.sam.gz : star/%.Aligned.sortedByCoord.out.bam
	@if [ -f '$(basename $@)' ]; then \
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"$(GZIP) $(basename $@)"); fi
endif

star/%.ReadsPerGene.out.tab : star/%.Aligned.sortedByCoord.out.bam
	

star/%.Log.final.out : star/%.Aligned.sortedByCoord.out.bam
	

star/%.Log.out.gz : star/%.Aligned.sortedByCoord.out.bam
	@if [ -f '$(basename $@)' ]; then \
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"$(GZIP) $(basename $@)"); fi

star/%.Log.progress.out.gz : star/%.Aligned.sortedByCoord.out.bam
	@if [ -f '$(basename $@)' ]; then \
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"$(GZIP) $(basename $@)"); fi

star/%.Unmapped.out.mate1.gz : star/%.Aligned.sortedByCoord.out.bam
	@if [ -f '$(basename $@)' ]; then \
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),,"$(GZIP) $(basename $@)"); fi

star/%.Unmapped.out.mate2.gz : star/%.Aligned.sortedByCoord.out.bam
	@if [ -f '$(basename $@)' ]; then \
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),,"$(GZIP) $(basename $@)"); fi

bam/%.bam : star/%.Aligned.sortedByCoord.out.bam
	$(INIT) ln -s -f ../$< $@

bam/%.bam.bai : star/%.Aligned.sortedByCoord.out.bam.bai
	$(INIT) ln -s -f ../$< $@

star/all$(PROJECT_PREFIX).ReadsPerGene.out.tab : $(foreach sample,$(SAMPLES),star/$(sample).ReadsPerGene.out.tab)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(RSCRIPT) $(STAR_PROCESS) --gtf $(GENCODE_GENE_GTF) --outputFile $@ --stranded $(STRAND_SPECIFICITY) $^")

star/all$(PROJECT_PREFIX).ReadsPerGene.out.tab.coding : $(foreach sample,$(SAMPLES),star/$(sample).ReadsPerGene.out.tab)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(RSCRIPT) $(STAR_PROCESS) --gtf $(GENCODE_GENE_GTF) --outputFile $@ \
	--stranded $(STRAND_SPECIFICITY) --geneBiotype protein_coding $^")

star/all$(PROJECT_PREFIX).alignment_stats.txt : $(foreach sample,$(SAMPLES),star/$(sample).Log.final.out)
	$(INIT) \
	{ \
	grep "|" $< | cut -f1 -d '|' | perl -p -e "s/ //g; s/\n/\t/g;"; echo ""; \
	for metrics in $^; do \
		samplename=$$(basename $${metrics%%.Log.final.out}); \
		echo -n $$samplename;\
		grep "|" $$metrics | cut -f2 -d '|' | perl -p -e "s/[ \t]//g; s/\n/\t/g;"; echo ""; \
	done; \
	} > $@



include usb-modules-v2/fastq_tools/fastq.mk
include usb-modules-v2/bam_tools/processBam.mk
include usb-modules-v2/aligners/align.mk
