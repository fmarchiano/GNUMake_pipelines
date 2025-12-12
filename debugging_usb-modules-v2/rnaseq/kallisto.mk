include usb-modules-v2/Makefile.inc

LOGDIR ?= log/kallisto.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: kallisto


kallisto : $(foreach sample,$(SAMPLES),kallisto/$(sample)/abundance.tsv.gz) \
kallisto/all$(PROJECT_PREFIX).est_count.txt.gz kallisto/all$(PROJECT_PREFIX).tpm.txt.gz \
$(if $(findstring true,$(PAIRED_END)),$(foreach sample,$(SAMPLES),pizzly/$(sample).fusions.txt.gz))


define kallisto
kallisto/$1/abundance.tsv.gz: fastq/$1.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),fastq/$1.2.fastq.gz)
	$$(call RUN,4,3G,$$(RESOURCE_REQ_MEDIUM),,"\
	$$(KALLISTO) quant -i $$(KALLISTO_INDEX) $$(KALLISTO_OPTIONS) \
	$$(if $$(findstring FIRST_READ_TRANSCRIPTION_STRAND,$$(STRAND_SPECIFICITY)),--fr-stranded) \
	$$(if $$(findstring SECOND_READ_TRANSCRIPTION_STRAND,$$(STRAND_SPECIFICITY)),--rf-stranded) \
	--threads $$(KALLISTO_NUM_CORES) \
	$$(if $$(findstring true,$$(KALLISTO_BIAS)),--bias) \
	$$(if $$(findstring true,$$(KALLISTO_GENOMEBAM)),--genomebam --gtf $$(KALLISTO_GTF) --chromosomes $$(KALLISTO_CHR)) \
	-o kallisto/$1 \
	$$(if $$(findstring false,$$(PAIRED_END)),--single --fragment-length $$(KALLISTO_FRAGMENT_LEN)) \
	$$^ && \
	$$(GZIP) kallisto/$1/abundance.tsv kallisto/$1/fusion.txt")
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call kallisto,$(sample))))


kallisto/all$(PROJECT_PREFIX).est_count.txt.gz : $(foreach sample,$(SAMPLES),kallisto/$(sample)/abundance.tsv.gz)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(RSCRIPT) usb-modules-v2/rnaseq/kallisto_merge.R --prefix kallisto/all$(PROJECT_PREFIX) $^ &&\
	$(GZIP) kallisto/all$(PROJECT_PREFIX).est_count.txt kallisto/all$(PROJECT_PREFIX).tpm.txt")

kallisto/all$(PROJECT_PREFIX).tpm.txt.gz : kallisto/all$(PROJECT_PREFIX).est_count.txt.gz


ifeq ($(findstring true,$(PAIRED_END)),true)
define pizzly
pizzly/$1.fusions.txt.gz : kallisto/$1/abundance.tsv.gz
	$$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_MEDIUM),$(PYTHON_MODULE),"\
	zcat kallisto/$1/fusion.txt.gz > pizzly/$1.fusion.txt && \
	$$(PIZZLY) -k $$(PIZZLY_K) \
	--gtf $$(KALLISTO_GTF) \
	--cache $$(PIZZLY_GTF_CACHE) \
	--fasta $$(PIZZLY_FASTA) \
	--output pizzly/$1 pizzly/$1.fusion.txt && \
	$$(PYTHON) $$(PIZZLY_FLATTEN_JSON) pizzly/$1.json pizzly/$1.fusions.txt && \
	$$(PYTHON) $$(PIZZLY_FLATTEN_JSON) pizzly/$1.unfiltered.json pizzly/$1.unfiltered.fusions.txt && \
	$$(GZIP) pizzly/$1.fusions.txt \
		 pizzly/$1.unfiltered.fusions.txt \
		 pizzly/$1.fusions.fasta \
		 pizzly/$1.unfiltered.fusions.fasta \
		 pizzly/$1.json \
		 pizzly/$1.unfiltered.json && \
		 rm pizzly/$1.fusion.txt")
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call pizzly,$(sample))))
endif


include usb-modules-v2/fastq_tools/fastq.mk
