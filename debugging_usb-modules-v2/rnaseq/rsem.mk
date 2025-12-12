include usb-modules-v2/Makefile.inc

LOGDIR ?= log/rsem.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: rsem

rsem : $(foreach type1,genes isoforms,$(foreach type2,expected_count TPM FPKM,rsem/all$(PROJECT_PREFIX).$(type1).$(type2).results)) \
$(foreach type1,genes isoforms,$(foreach type2,expected_count TPM FPKM,rsem/all$(PROJECT_PREFIX).$(type1).$(type2).results_coding)) \
rsem/all$(PROJECT_PREFIX).genes.expected_count.results_coding_uq

define rsem-calc-expression
rsem/$1.genes.results : star/$1.Aligned.toTranscriptome.out.bam 
	$$(call RUN,$$(RSEM_NUM_CORES),$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(RSEM_MODULE),"\
	$$(RSEM_CALC_EXPR) $$(RSEM_OPTIONS) $$(if $$(findstring true,$$(PAIRED_END)),--paired-end) \
	--temporary-folder $$(TMPDIR)/rsem_$1.temp \
	$$< $$(RSEM_INDEX) $$(@D)/$1")
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call rsem-calc-expression,$(sample))))

rsem/%.isoforms.results : rsem/%.genes.results
	
define rsem-gen-dat-matrix
rsem/all$$(PROJECT_PREFIX).$1.$2.results : $$(foreach sample,$$(SAMPLES),rsem/$$(sample).$1.results)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(PERL_MODULE) $$(R_MODULE),"\
	$$(RSEM_GEN_DATA_MATRIX) $$(word 3,$$(subst .,$$(space),$$@)) $$^ | sed 's/rsem\///g;' | \
	sed \"s/\.$$(word 2,$$(subst .,$$(space),$$@))\.results//g\" | tr -d \"\\\"\" > $$@.tmp && \
	$$(RSCRIPT) $$(RSEM_PROCCESS) --inputRSEMFile $$@.tmp --gtf $$(GENCODE_GTF) --outputFile $$@ && $$(RM) $$@.tmp")
endef
$(foreach type1,genes isoforms,\
	$(foreach type2,expected_count FPKM TPM,\
		$(eval $(call rsem-gen-dat-matrix,$(type1),$(type2)))))


define rsem-gen-dat-matrix-coding
rsem/all$$(PROJECT_PREFIX).$1.$2.results_coding : $$(foreach sample,$$(SAMPLES),rsem/$$(sample).$1.results)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(PERL_MODULE) $$(R_MODULE),"\
	$$(RSEM_GEN_DATA_MATRIX) $$(word 3,$$(subst .,$$(space),$$@)) $$^ | sed 's/rsem\///g;' | \
	sed \"s/\.$$(word 2,$$(subst .,$$(space),$$@))\.results//g\" | tr -d \"\\\"\" > $$@.tmp && \
	$$(RSCRIPT) $$(RSEM_PROCCESS) --inputRSEMFile $$@.tmp --gtf $$(GENCODE_GTF) --outputFile $$@ --geneBiotype protein_coding \
	&& $$(RM) $$@.tmp")
endef
$(foreach type1,genes isoforms,\
	$(foreach type2,expected_count FPKM TPM,\
	$(eval $(call rsem-gen-dat-matrix-coding,$(type1),$(type2)))))

rsem/all$(PROJECT_PREFIX).genes.expected_count.results_coding_uq : $(foreach sample,$(SAMPLES),rsem/$(sample).genes.results)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(PERL_MODULE) $(R_MODULE),"\
	$(RSEM_GEN_DATA_MATRIX) expected_count $^ | sed 's/rsem\///g; s/\.genes\.results//g' | tr -d \"\\\"\" > $@.tmp && \
	$(RSCRIPT) $(RSEM_PROCCESS) --inputRSEMFile $@.tmp --gtf $(GENCODE_GTF) --outputFile $@ --geneBiotype protein_coding \
	--normalizationMethod \"uq\" --threshold_for_uq 1000 && $(RM) $@.tmp")	

include usb-modules-v2/bam_tools/processBam.mk
include usb-modules-v2/aligners/align.mk
include usb-modules-v2/aligners/starAligner.mk

