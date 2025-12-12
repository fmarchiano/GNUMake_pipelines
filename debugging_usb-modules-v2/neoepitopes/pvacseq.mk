include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/pvacseq.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: pvacseq

pvacseq : $(foreach pair,$(SAMPLE_PAIRS),$(foreach prefix,$(CALLER_PREFIX),pvacseq/results/MHC_Class_I/$(pair).$(prefix).final.tsv))

CHECK_HLA_GENOTYPE = if [ $1 = "HLA-,HLA-,HLA-,HLA-,HLA-,HLA-" ]; then $2; else $3; fi
CHECK_VCF_CMD = if [ `grep -v '^\#' $1 | wc -l` -eq 0 ] && [ `grep '^\#CHROM' $1 | wc -l` -eq 1 ]

PVACSEQ_HEADER = Chromosome Start Stop Reference Variant Transcript Ensembl_Gene_ID Variant_Type Mutation Protein_Position Gene_Name HLA_Allele Peptide_Length Sub-peptide_Position Mutation_Position MT_Epitope_Seq WT_Epitope_Seq Best_MT_Score_Method Best_MT_Score Corresponding_WT_Score Corresponding_Fold_Change Tumor_DNA_Depth Tumor_DNA_VAF Tumor_RNA_Deptg Tumor_RNA_VAF Normal_Depth Normal_VAF Gene_Expression Transcript_Expression Median_MT_Score Median_WT_Score Median_Fold_Change NetMHCpan WT_Score NetMHCpan MT_Score

pvacseq/input/%.vep.vcf : pvacseq/input/%.vcf
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(VEP_MODULE),"\
	$(VEP) --plugin Downstream --plugin Wildtype --dir_plugins $(VEP_PLUGIN_DIR) --assembly GRCh37 \
	--input_file $< --output_file $@")

define prepare-pvacseq
pvacseq/input/$1_$2.$3.vcf : vcf/$1_$2.$$(call VCF_SUFFIXES,$3).vcf
	$$(INIT) $$(PVACSEQ_PREPARE_VCF) < $$< > $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach prefix,$(CALLER_PREFIX),\
		$(eval $(call prepare-pvacseq,$(tumor.$(pair)),$(normal.$(pair)),$(prefix)))))

define run-pvacseq
pvacseq/results/MHC_Class_I/$1_$2.$3.final.tsv : pvacseq/input/$1_$2.$3.vep.vcf pvacseq/optitype/$2.results.tsv
	$$(RM) pvacseq/tmp/$1_$2_3; \
	HLA_GENOTYPE=`tail -1 $$(<<) | cut -f2-7 | sed 's/^/HLA-/; s/\t/,HLA-/g'`; \
	$$(call CHECK_HLA_GENOTYPE,$$$$HLA_GENOTYPE,echo $$(PVACSEQ_HEADER) > $$@,\
		$$(call CHECK_VCF_HIGH_MODERATE_CMD,$$<,echo $$(PVACSEQ_HEADER) > $$@,\
		$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_SHORT),$$(SINGULARITY_MODULE),"\
		$$(SINGULARITY_EXEC) $$(PVACTOOLS_IMG) $$(PVACSEQ) run ./$$< $1_$2 $$$$HLA_GENOTYPE $$(PVACSEQ_ALGO) \
		./pvacseq/tmp/$1_$2.$3 $$(PVACSEQ_OPTS) && \
		mv pvacseq/tmp/$1_$2.$3/MHC_Class_I/$1_$2.final.tsv $$@")))
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach prefix,$(CALLER_PREFIX),\
		$(eval $(call run-pvacseq,$(tumor.$(pair)),$(normal.$(pair)),$(prefix)))))

#define run-optitype
#pvacseq/optitype/$1.results.tsv : fastq/$1.1.fastq $$(if $$(PAIRED_END),fastq/$1.2.fastq,)
#	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
#		$$(SINGULARITY_RUN) $$(OPTITYPE_IMG) -i $^ $$(if NOT RNA,--dna,--rna) -o pvacseq/optitype/$1 && \
#		mv $(shell ls) && $(RMR) xx")
#endef
#$(foreach normal,$(NORMAL_SAMPLES),\
#	$(eval $(call run-optitype,$(normal.$(pair)))))



