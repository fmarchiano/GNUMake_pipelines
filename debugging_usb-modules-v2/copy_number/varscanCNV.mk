	# Run VarScan to detect copynumber
# Detect copy number
##### DEFAULTS ######

LOGDIR ?= log/varscanCNV.$(NOW)

##### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all

all : varscan/segment/all$(PROJECT_PREFIX).geneCN.GL_LRR.pdf varscan/segment/all$(PROJECT_PREFIX).geneCN.log2Ratio.pdf

ifeq ($(CAPTURE_METHOD),PCR)
define varscan-copynum-tumor-normal
varscan/copynum/$1_$2.$$(notdir $3).copynumber : bam/$1.bam bam/$2.bam $3
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(SAMTOOLS_MODULE) $$(JAVA7_MODULE) $$(BEDTOOLS_MODULE),"\
	TMP1=`mktemp`.bam && $$(BEDTOOLS) intersect -wa -F 1 -a $$< -b $$(word 3, $$^) > \$$$$TMP1 && \
	$$(SAMTOOLS) index \$$$$TMP1 && \
	TMP2=`mktemp`.bam && $$(BEDTOOLS) intersect -wa -F 1 -a $$(word 2,$$^) -b $$(word 3, $$^) > \$$$$TMP2 && \
	$$(SAMTOOLS) index \$$$$TMP2 && \
	$$(SAMTOOLS) mpileup $$(CBS_MPILEUP_OPTS) -l $$(word 3,$$^) -f $$(REF_FASTA) \$$$$TMP2 \$$$$TMP1 | \
	awk 'NF == 9 { print }' |  \
	$$(VARSCAN) copynumber - $$(basename $$@) --mpileup 1 --max-segment-size $$(VARSCAN_CNV_MAX_SEG_SIZE) && \
	$$(RM) \$$$$TMP1 \$$$$TMP2")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach pool,$(TARGETS_FILE_INTERVALS_POOLS_CNA),\
		$(eval $(call varscan-copynum-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)),$(pool)))))
else
define varscan-copynum-tumor-normal
varscan/copynum/$1_$2.copynumber : bam/$1.bam bam/$2.bam
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_SHORT),$$(SAMTOOLS_MODULE) $$(JAVA7_MODULE),"\
	$$(SAMTOOLS) mpileup $$(CBS_MPILEUP_OPTS) -l $$(TARGETS_FILE_INTERVALS) -f $$(REF_FASTA) $$(word 2,$$^) $$< \
	| awk 'NF == 9 { print }' |  \
	$$(VARSCAN) copynumber - $$(basename $$@) --mpileup 1")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call varscan-copynum-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif


varscan/copycall/%.copycall : varscan/copynum/%.copynumber
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA7_MODULE),"\
	n=`awk '{ total += \$$7 } END { printf \"\%.6f\"$(,) total / NR }' $<`; \
	if [ \$$(bc <<< \"\$$n > 0\") -eq 1 ]; then \
		recenter_opt=\"--recenter-up \$$n\"; \
	else \
		n=\$$(bc <<< \"\$$n*-1\"); \
		recenter_opt=\"--recenter-down \$$n\"; \
	fi; \
	$(VARSCAN) copyCaller $< --output-file $@ \$$recenter_opt")

varscan/rna_logratio/control.rna_logratio.txt : $(foreach normal,$(PANEL_OF_NORMAL_SAMPLES),star/$(normal).ReadsPerGene.out.tab)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(VARSCAN_CONTROL_PROFILE_FROM_STAR) --outFile $@ --gtf $(GENCODE_GENE_GTF) $^")

define rna-logratio
varscan/rna_logratio/$1.rna_logratio.txt : star/$1.ReadsPerGene.out.tab varscan/rna_logratio/control.rna_logratio.txt
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(VARSCAN_RNA_LOGRATIO) --tumor_file $$< --normal_file $$(<<) --gtf $$(GENCODE_GENE_GTF) --outfile $$@")
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call rna-logratio,$(sample))))


ifeq ($(CAPTURE_METHOD),PCR)
define varscan-segment
varscan/segment/$1_$2.segment.Rdata : $$(foreach pool,$$(TARGETS_FILE_INTERVALS_POOLS_CNA),varscan/copycall/$1_$2.$$(notdir $$(pool)).copycall)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(CBS_SEGMENTCNV) --alpha $$(CBS_SEG_ALPHA) --smoothRegion $$(CBS_SEG_SMOOTH) \
	--trim $$(CBS_TRIM) --clen $$(CBS_CLEN) --undoSD $$(CBS_SEG_SD) --separate_arm_seg TRUE \
	--excl_N_outlier_pc $$(CBS_EXCL_N_OUTLIER_PC) --minNdepth $$(CBS_MIN_N_DEPTH) \
	$$(if $$(CENTROMERE_TABLE),--centromereFile=$$(CENTROMERE_TABLE)) --prefix=$$(@D)/$1_$2 $$^")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call varscan-segment,$(tumor.$(pair)),$(normal.$(pair)))))
else
ifeq ($(CAPTURE_METHOD),RNA)
varscan/segment/%.segment.Rdata : varscan/rna_logratio/%.rna_logratio.txt
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(CBS_SEGMENTCNV) --alpha $(CBS_SEG_ALPHA) --smoothRegion $(CBS_SEG_SMOOTH) \
	--trim $(CBS_TRIM) --clen $(CBS_CLEN) --undoSD $(CBS_SEG_SD) \
	--excl_N_outlier_pc $(CBS_EXCL_N_OUTLIER_PC) --outlierSDscale $(CBS_OUTLIER_SD_SCALE) \
	--minNdepth $(CBS_MIN_N_DEPTH) --maxNdepth $(CBS_MAX_N_DEPTH) --minTdepth $(CBS_MIN_T_DEPTH) \
	$(if $(CENTROMERE_TABLE),--centromereFile=$(CENTROMERE_TABLE)) --prefix=$(@D)/$* $^")
else
varscan/segment/%.segment.Rdata : varscan/copycall/%.copycall
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(CBS_SEGMENTCNV) --alpha $(CBS_SEG_ALPHA) --smoothRegion $(CBS_SEG_SMOOTH) \
	--trim $(CBS_TRIM) --clen $(CBS_CLEN) --undoSD $(CBS_SEG_SD) \
	$(if $(CENTROMERE_TABLE),--centromereFile=$(CENTROMERE_TABLE)) --prefix=$(@D)/$* $^")
endif
endif

varscan/segment/%.collapsed_seg.txt : varscan/segment/%.segment.Rdata
	
ifeq ($(CAPTURE_METHOD),RNA)
varscan/segment/all$(PROJECT_PREFIX).geneCN.GL_LRR.txt : $(foreach sample,$(SAMPLES),varscan/segment/$(sample).collapsed_seg.txt)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(VARSCAN_GENE_CN) $(VARSCAN_GENE_CN_OPTS) --genesFile $(TARGETS_FILE_GENES) --outFile $(@D)/all$(PROJECT_PREFIX).geneCN $^")
else
varscan/segment/all$(PROJECT_PREFIX).geneCN.GL_LRR.txt : $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).collapsed_seg.txt)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(VARSCAN_GENE_CN) $(VARSCAN_GENE_CN_OPTS) --genesFile $(TARGETS_FILE_GENES) --outFile $(@D)/all$(PROJECT_PREFIX).geneCN $^")	
endif

varscan/segment/all$(PROJECT_PREFIX).geneCN.log2Ratio.txt : varscan/segment/all$(PROJECT_PREFIX).geneCN.GL_LRR.txt
	

varscan/segment/all$(PROJECT_PREFIX).geneCN.%.pdf : varscan/segment/all$(PROJECT_PREFIX).geneCN.%.txt
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(FACETS_GENE_CN_PLOT) $(FACETS_GENE_CN_PLOT_OPTS) $< $@")

define varscan-segment-sd-alpha-smooth
varscan/segment_sd$1_alpha$2_smooth$3/%.segment.Rdata : varscan/copycall/%.copycall
	$$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(CBS_SEGMENTCNV) --undoSD $1 --alpha $2 --smoothRegion $3 \
	$$(if $$(CENTROMERE_TABLE),--centromereFile=$$(CENTROMERE_TABLE)) --prefix=$$(@D)/$$* $$<")
endef
$(foreach sd,$(SEG_SDS),\
	$(foreach alpha,$(SEG_ALPHAS),\
		$(foreach smooth,$(SEG_SMOOTHS),\
			$(eval $(call varscan-segment-sd-alpha-smooth,$(sd),$(alpha),$(smooth))))))


