.DEFAULT_GOAL := brass

# run brass
MUT_CALLER = brass

include usb-modules-v2/Makefile.inc
include usb-modules-v2/vcf_tools/vcftools.mk

LOGDIR = log/brass.$(NOW)
PHONY += brass

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

HIGH_DEPTH_BED = brass_files/merged/depth_mask.bed
VIRAL = $(IMG_DIR)/BRASS/ref/viral/2bit/viral.genomic.fa.2bit
MICROBIAL = $(IMG_DIR)/BRASS/ref/bacterial/2bit/all_ncbi_bacteria
CENTTEL = $(IMG_DIR)/BRASS/ref/centtel.tsv
GCBINS = $(IMG_DIR)/BRASS/ref/gcBins.bed.gz
CYTOBAND = $(IMG_DIR)/BRASS/ref/cytoband.bed
GENOME_CACHE = $(GROUP_DIR)/ref/genomes/hg38/cache/vagrent.human.hg38.homo_sapiens_core_80_38.cache.gz
BRASS_PROTOCOL = WGS# WGS|WXS|RNA

BRASS_OPTS = -g $(REF_FASTA) -s HUMAN -as 38 -pr $(BRASS_PROTOCOL) -gc $(GENOME_CACHE) -d $(HIGH_DEPTH_BED) -vi $(VIRAL) -mi $(MICROBIAL) -ct $(CENTTEL) -b $(GCBINS) -cb $(CYTOBAND)

brass : brass_files/merged/depth_mask.bed $(foreach pair,$(SAMPLE_PAIRS),brass/$(pair).brass_wrap.log)

bam/%.bam.bas : bam/%.bam
	$(call RUN,1,6G,$(RESOURCE_REQ_MEDIUM),$(SINGULARITY_MODULE),"\
	$(BRASS) bam_stats \
	-i $< \
	-o $@ \
	-r $(REF_FASTA).fai")

brass_files/%.bw : bam/%.bam
	$(MKDIR) $(@D)
	$(call RUN,1,6G,$(RESOURCE_REQ_MEDIUM),$(SINGULARITY_MODULE),"\
	$(CPGBIGWIG) bam2bw \
	-i $< \
	-o $@")

brass_files/depth/%.bed : brass_files/%.bw
	$(MKDIR) $(@D)
	$(call RUN,1,6G,$(RESOURCE_REQ_MEDIUM),$(SINGULARITY_MODULE),"\
	$(CPGBIGWIG) detectExtremeDepth \
	-b $< \
	-o brass_files/depth/ ")

brass_files/sorted/%.bed: brass_files/depth/%.bed
	$(MKDIR) $(@D)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	cat $< | sort -k1,1 -k2,2n > $@")

brass_files/merged/intersect.bed: $(foreach pair,$(SAMPLE_PAIRS),brass_files/sorted/$(tumor.$(pair)).bed)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BEDTOOLS_MODULE),"\
	$(BEDTOOLS) multiinter -i $^ > $@")

brass_files/merged/intersect.filtered.bed: brass_files/merged/intersect.bed
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(PERL_MODULE),"\
	perl -ane 'next if ($$F[3] < 3); print $$F[0]$(comma)qq{\t}$(comma)$$F[1]$(comma)qq{\t}$(comma)$$F[2]$(comma)qq{\n};' $< \
	| sort -k1$(comma)1 -k2$(comma)2n \
	> $@")

brass_files/merged/depth_mask.bed.gz: brass_files/merged/intersect.filtered.bed
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BEDTOOLS_MODULE),"\
	$(BEDTOOLS) merge -i $< -d 250 \
	| perl -ane 'next if ($$F[2]-$$F[1] < 500); print;' \
	| $(BGZIP) -c > $@")

brass_files/merged/depth_mask.bed.gz.tbi: brass_files/merged/depth_mask.bed.gz
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	tabix -p bed $<")

brass_files/merged/depth_mask.bed: brass_files/merged/depth_mask.bed.gz
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	$(GUNZIP) $< > $@")


brass_files/sampstat/%.txt: facets/cncf/all.summary.txt
	$(MKDIR) $(@D)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	awk -v sample="$*" 'BEGIN{FS="\t"} NR > 1 && $$3 == sample { \
		if ($$9 != "NA" && $$10 != "NA") { \
			printf("rho\t%s\nPloidy\t%s\nGenderChr\tY\nGenderChrFound\tY\n", $$9, $$10); \
		} \
	}' $< > $@")

define brass-tumor-normal
brass/$1_$2.brass_wrap.log : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass_files/merged/depth_mask.bed brass_files/sampstat/$1.txt
	$$(call RUN,2,32G,$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(BRASS) brass.pl $$(BRASS_OPTS) \
	-ss $$(<<<<<<) \
	-t $$< \
	-n $$(<<) \
	-tn $1 \
	-nn $2 \
	-c 2 \
	-o brass/$1_$2 && touch $$@")

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call brass-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))