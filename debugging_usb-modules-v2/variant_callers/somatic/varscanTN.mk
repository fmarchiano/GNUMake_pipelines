# Run VarScan on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######

MUT_CALLER := varscan


##### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/varscanTN.$(NOW)

VPATH ?= bam

VARIANT_TYPES = varscan_indels 
#varscan_snps

PHONY += varscan varscan_vcfs varscan_tables
varscan : varscan_vcfs varscan_tables
varscan_vcfs : $(foreach type,$(VARIANT_TYPES),$(call MAKE_VCF_FILE_LIST,$(type)))
varscan_tables : $(foreach type,$(VARIANT_TYPES),$(call MAKE_TABLE_FILE_LIST,$(type)))

define varscan-somatic-tumor-normal-chr
varscan/chr_vcf/$1_$2.$3.varscan_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SAMTOOLS_MODULE) $$(JAVA8_MODULE),"\
	$$(SAMTOOLS) mpileup -r $3 -q $$(MIN_MQ) -d 10000 -f $$(REF_FASTA) $$(word 2,$$^) $$< | \
	$$(VARSCAN) somatic - varscan/chr_vcf/$1_$2.$3 --output-vcf 1 --mpileup 1 --min-var-freq $(MIN_AF_SNP) && \
	$$(VARSCAN) processSomatic varscan/chr_vcf/$1_$2.$3.snp.vcf --min-tumor-freq $$(MIN_AF_SNP) && \
	$$(VARSCAN) processSomatic varscan/chr_vcf/$1_$2.$3.indel.vcf --min-tumor-freq $$(MIN_AF_INDEL) && \
	$$(RM) varscan/chr_vcf/$1_$2.$3.snp.vcf varscan/chr_vcf/$1_$2.$3.indel.vcf && \
	$$(RM) varscan/chr_vcf/$1_$2.$3.indel.Germline* varscan/chr_vcf/$1_$2.$3.indel.LOH* && \
	$$(RM) varscan/chr_vcf/$1_$2.$3.snp.Germline* varscan/chr_vcf/$1_$2.$3.snp.LOH* && touch $$@")
varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf : varscan/chr_vcf/$1_$2.$3.varscan_timestamp
varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf : varscan/chr_vcf/$1_$2.$3.varscan_timestamp

varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.fpft.vcf : varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf bam/$1.bam bam/$1.bam.bai
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(JAVA8_MODULE),"\
	awk '! /\#/' $$< | \
	awk '{if(length(\$$$$4) > length(\$$$$5)) print \$$$$1\"\t\"(\$$$$2-1)\"\t\"(\$$$$2+length(\$$$$4)-1); \
	else print \$$$$1\"\t\"(\$$$$2-1)\"\t\"(\$$$$2+length(\$$$$5)-1)}' > varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf.region && \
	$$(BAM_READCOUNT) -f $$(REF_FASTA) -l varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf.region $$(word 2,$$^) > \
	varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf.bamrc && \
	$$(VARSCAN) fpfilter $$< varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf.bamrc \
	--output-file $$@ --filtered-file varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.fpfail.vcf \
	--min-var-freq $$(MIN_AF_SNP) --min-ref-readpos 0 --min-var-readpos 0 --min-ref-dist3 0 --min-var-dist3 0 && \
	rm varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf.region varscan/chr_vcf/$1_$2.$3.snp.Somatic.hc.vcf.bamrc")

varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.fpft.vcf : varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf bam/$1.bam bam/$1.bam.bai
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(JAVA8_MODULE),"\
	awk '! /\#/' $$< | \
	awk '{if(length(\$$$$4) > length(\$$$$5)) print \$$$$1\"\t\"(\$$$$2-1)\"\t\"(\$$$$2+length($$$$4)-1); \
	else print \$$$$1\"\t\"(\$$$$2-1)\"\t\"(\$$$$2+length(\$$$$5)-1)}' > varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf.region && \
	$$(BAM_READCOUNT) -f $$(REF_FASTA) -l varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf.region $$(word 2,$$^) > \
	varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf.bamrc && \
	$$(VARSCAN) fpfilter $$< varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf.bamrc \
	--output-file $$@ --filtered-file varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.fpfail.vcf \
	--min-var-freq $$(MIN_AF_INDEL) --min-ref-readpos 0 --min-var-readpos 0 --min-ref-dist3 0 --min-var-dist3 0 && \
	rm varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf.region varscan/chr_vcf/$1_$2.$3.indel.Somatic.hc.vcf.bamrc")

endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call varscan-somatic-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

define varscan-somatic-tumor-normal-merge
varscan/vcf/$1_$2.%.vcf : $$(foreach chr,$$(CHROMOSOMES),varscan/chr_vcf/$1_$2.$$(chr).%.vcf)
	$$(INIT) module load $$(PERL_MODULE); grep '^#' $$< > $$@; cat $$^ | grep -v '^#' >> $$@ 2> $$(LOG)

endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call varscan-somatic-tumor-normal-merge,$(tumor.$(pair)),$(normal.$(pair)))))

define varscan-somatic-fix
vcf/$1_$2.varscan_indels.vcf : varscan/vcf/$1_$2.indel.Somatic.hc.fpft.vcf
	$$(INIT) $$(FIX_VARSCAN_VCF) -t $1 -n $2 $$< | $$(VCF_SORT) $$(REF_DICT) - > $$@

vcf/$1_$2.varscan_snps.vcf : varscan/vcf/$1_$2.snp.Somatic.hc.fpft.vcf
	$$(INIT) $$(FIX_VARSCAN_VCF) -t $1 -n $2 $$< | $$(VCF_SORT) $$(REF_DICT) - > $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call varscan-somatic-fix,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/vcf_tools/vcftools.mk

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: $(PHONY)

