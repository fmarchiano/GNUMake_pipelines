# This module runs ContEst on snp vcf files from gatk
# Author: inodb

##### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc

LOGDIR ?= log/contest.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: contest

contest : contest/all$(PROJECT_PREFIX).contest.txt 

# ContEst on gatk snp_vcf folder
ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
define contest-tumor-normal
contest/$1_$2.contest.txt : bam/$1.bam bam/$2.bam
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
	$$(call GATK,ContEst,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $$(REF_FASTA) \
	-I:eval $$(<) -I:genotype $$(<<) --popfile $$(CONTEST_POP_FILE) \
	--min_mapq $$(MIN_MQ) --min_qscore 20 \
	$$(if $(TARGETS_FILE_INTERVALS),-L $$(TARGETS_FILE_INTERVALS)) \
	-isr INTERSECTION -o $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call contest-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

#contest/snp_vcf/%_contamination.txt : bam/%.bam gatk/dbsnp/%.gatk_snps.vcf
#	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(JAVA8_MODULE),"\
#	$(MKDIR) contest/snp_vcf; \
#	$(call GATK,ContEst,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -I $(<) -B:genotypes$(,)vcf $(<<) -o $(@)")
endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
contest/snp_vcf/%_contamination.txt : bam/%.bam tvc/dbsnp/%/TSVC_variants.vcf
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(JAVA8_MODULE),"\
	$(MKDIR) contest/snp_vcf; \
	$(call GATK,ContEst,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -I $(<) -B:genotypes$(,)vcf $(<<) -o $(@)")
endif

contest/all$(PROJECT_PREFIX).contest.txt : $(foreach pair,$(SAMPLE_PAIRS),contest/$(pair).contest.txt)
	( \
		head -1 $< | sed "s/^/sample\t/"; \
		for s in $(^); do \
			grep -P "META\t" $$s | sed "s/^/`basename $$s _contest.txt`/"; \
		done | sort -rnk5,5; \
	) > $@


include usb-modules-v2/vcf_tools/vcftools.mk
include usb-modules-v2/variant_callers/TVC.mk
include usb-modules-v2/variant_callers/gatk.mk
