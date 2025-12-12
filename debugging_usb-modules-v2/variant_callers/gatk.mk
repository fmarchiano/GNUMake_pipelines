
# common GATK steps for variant calling
# Author: Raymond Lim <raylim@mm.st> & Fong Chun Chan <fongchunchan@gmail.com>
#
# compared to the previous iteration
# 1, moved to GVCF mode in GATK4 HaplotypeCaller (with padding) -> 
#        GenomicsDBImport (no padding) -> GenotypeGVCFs ->
#        split into SNPs/INDELs -> filter (see below)
# 2, GVCF means sample set pairs no longer required or supported
# 3, Current version of GenomicsDB only takes one contig so SPLIT_CHR is no longer an option
# 4, finally really implemented the variantRecalibration properly
# Note: no longer splitting by chromosomes, rather split the genomic into 
#        $GATK_INTERVALS_COUNT number of intervals
# the options are
# 1, GATK_HARD_FILTER_SNPS/GATK_HARD_FILTER_INDELS: if false, use VariantRecalibration



ifndef GATK_MK

include usb-modules-v2/Makefile.inc

HAPLOTYPE_CALLER_OPTS = --dbsnp $(DBSNP) -ERC GVCF -R $(REF_FASTA)

###### RECIPES #######

%.intervals : %.vcf
	$(INIT) sed '/^#/d' $< | awk '{print $$1":"$$2 }' > $@

#gatk/dbsnp/%.gatk_snps.vcf.gz : bam/%.bam bam/%.bai
#	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(GATK42_MODULE),"\
#	$(call GATK42,HaplotypeCaller,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
#	-R $(REF_FASTA) -I $< -O $@ \
#	--dbsnp $(DBSNP_TARGETS_INTERVALS) --alleles $(DBSNP_TARGETS_INTERVALS) \
#	--genotyping-mode GENOTYPE_GIVEN_ALLELES \
#	--output-mode EMIT_ALL_SITES -ERC GVCF")

gatk/dbsnp/%.gatk_snps.vcf : bam/%.bam bam/%.bai
	$(call RUN,4,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call GATK,UnifiedGenotyper,$(RESOURCE_REQ_LOW_MEM)) \
	-nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP_TARGETS_INTERVALS) $(foreach bam,$(filter %.bam,$<),-I $(bam) ) \
	--genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles $(DBSNP_TARGETS_INTERVALS_COMMON) -o $@ --output_mode EMIT_ALL_SITES")

gatk/vcf_ug/%.variants.vcf : bam/%.bam bam/%.bai
	$(call RUN,4,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call GATK,UnifiedGenotyper,$(RESOURCE_REQ_LOW_MEM)) -nt 4 -R $(REF_FASTA) \
	--dbsnp $(DBSNP_TARGETS_INTERVALS) $(foreach bam,$(filter %.bam,$<),-I $(bam) ) -o $@")

#### Section for Haplotype caller for variant calling

MAX_INTERVAL_IDX = $(shell expr $(GATK_INTERVALS_COUNT) - 1)

gatk/intervals/$(shell printf "%04d" $(MAX_INTERVAL_IDX))-scattered.intervals : $(TARGETS_FILE_INTERVALS)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	$(call GATK42,SplitIntervals,$(RESOURCE_REQ_LOW_MEM_JAVA)) \
	-L $< -R $(REF_FASTA) -O $(dir $@) --scatter-count $(GATK_INTERVALS_COUNT)") 

define split-intervals
gatk/intervals/$1-scattered.intervals : gatk/intervals/$$(shell printf "%04d" $$(MAX_INTERVAL_IDX))-scattered.intervals
	

endef
$(foreach interval,$(shell seq 0 $(shell expr $(MAX_INTERVAL_IDX) - 1)),$(eval $(call split-intervals,$(shell printf "%04d" $(interval)))))

define hapcall-interval
gatk/intervals_gvcf/$1/%.variants.vcf.gz : bam/%.bam bam/%.bai gatk/intervals/$1-scattered.intervals 
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),,"\
	$$(MKDIR) gatk/intervals_gvcf/; \
	$$(call GATK42,HaplotypeCaller,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) $$(HAPLOTYPE_CALLER_OPTS) \
	-I $$< -O $$@ -L $$(<<<) -ip 150")
endef
$(foreach interval,$(shell seq 0 $(MAX_INTERVAL_IDX)),$(eval $(call hapcall-interval,$(shell printf "%04d" $(interval)))))

define genomicsdb-interval
gatk/intervals_db/$1.db.vcf.gz : $$(foreach sample,$$(SAMPLES),gatk/intervals_gvcf/$1/$$(sample).variants.vcf.gz) gatk/intervals/$1-scattered.intervals
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_SHORT),,"\
	$$(MKDIR) gatk/intervals_db/; $$(RMR) $$(@D)/$1; \
	$$(call GATK42,GenomicsDBImport,$$(RESOURCE_REQ_HIGH_MEM_JAVA)) \
	$$(foreach vcf,$$(filter %.vcf.gz,$$^),-V $$(vcf) ) --genomicsdb-workspace-path $$(@D)/$1 \
	-L $$(lastword $$^) && \
	$$(call GATK42,GenotypeGVCFs,$$(RESOURCE_REQ_HIGH_MEM_JAVA)) -R $$(REF_FASTA) \
	-V gendb://$$(@D)/$1 -O $$@ && $$(RMR) $$(@D)/$1")
endef
$(foreach interval,$(shell seq 0 $(MAX_INTERVAL_IDX)),$(eval $(call genomicsdb-interval,$(shell printf "%04d" $(interval)))))

gatk/gvcf/genomics.variants.vcf.gz : $(foreach interval,$(shell seq 0 $(MAX_INTERVAL_IDX)),gatk/intervals_db/$(shell printf "%04d" $(interval)).db.vcf.gz)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),,"\
	$(call GATK42,SortVcf,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	$(foreach vcf,$^,-I $(vcf) ) -O $@")
	


####### END HaplotypeCaller


gatk/gatk_snps.vcf.gz : gatk/gvcf/genomics.variants.vcf.gz gatk/gvcf/genomics.variants.vcf.gz.tbi
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),,"\
	$(call GATK42,SelectVariants,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $(REF_FASTA) -V $< -O $@ -select-type SNP")

gatk/gatk_indels.vcf.gz : gatk/gvcf/genomics.variants.vcf.gz gatk/gvcf/genomics.variants.vcf.gz.tbi
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),,"\
	$(call GATK42,SelectVariants,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $(REF_FASTA) -V $< -O $@ -select-type INDEL")


###### Here are the options for filtering HaplotypeCaller Results
##### Different for snps and indels
# a, use variant Recalibration (recommended)
# b, use hard filters (for very small sample size)

HAPCALL_SNP_FILTERS := --filterName 'QD' --filterExpression 'QD < $(HAPCALL_SNP_QD_THRESHOLD)' \
--filterName 'MQ' --filterExpression 'MQ < $(HAPCALL_SNP_MQ_THRESHOLD)' \
--filterName 'FS' --filterExpression 'FS > $(HAPCALL_SNP_FS_THRESHOLD)' \
--filterName 'HapScore' --filterExpression 'HaplotypeScore > $(HAPCALL_SNP_HAP_SCORE_THRESHOLD)' \
--filterName 'MQRankSum' --filterExpression 'MQRankSum < $(HAPCALL_SNP_MQ_RANKSUM_THRESHOLD)' \
--filterName 'ReadPosRankSum' --filterExpression 'ReadPosRankSum < $(HAPCALL_SNP_READPOS_RANKSUM_THRESHOLD)' \
--filterName 'Depth' --filterExpression 'DP < 5'

HAPCALL_INDEL_FILTERS = --filterName 'QD' --filterExpression 'QD < $(HAPCALL_INDEL_QD_THRESHOLD)' \
--filterName 'ReadPosRankSum' --filterExpression 'ReadPosRankSum < ($HAPCALL_INDEL_MQ_RANKSUM_THRESHOLD)' \
--filterName 'InbreedingCoeff' --filterExpression 'InbreedingCoeff < $(HAPCALL_INDEL_INBREED_COEFF_THRESHOLD)'  \
--filterName 'FS' --filterExpression 'FS > $(HAPCALL_INDEL_FS_THRESHOLD)' \
--filterName 'DP' --filterExpression 'DP < 5'

# not used InbreedingCoeff
VARIANT_RECAL_ANNOTATIONS_SNPS = QD MQ MQRankSum ReadPosRankSum FS SOR DP 
VARIANT_RECAL_ANNOTATIONS_INDELS = QD DP FS SOR ReadPosRankSum MQRankSum 

define SELECT_SAMPLE
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),,"\
	$(call GATK42,SelectVariants,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $(REF_FASTA) -V $1 -O $2 -sn $3")
endef
	
define VARIANT_RECAL_SNPS
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(R_MODULE),"\
	$(call GATK42,VariantRecalibrator,$(RESOURCE_REQ_LOW_MEM_JAVA)) \
	-R $(REF_FASTA) -V $1 -O $2 -mode SNP \
	--resource hapmap$(,)known=false$(,)training=true$(,)truth=true$(,)prior=15.0:$(HAPMAP) \
	--resource omni$(,)known=false$(,)training=true$(,)truth=true$(,)prior=12.0:$(OMNI) \
	--resource 1000G$(,)known=false$(,)training=true$(,)truth=false$(,)prior=10.0:$(R1000G) \
	--resource dbsnp$(,)known=true$(,)training=false$(,)truth=false$(,)prior=2.0:$(DBSNP) \
	$(foreach i,$(VARIANT_RECAL_ANNOTATIONS_SNPS), -an $i) \
	-tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	--tranches-file $(basename $2).tranches")
endef

define VARIANT_RECAL_INDELS
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(call GATK42,VariantRecalibrator,$(RESOURCE_REQ_LOW_MEM_JAVA)) \
	-R $(REF_FASTA) -V $1 -O $2 -mode INDEL --max-gaussians 4 \
	--resource mills$(,)known=false$(,)training=true$(,)truth=true$(,)prior=12.0:$(KNOWN_INDELS) \
	--resource dbsnp$(,)known=true$(,)training=false$(,)truth=false$(,)prior=2.0:$(DBSNP) \
	$(foreach i,$(VARIANT_RECAL_ANNOTATIONS_INDELS), -an $i) \
	-tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
	--tranches-file $(basename $2).tranches")
endef

define APPLY_VARIANT_RECAL_SNPS
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(R_MODULE),"\
	$(call GATK42,ApplyVQSR,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $(REF_FASTA) -V $1 --recal-file $2 -O $3 -mode SNP \
	--truth-sensitivity-filter-level $(VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL_SNPS) \
	--tranches-file $(basename $2).tranches")
endef

define APPLY_VARIANT_RECAL_INDELS
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(call GATK42,ApplyVQSR,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $(REF_FASTA) -V $1 --recal-file $2 -O $3 -mode INDEL\
	--truth-sensitivity-filter-level $(VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL_INDELS) \
	--tranches-file $(basename $2).tranches")
endef

## recipes for SNPs
ifeq ($(GATK_HARD_FILTER_SNPS),true)
define SELECT_SAMPLE_SNPS
gatk/vcf/$1.gatk_snps.vcf.gz : gatk/gatk_snps.vcf.gz
	$(call SELECT_SAMPLE,$^,$@,$1)
endef
$(foreach sample,$(SAMPLES),$(eval $(call SELECT_SAMPLE_SNPS,$(sample))))

gatk/vcf/%.gatk_snps.filtered.vcf.gz : gatk/vcf/%.gatk_snps.vcf.gz
	$(call RUN,4,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),,"\
	$(call GATK42,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $(REF_FASTA) \
	$(HAPCALL_SNP_FILTERS) -V $< -O $@")

else # else = no hard filter so use variantrecal

gatk/gatk_snps.recal.vcf.gz : gatk/gatk_snps.vcf.gz
	$(call VARIANT_RECAL_SNPS,$^,$@)

gatk/gatk_snps.filtered.vcf.gz : gatk/gatk_snps.vcf.gz gatk/gatk_snps.recal.vcf.gz
	$(call APPLY_VARIANT_RECAL_SNPS,$<,$(word 2,$^),$@)
endif # end GATK_HARD_FILTER_SNPS


ifeq ($(GATK_HARD_FILTER_INDELS),true)
define SELECT_SAMPLE_INDELS
gatk/vcf/$1.gatk_indels.vcf.gz : gatk/gatk_indels.vcf.gz
	$(call SELECT_SAMPLE,$^,$@,$1)
endef
$(foreach sample,$(SAMPLES),$(eval $(call SELECT_SAMPLE_INDELS,$(sample))))

gatk/vcf/%.gatk_indels.filtered.vcf.gz : gatk/vcf/%.gatk_indels.vcf.gz
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),,"\
	$(call GATK42,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $(REF_FASTA) \
	$(HAPCALL_INDEL_FILTERS) -V $< -O $@")
	
else # else = no hard filter so use variantrecal

gatk/gatk_indels.recal.vcf.gz : gatk/gatk_indels.vcf.gz
	$(call VARIANT_RECAL_INDELS,$^,$@)

gatk/gatk_indels.filtered.vcf.gz : gatk/gatk_indels.vcf.gz gatk/gatk_indels.recal.vcf.gz
	$(call APPLY_VARIANT_RECAL_INDELS,$<,$(word 2,$^),$@)
endif # end GATK_HARD_FILTER_INDELS



##### END cleaning GATK results with either filtering or recalibration


# filter for only novel snps/indels
%.novel.txt : %.txt
	$(INIT) /bin/awk 'NR == 1 || $$4 == "."' $< > $@

vcf/all$(PROJECT_PREFIX).gatk_snps.vcf : gatk/gatk_snps.filtered.vcf.gz
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),,"\
	$(GUNZIP) $< > $@")
#	$(INIT) $(FIX_GATK_VCF) $< > $@

vcf/all$(PROJECT_PREFIX).gatk_indels.vcf : gatk/gatk_indels.filtered.vcf.gz
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),,"\
	$(GUNZIP) $< > $@")
#	$(INIT) $(FIX_GATK_VCF) $< > $@


reports/%/index.html : reports/%.dp_ft.grp metrics/hs_metrics.txt
	$(call RUN,1,"$(VARIANT_EVAL_GATK_REPORT) --metrics $(word 2,$^) --outDir $(@D) $<")

$(REF_FASTA).fai : $(REF_FASTA)
	$(call RUN,1,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_MEDIUM),$(SAMTOOLS_MODULE),"\
	$(SAMTOOLS) faidx $<")

$(REF_FASTA:.fasta=.dict) : $(REF_FASTA)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	$(call PICARD,CreateSequenceDictionary,$(RESOURCE_REQ_LOW_MEM_JAVA)) REFERENCE=$< OUTPUT=$@")


# merge variants 
include usb-modules-v2/bam_tools/processBam.mk
include usb-modules-v2/vcf_tools/vcftools.mk

endif
GATK_MK = true
