
# common GATK steps for variant calling
# Author: Raymond Lim <raylim@mm.st> & Fong Chun Chan <fongchunchan@gmail.com>
#

ifndef GATK_MK

include usb-modules-v2/Makefile.inc

HAPLOTYPE_CALLER_OPTS = --dbsnp $(DBSNP) -stand_call_conf $(HAPCALL_CALL_THRESHOLD) -R $(REF_FASTA)

###### RECIPES #######

%.intervals : %.vcf
	$(INIT) sed '/^#/d' $< | awk '{print $$1":"$$2 }' > $@

gatk/dbsnp/%.gatk_snps.vcf : bam/%.bam bam/%.bai
	$(call RUN,4,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call GATK,UnifiedGenotyper,$(RESOURCE_REQ_LOW_MEM)) \
	-nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP_TARGETS_INTERVALS) $(foreach bam,$(filter %.bam,$<),-I $(bam) ) \
	--genotyping_mode GENOTYPE_GIVEN_ALLELES -alleles $(DBSNP_TARGETS_INTERVALS) -o $@ --output_mode EMIT_ALL_SITES")

#### Section for Haplotype caller for variant calling
## options are SPLIT_CHR and single vs sample sets

ifeq ($(SPLIT_CHR),true)
ifdef SAMPLE_SET_PAIRS
define hapcall-vcf-sets-chr
gatk/chr_vcf/$1.$2.variants.vcf : $$(foreach sample,$$(samples.$1),gatk/chr_vcf/$$(sample).$2.variants.intervals) $$(foreach sample,$$(samples.$1),bam/$$(sample).bam bam/$$(sample).bai)
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
	$$(call GATK,HaplotypeCaller,$$(RESOURCE_REQ_MEDIUM_MEM)) $$(HAPLOTYPE_CALLER_OPTS) \
	$$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) $$(foreach intervals,$$(filter %.intervals,$$^),-L $$(intervals) ) -o $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call hapcall-vcf-sets-chr,$(set),$(chr)))))

define merge-chr-variants-sets
gatk/vcf/$1.variants.vcf : $$(foreach chr,$$(CHROMOSOMES),gatk/chr_vcf/$1.$$(chr).variants.vcf)
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) \
	--assumeIdenticalSamples $$(foreach i,$$^, --variant $$i) -R $$(REF_FASTA) -o $$@")
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call merge-chr-variants-sets,$(set))))
endif # def SAMPLE_SET_PAIRS

define chr-variants
gatk/chr_vcf/%.$1.variants.vcf : bam/%.bam bam/%.bai
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
	$$(call GATK37,HaplotypeCaller,$$(RESOURCE_REQ_MEDIUM_MEM)) $$(HAPLOTYPE_CALLER_OPTS) \
	-L $1 -I $$< -o $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-variants,$(chr))))

define merge-chr-variants
gatk/vcf/$1.variants.vcf : $$(foreach chr,$$(CHROMOSOMES),gatk/chr_vcf/$1.$$(chr).variants.vcf)
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_MEDIUM_MEM)) \
	--assumeIdenticalSamples $$(foreach i,$$^, --variant $$i) -R $$(REF_FASTA) -o $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call merge-chr-variants,$(sample))))

else #### no splitting by chr ####
## call sample sets
ifdef SAMPLE_SETS
define hapcall-vcf-sets
gatk/vcf/$1.variants.vcf : $$(foreach sample,$$(samples.$1),gatk/vcf/$$(sample).variants.vcf) $$(foreach sample,$$(samples.$1),bam/$$(sample).bam bam/$$(sample).bai)
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
	$$(call GATK,HaplotypeCaller,$$(RESOURCE_REQ_MEDIUM_MEM)) $$(HAPLOTYPE_CALLER_OPTS) \
	$$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) $$(foreach vcf,$$(filter %.vcf,$$^),-L $$(vcf) ) -o $$@")
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call hapcall-vcf-sets,$(set))))
endif

define hapcall-vcf
gatk/vcf/$1.variants.vcf : bam/$1.bam bam/$1.bai
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_MEDIUM),$$(JAVA8_MODULE),"\
	$$(call GATK,HaplotypeCaller,$$(RESOURCE_REQ_MEDIUM_MEM)) $$(HAPLOTYPE_CALLER_OPTS) -I $$< -o $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call hapcall-vcf,$(sample))))

endif # split by chr

####### END HaplotypeCaller

gatk/vcf/%.variants.snps.vcf : gatk/vcf/%.variants.vcf gatk/vcf/%.variants.vcf.idx
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(JAVA8_MODULE),"\
	$(call GATK,SelectVariants,$(RESOURCE_REQ_MEDIUM_MEM)) \
	-R $(REF_FASTA) --variant $<  -o $@ -selectType SNP")

gatk/vcf/%.variants.indels.vcf : gatk/vcf/%.variants.vcf gatk/vcf/%.variants.vcf.idx
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(JAVA8_MODULE),"\
	$(call GATK,SelectVariants,$(RESOURCE_REQ_MEDIUM_MEM)) \
	-R $(REF_FASTA) --variant $<  -o $@ -selectType INDEL")


###### Here are the options for filtering HaplotypeCaller Results
##### Different for snps and indels
# a, use variant Recalibration (recommended)
# b, use hard filters

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

define VARIANT_RECAL_SNPS
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call GATK,VariantRecalibrator,$(RESOURCE_REQ_LOW_MEM)) \
	-R $(REF_FASTA) -nt 1 \
	-resource:hapmap$(,)known=false$(,)training=true$(,)truth=true$(,)prior=15.0 $(HAPMAP) \
	-resource:omni$(,)known=false$(,)training=true$(,)truth=true$(,)prior=12.0 $(OMNI) \
	-resource:1000G$(,)known=false$(,)training=true$(,)truth=false$(,)prior=10.0 $(R1000G) \
	-resource:dbsnp$(,)known=true$(,)training=false$(,)truth=false$(,)prior=2.0 $(DBSNP) \
	$(foreach i,$(VARIANT_RECAL_ANNOTATIONS_SNPS), -an $i) -mode SNP \
	$(foreach i,$(filter %.vcf,$2), -input $i) \
	-recalFile $1 -tranchesFile $(basename $1).tranches")
endef
# -rscriptFile $$(basename $1).snps.plots.R")

define VARIANT_RECAL_INDELS
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(JAVA8_MODULE),"\
	$(call GATK,VariantRecalibrator,$(RESOURCE_REQ_LOW_MEM)) \
	-R $(REF_FASTA) -nt 1 \
	--maxGaussians 4 \
	-resource:mills$(,)known=false$(,)training=true$(,)truth=true$(,)prior=12.0 $(KNOWN_INDELS) \
	-resource:dbsnp$(,)known=true$(,)training=false$(,)truth=false$(,)prior=2.0 $(DBSNP) \
	$(foreach i,$(VARIANT_RECAL_ANNOTATIONS_INDELS), -an $i) -mode INDEL \
	$(foreach i,$(filter %.vcf,$2), -input $i) \
	-recalFile $1 -tranchesFile $(basename $1).tranches")
endef

define APPLY_VARIANT_RECAL_SNPS
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(JAVA8_MODULE),"\
	$(call GATK,ApplyRecalibration,$(RESOURCE_REQ_MEDIUM_MEM)) \
	-R $(REF_FASTA) -input $2 -recalFile $3 -mode SNP \
	--ts_filter_level $(VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL_SNPS) \
	-tranchesFile $(basename $3).tranches -o $1")
endef

define APPLY_VARIANT_RECAL_INDELS
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(JAVA8_MODULE),"\
	$(call GATK,ApplyRecalibration,$(RESOURCE_REQ_MEDIUM_MEM)) \
	-R $(REF_FASTA) -input $2 -recalFile $3 -mode INDEL\
	--ts_filter_level $(VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL_INDELS) \
	-tranchesFile $(basename $3).tranches -o $1")
endef

## recipes for SNPs
ifeq ($(GATK_HARD_FILTER_SNPS),true)
gatk/vcf/%.variants.snps.filtered.vcf : gatk/vcf/%.variants.snps.vcf gatk/vcf/%.variants.snps.vcf.idx
	$(call RUN,4,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(JAVA8_MODULE),"\
	$(call GATK,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM)) -R $(REF_FASTA) \
	$(HAPCALL_SNP_FILTERS) -o $@ --variant $<")
else # else = no hard filter so use variantrecal

# pool sample vcfs for recalibration
ifeq ($(GATK_POOL_VARIANT_RECAL),true)
gatk/vcf/samples.snps.recal.vcf : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.snps.vcf) $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.snps.vcf.idx)
	$(call VARIANT_RECAL_SNPS,$@,$^)

define SAMPLE_APPLY_RECAL_SNPS
gatk/vcf/$1.variants.snps.filtered.vcf : gatk/vcf/$1.variants.snps.vcf gatk/vcf/samples.snps.recal.vcf gatk/vcf/samples.snps.recal.vcf.idx gatk/vcf/$1.variants.snps.vcf.idx 
	$$(call APPLY_VARIANT_RECAL_SNPS,$$@,$$<,$$(word 2,$$^))
endef
$(foreach sample,$(SAMPLES),$(eval $(call SAMPLE_APPLY_RECAL_SNPS,$(sample))))

ifdef SAMPLE_SETS
gatk/vcf/sets.snps.recal.vcf : $(foreach set,$(SAMPLE_SET_PAIRS),gatk/vcf/$(set).variants.snps.vcf gatk/vcf/$(set).variants.snps.vcf.idx)
	$(call VARIANT_RECAL_SNPS,$@,$^)
			
define SETS_APPLY_RECAL_SNPS
gatk/vcf/$1.variants.snps.filtered.vcf : gatk/vcf/$1.variants.snps.vcf gatk/vcf/sets.snps.recal.vcf gatk/vcf/sets.snps.recal.vcf.idx gatk/vcf/$1.variants.snps.vcf.idx 
	$$(call APPLY_VARIANT_RECAL_SNPS,$$@,$$<,$$(word 2,$$^))
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call SETS_APPLY_RECAL_SNPS,$(set))))
endif #end ifdef SAMPLE_SETS

else  #else GATK_POOL_VARIANT_RECAL
gatk/vcf/%.variants.snps.recal.vcf : gatk/vcf/%.variants.snps.vcf gatk/vcf/%.variants.snps.vcf.idx
	$(call VARIANT_RECAL_SNPS,$@,$^)

gatk/vcf/%.variants.snps.filtered.vcf : gatk/vcf/%.variants.snps.vcf gatk/vcf/%.variants.snps.recal.vcf gatk/vcf/%.variants.snps.vcf.idx gatk/vcf/%.variants.snps.recal.vcf.idx
	$(call APPLY_VARIANT_RECAL_SNPS,$@,$<,$(word 2,$^))
endif # end GATK_POOL_SNP_RECAL
endif # end GATK_HARD_FILTER_SNPS


ifeq ($(GATK_HARD_FILTER_INDELS),true)
gatk/vcf/%.variants.indels.filtered.vcf : gatk/vcf/%.variants.indels.vcf gatk/vcf/%.variants.indels.vcf.idx
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(JAVA8_MODULE),"\
	$(call GATK,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM)) -R $(REF_FASTA) \
	$(HAPCALL_INDEL_FILTERS) -o $@ --variant $<")
	
else # else = no hard filter so use variantrecal

# pool sample vcfs for recalibration
ifeq ($(GATK_POOL_VARIANT_RECAL),true)
gatk/vcf/samples.indels.recal.vcf : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.indels.vcf) $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.indels.vcf.idx)
	$(call VARIANT_RECAL_INDELS,$@,$^)

define SAMPLE_APPLY_RECAL_INDELS
gatk/vcf/$1.variants.indels.filtered.vcf : gatk/vcf/$1.variants.indels.vcf gatk/vcf/samples.indels.recal.vcf gatk/vcf/samples.indels.recal.vcf.idx gatk/vcf/$1.variants.indels.vcf.idx 
	$$(call APPLY_VARIANT_RECAL_INDELS,$$@,$$<,$$(word 2,$$^))
endef
$(foreach sample,$(SAMPLES),$(eval $(call SAMPLE_APPLY_RECAL_INDELS,$(sample))))

ifdef SAMPLE_SETS
gatk/vcf/sets.indels.recal.vcf : $(foreach set,$(SAMPLE_SET_PAIRS),gatk/vcf/$(set).variants.indels.vcf gatk/vcf/$(set).variants.indels.vcf.idx)
	$(call VARIANT_RECAL_INDELS,$@,$^)
			
define SETS_APPLY_RECAL_INDELS
gatk/vcf/$1.variants.indels.filtered.vcf : gatk/vcf/$1.variants.indels.vcf gatk/vcf/sets.indels.recal.vcf gatk/vcf/sets.indels.recal.vcf.idx gatk/vcf/$1.variants.indels.vcf.idx 
	$$(call APPLY_VARIANT_RECAL_INDELS,$$@,$$<,$$(word 2,$$^))
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call SETS_APPLY_RECAL_INDELS,$(set))))
endif #end ifdef SAMPLE_SETS

else  #else GATK_POOL_VARIANT_RECAL
gatk/vcf/%.variants.indels.recal.vcf : gatk/vcf/%.variants.indels.vcf gatk/vcf/%.variants.indels.vcf.idx
	$(call VARIANT_RECAL_INDELS,$@,$^)

gatk/vcf/%.variants.indels.filtered.vcf : gatk/vcf/%.variants.indels.vcf gatk/vcf/%.variants.indels.recal.vcf gatk/vcf/%.variants.indels.vcf.idx gatk/vcf/%.variants.indels.recal.vcf.idx
	$(call APPLY_VARIANT_RECAL_INDELS,$@,$<,$(word 2,$^))
endif # end GATK_POOL_VARIANT_RECAL
endif # end GATK_HARD_FILTER_INDELS



##### END cleaning GATK results with either filtering or recalibration


# filter for only novel snps/indels
%.novel.txt : %.txt
	$(INIT) /bin/awk 'NR == 1 || $$4 == "."' $< > $@

vcf/%.gatk_snps.vcf : gatk/vcf/%.variants.snps.filtered.vcf
	$(INIT) $(FIX_GATK_VCF) $< > $@

vcf/%.gatk_indels.vcf : gatk/vcf/%.variants.indels.filtered.vcf
	$(INIT) $(FIX_GATK_VCF) $< > $@


reports/%/index.html : reports/%.dp_ft.grp metrics/hs_metrics.txt
	$(call RUN,1,"$(VARIANT_EVAL_GATK_REPORT) --metrics $(word 2,$^) --outDir $(@D) $<")

$(REF_FASTA).fai : $(REF_FASTA)
	$(call RUN,1,$(RESOURCE_REQ_LOWMEM),$(RESOURCE_REQ_MEDIUM),$(SAMTOOLS_MODULE),"\
	$(SAMTOOLS) faidx $<")

$(REF_FASTA:.fasta=.dict) : $(REF_FASTA)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT)),$(JAVA8_MODULE),"\
	$(call PICARD,CreateSequenceDictionary,$(RESOURCE_REQ_LOW_MEM)) REFERENCE=$< OUTPUT=$@")


# merge variants 
include usb-modules-v2/bam_tools/processBam.mk
include usb-modules-v2/vcf_tools/vcftools.mk

endif
GATK_MK = true
