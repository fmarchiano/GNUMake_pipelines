ifndef VCFTOOLS_MK_NONTVC

%.dp_ft.vcf : %.vcf
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
		$(call GATK,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
		-R $(REF_FASTA) -V $< -o $@ --filterExpression \
		'! vc.hasAttribute(\"DP\") || DP < $(MIN_NORMAL_DEPTH) || DP == \".\"' --filterName DepthO")

%.lowFreq_ft.vcf : %.vcf
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
		$(call GATK,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $$(REF_FASTA) -V $< -o $@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 * 1.0 > 0' --filterName interrogation \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 * 1.0 == 0' --filterName interrogation_Absent")

#%.pass.vcf : %.vcf
#	$(call CHECK_VCF,$<,$@,\
#	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
#		$(VCF_PASS) $< $@ $(subst pass,fail,$@)"))

ifdef SAMPLE_PAIRS
define som-ad-ft-tumor-normal
vcf/$1_$2.%.som_ad_ft.vcf : vcf/$1_$2.%.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 < $(MIN_TUMOR_AD)' \
		--filterName tumorVarAlleleDepth \
		--filterExpression 'if (vc.getGenotype(\"$2\").getDP() > $(MIN_NORMAL_DEPTH)) { \
		( vc.getGenotype(\"$2\").getAD().1 * 1.0 / vc.getGenotype(\"$2\").getDP()) > \
			(vc.getGenotype(\"$1\").getAD().1 * 1.0 / vc.getGenotype(\"$1\").getDP()) / $(MIN_TN_AD_RATIO) \
			} else { ( vc.getGenotype(\"$2\").getAD().1 * 1.0 / vc.getGenotype(\"$2\").getDP()) > \
			(vc.getGenotype(\"$1\").getAD().1 * 1.0 / vc.getGenotype(\"$1\").getDP()) / $(MIN_TN_AD_RATIO) && \
			vc.getGenotype(\"$1\").getAD().1 * 1.0 < 1 && vc.getGenotype(\"$2\").getAD().1 > 1 }' \
		--filterName somaticAlleleDepth \
		--filterExpression 'vc.getGenotype(\"$1\").getDP() <= $$(MIN_TUMOR_DEPTH) || vc.getGenotype(\"$2\").getDP() <= $$(MIN_NORMAL_DEPTH)' \
		--filterName depthFilter && sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call som-ad-ft-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define som-ad-ft-tumor-normal_IGR
vcf/$1_$2.%.som_ad_ft_igr.vcf : vcf/$1_$2.%.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $$(REF_FASTA) -V $$< -o $$@ \
			--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"AC\") * 1.0 < $(MIN_TUMOR_AD)' \
			--filterName tumorVarAlleleDepth \
			--filterExpression 'if ((vc.getGenotype(\"$2\").getAnyAttribute(\"AC\") * 1.0 + \
				vc.getGenotype(\"$2\").getAnyAttribute(\"RC\") * 1.0) > $(MIN_NORMAL_DEPTH)) { \
					(vc.getGenotype(\"$2\").getAnyAttribute(\"AC\") * 1.0 / \
					(vc.getGenotype(\"$2\").getAnyAttribute(\"AC\") * 1.0 + vc.getGenotype(\"$2\").getAnyAttribute(\"RC\") * 1.0)) > \
					(vc.getGenotype(\"$1\").getAnyAttribute(\"AC\") * 1.0 * 1.0 / \
					(vc.getGenotype(\"$1\").getAnyAttribute(\"AC\") * 1.0 + vc.getGenotype(\"$1\").getAnyAttribute(\"RC\") * 1.0)) / $(MIN_TN_AD_RATIO) \
				} else { ( vc.getGenotype(\"$2\").getAnyAttribute(\"AC\") * 1.0 / \
					(vc.getGenotype(\"$2\").getAnyAttribute(\"AC\") * 1.0 + vc.getGenotype(\"$2\").getAnyAttribute(\"RC\") * 1.0)) > \
					(vc.getGenotype(\"$1\").getAnyAttribute(\"AC\") * 1.0 / \
					(vc.getGenotype(\"$1\").getAnyAttribute(\"AC\") * 1.0 + vc.getGenotype(\"$1\").getAnyAttribute(\"RC\") * 1.0)) / $(MIN_TN_AD_RATIO) && \
					vc.getGenotype(\"$1\").getAnyAttribute(\"AC\") * 1.0 < 1 && \
					vc.getGenotype(\"$2\").getAnyAttribute(\"AC\") * 1.0 > 1 }' \
			--filterName somaticAlleleDepth \
			--filterExpression '(vc.getGenotype(\"$1\").getAnyAttribute(\"AC\") * 1.0 + \
				vc.getGenotype(\"$1\").getAnyAttribute(\"RC\") * 1.0) <= $$(MIN_TUMOR_DEPTH) || \
				(vc.getGenotype(\"$2\").getAnyAttribute(\"AC\") * 1.0 + \
				vc.getGenotype(\"$2\").getAnyAttribute(\"RC\") * 1.0) <= $$(MIN_NORMAL_DEPTH)' \
			--filterName depthFilter && sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call som-ad-ft-tumor-normal_IGR,$(tumor.$(pair)),$(normal.$(pair)))))

define sufam
ifeq ($3,$2)
vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.vcf
	$$(INIT) ln -f $$< $$@
else
# Note, if sufam.tmp.tmp is empty, <grep -v '^#'> will finish with exit code 1, which will make run.py stop.
# Using <$$(call CHECK_VCF_CMD,$$@.tmp,cp $$@.tmp $$@,> or hardcoding the same check in the code below works in some cases, but breaks when sample sets are larger.
# The error complains about grep (from the conditional) not finding $$@.tmp. But $$@.tmp should be created (and the log says it is). I don't get it.
# An alternative to hide grep's exit code 1 is to wrap it in a <{grep ... || true; }> construct (https://stackoverflow.com/questions/6550484/prevent-grep-returning-an-error-when-input-doesnt-match)
vcf/$3.%.sufam.tmp : $$(foreach tumor,$$(wordlist 1,$$(shell expr $$(words $$(subst _,$$(space),$3)) - 1),$$(subst _,$$(space),$3)),vcf/$$(tumor)_$$(lastword $$(subst _,$$(space),$3)).%.vcf)
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
		$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGH_MEM_JAVA)) \
		$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@.tmp --genotypemergeoption UNSORTED \
		-sites_only -R $$(REF_FASTA) && \
		grep \"#\" $$@.tmp > $$@ && \
		 { grep -v \"#\" $$@.tmp || true; } | awk '{OFS=\"\t\"; \$$$$7=\"PASS\" ; print ;}' >> $$@")

ifeq ($$(findstring varscan,$$(MUT_CALLER)),varscan)
vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.vcf vcf/$3.%.sufam.tmp bam/$1.bam bam/$2.bam
	$$(call CHECK_VCF_CMD,$$(word 2,$$^),cp $$(word 1,$$^) $$@,\
		$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(JAVA8_MODULE) $$(SAMTOOLS_MODULE),"\
		$$(call GATK,SelectVariants,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $$(REF_FASTA) \
		--variant $$(word 2,$$^) --discordance $$(word 1,$$^) -sites_only | grep -v \"#\" | \
		awk 'BEGIN { OFS=\"\t\"} { print \$$$$1$$(,)\$$$$2-10$$(,)\$$$$2+10 }' > $$@.tmp1 && \
		$$(SAMTOOLS) mpileup -l $$@.tmp1 -d 10000 -q $$(MIN_MQ) -f $$(REF_FASTA) $$(word 4,$$^) $$(word 3,$$^) | \
		$$(VARSCAN) somatic - $$@.tmp2 --output-vcf 1 --mpileup 1 --validation && \
		$$(FIX_VARSCAN_VCF) -t $1 -n $2 $$@.tmp2.validation > $$@.tmp3 && \
		$$(call GATK,SelectVariants,$$(RESOURCE_REQ_HIGH_MEM_JAVA)) -R $$(REF_FASTA) \
		--variant $$@.tmp3 --concordance $$(word 2,$$^) > $$@.tmp4 && \
		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $$(REF_FASTA) -V $$@.tmp4 -o $$@.tmp5 \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 * 1.0 > 0' --filterName interrogation && \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 * 1.0 == 0' --filterName interrogation_Absent && \
		$$(call PURGE_AND_LOAD, $$(SNP_EFF_MODULE)) && \
		$$(call SNP_SIFT,$$(RESOURCE_REQ_LOW_MEM_JAVA)) filter $$(SNP_SIFT_OPTS) -f $$@.tmp5 \
		\"(FILTER has 'interrogation') | (FILTER has 'interrogation_Absent')\" > $$@.tmp6 && \
		$$(call PURGE_AND_LOAD, $$(JAVA8_MODULE)) && \
		$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGH_MEM_JAVA)) --variant $$< --variant $$@.tmp6 -o $$@ \
		--genotypemergeoption UNSORTED -R $$(REF_FASTA) && \
		$$(RM) $$@.tmp* $$(word 1,$$^)*.idx $$(word 2,$$^)* $$(word 2,$$^)*.idx"))
else
vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.vcf vcf/$3.%.sufam.tmp bam/$1.bam bam/$2.bam
	$$(call CHECK_VCF_CMD,$$(word 2,$$^),cp $$(word 1,$$^) $$@,\
		$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(JAVA8_MODULE),"\
		$$(call GATK,SelectVariants,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $$(REF_FASTA) \
		--variant $$(word 2,$$^) --discordance $$(word 1,$$^) -sites_only -o $$@.tmp1 && \
		$$(call GATK4130,HaplotypeCaller,$$(RESOURCE_REQ_HIGH_MEM_JAVA)) -VS LENIENT -R $$(REF_FASTA) -I $$(word 3,$$^) -I $$(word 4,$$^) \
		--dbsnp $$(DBSNP_TARGETS_INTERVALS) -L $$@.tmp1 \
		--genotyping-mode GENOTYPE_GIVEN_ALLELES --output-mode EMIT_ALL_SITES --alleles $$@.tmp1 -O $$@.tmp2 && \
		$$(FIX_GATK_VCF) $$@.tmp2 | sed -e 's/\t\*$$(,)/\t/g' -e 's/$$(,)\*\t/\t/g' > $$@.tmp2fixed && \
		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $$(REF_FASTA) -V $$@.tmp2fixed -o $$@.tmp3 \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 * 1.0 > 0' --filterName interrogation \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 * 1.0 == 0' --filterName interrogation_Absent && \
		$$(call PURGE_AND_LOAD, $$(SNP_EFF_MODULE)) && \
		$$(call SNP_SIFT,$$(RESOURCE_REQ_LOW_MEM_JAVA)) filter $$(SNP_SIFT_OPTS) -f $$@.tmp3 \
		\"(FILTER has 'interrogation') | (FILTER has 'interrogation_Absent')\" > $$@.tmp4 && \
		$$(call PURGE_AND_LOAD, $$(JAVA8_MODULE)) && \
		$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGH_MEM_JAVA)) --variant $$< --variant $$@.tmp4 -o $$@ \
		--genotypemergeoption UNSORTED -R $$(REF_FASTA) && \
		$$(RM) $$@.tmp* $$(word 1,$$^)*.idx $$(word 2,$$^)* $$(word 2,$$^)*.idx"))

# && \
#		$$(RM) $$@.tmp1 $$@.tmp2 $$@.tmp3 $$(word 2,$$^) $$@.tmp1.idx $$@.tmp2.idx $$@.tmp3.idx $$(word 2,$$^).idx \
#		$$@.tmp2fixed $$@.tmp2fixed.idx"))
endif
endif #ifeq ($3,$2)
endef #define sufam

$(foreach set,$(SAMPLE_SETS),\
	$(foreach tumor,$(wordlist 1,$(shell expr $(words $(subst _,$(space),$(set))) - 1),$(subst _,$(space),$(set))),\
		$(eval $(call sufam,$(tumor),$(lastword $(subst _,$(space),$(set))),$(subst $(tumor)_,,$(set))))))

include usb-modules-v2/variant_callers/somatic/pon.mk

endif #ifdef SAMPLE_PAIRS

ifeq ($(findstring true,$(TUMOR_ONLY)),true)
define som-ad-ft-tumor-only
vcf/$1.%.som_ad_ft.vcf : vcf/$1.%.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 < $(MIN_TUMOR_AD)' \
		--filterName tumorVarAlleleDepth \
		--filterExpression 'vc.getGenotype(\"$1\").getDP() <= $$(MIN_TUMOR_DEPTH)' \
		--filterName depthFilter && sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach sample,$(SAMPLES),$(eval $(call som-ad-ft-tumor-only,$(sample))))
endif

VCFTOOLS_MK_NONTVC = true
endif
