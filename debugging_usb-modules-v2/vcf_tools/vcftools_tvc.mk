ifndef VCFTOOLS_MK_TVC



%.dp_ft.vcf : %.vcf
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
		$(call GATK,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
		-R $(REF_FASTA) -V $< -o $@ \
		--filterExpression 'DP <= $(MIN_NORMAL_DEPTH) || FDP <= $(MIN_NORMAL_DEPTH)' \
		--filterName Depth && $(RM) $< $<.idx")

#%.pass.vcf : %.vcf
#	$(call CHECK_VCF,$<,$@,\
#	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(SNP_EFF_MODULE),"\
#		$(call SNP_SIFT,$(RESOURCE_REQ_LOW_MEM_JAVA)) filter $(SNP_SIFT_OPTS) -f $< \
#		\"( na FILTER ) | (FILTER = 'PASS') | (FILTER has 'HOTSPOT')\" > $@"))

ifdef SAMPLE_PAIRS
define som-ad-ft-tumor-normal
vcf/$1_$2.%.som_ad_ft.vcf : vcf/$1_$2.%.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") < $(MIN_TUMOR_AD)' \
		--filterName tumorVarAD_flow \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") < $(MIN_TUMOR_AD)' \
		--filterName tumorVarAD_raw \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"DP\") * 1.0 <= $(MIN_TUMOR_DEPTH)' \
		--filterName tumorDepthFilter_raw \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FDP\") * 1.0 <= $(MIN_TUMOR_DEPTH)' \
		--filterName tumorDepthFilter_flow \
		--filterExpression 'vc.getGenotype(\"$2\").getAnyAttribute(\"DP\") * 1.0 <= $(MIN_NORMAL_DEPTH)' \
		--filterName normalDepthFilter_raw \
		--filterExpression '( vc.getGenotype(\"$2\").getAnyAttribute(\"AO\") * 1.0 / vc.getGenotype(\"$2\").getAnyAttribute(\"DP\") * 1.0  ) \
			> ( vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") * 1.0 / vc.getGenotype(\"$1\").getAnyAttribute(\"DP\") * 1.0  ) / $(MIN_TN_AD_RATIO)' \
		--filterName somAD_raw \
		--filterExpression '( vc.getGenotype(\"$2\").getAnyAttribute(\"AO\") * 1.0 / vc.getGenotype(\"$2\").getAnyAttribute(\"DP\") * 1.0  ) \
			> ( vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") * 1.0 / vc.getGenotype(\"$1\").getAnyAttribute(\"FDP\") * 1.0  ) / $(MIN_TN_AD_RATIO)' \
		--filterName somAD_flow \
		&& sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call som-ad-ft-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
#               --filterExpression 'QUAL < 60' \
#               --filterName lowQual \

define sufam
ifeq ($3,$2)
vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.vcf
	$$(INIT) ln -f $$< $$@
else

vcf/$3.%.sufam.tmp1.vcf : $$(foreach tumor,$$(wordlist 1,$$(shell expr $$(words $$(subst _,$$(space),$3)) - 1),$$(subst _,$$(space),$3)),vcf/$$(tumor)_$$(lastword $$(subst _,$$(space),$3)).%.vcf)
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_VSHORT),$$(JAVA8_MODULE),"\
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGH_MEM_JAVA)) \
	$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED -R $$(REF_FASTA)")

vcf/$1_$2.%/isec/0001.vcf : vcf/$1_$2.%.vcf vcf/$3.%.sufam.tmp1.vcf vcf/$1_$2.%.vcf.gz vcf/$3.%.sufam.tmp1.vcf.gz vcf/$1_$2.%.vcf.gz.tbi vcf/$3.%.sufam.tmp1.vcf.gz.tbi
	$(MKDIR) $$(dir $$@); $$(call CHECK_VCF,$$(word 2,$$^),$$@,\
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE),"\
	$$(BCFTOOLS) isec -O v -p $$(dir $$@) $$(word 3,$$^) $$(word 4,$$^)"))

vcf/$1_$2.%/TSVC_variants.vcf.gz : vcf/$1_$2.%/isec/0001.sorted.vcf bam/$1.bam bam/$1.bam.bai
	$$(call RUN,$$(TVC_NUM_CORES),$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$(OPENBLAS_MODULE) $(PYTHON_MODULE),"\
	$$(TVC) -s $$< -i $$(word 2,$$^) -r $$(REF_FASTA) -o $$(@D) -N $$(TVC_NUM_CORES) \
	$$(if $$(TARGETS_FILE_INTERVALS),-b $$(TARGETS_FILE_INTERVALS)) -m $$(TVC_MOTIF) \
	-t $$(TVC_ROOT_DIR) --primer-trim-bed $$(PRIMER_TRIM_BED) -p $(TVC_SENSITIVE_JSON)")

vcf/$1_$2.%/isec2/0002.vcf : vcf/$1_$2.%/TSVC_variants.biallelic_ft.vcf vcf/$1_$2.%/isec/0001.sorted.vcf vcf/$1_$2.%/TSVC_variants.biallelic_ft.vcf.gz vcf/$1_$2.%/isec/0001.sorted.vcf.gz vcf/$1_$2.%/TSVC_variants.biallelic_ft.vcf.gz.tbi vcf/$1_$2.%/isec/0001.sorted.vcf.gz.tbi
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE),"\
	$$(BCFTOOLS) isec -O v -p $$(dir $$@) $$(word 3,$$^) $$(word 4,$$^) && \
	grep -v \"##contig\" $$@ > tmp && mv tmp $$@" )

vcf/$1_$2.%/isec2/0001.vcf : vcf/$1_$2.%/isec2/0002.vcf
	

vcf/$1_$2.%/isec3/0002.vcf : vcf/$1_$2.%/TSVC_variants.multiallelic_ft.norm.left_align.vcf vcf/$1_$2.%/isec2/0001.vcf vcf/$1_$2.%/TSVC_variants.multiallelic_ft.norm.left_align.vcf.gz vcf/$1_$2.%/isec2/0001.vcf.gz vcf/$1_$2.%/TSVC_variants.multiallelic_ft.norm.left_align.vcf.gz.tbi vcf/$1_$2.%/isec2/0001.vcf.gz.tbi
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE),"\
	$$(BCFTOOLS) isec -O v -p $$(dir $$@) $$(word 3,$$^) $$(word 4,$$^)")

vcf/$1_$2.%.sufam.tmp2.vcf : vcf/$1_$2.%/isec2/0002.post_bcftools.vcf vcf/$1_$2.%/isec3/0002.post_bcftools.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(JAVA8_MODULE),"\
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_LOW_MEM_JAVA)) -R $$(REF_FASTA) \
	$$(foreach vcf,$$^,-V $$(vcf)) -o $$@ --assumeIdenticalSamples && $$(RMR) $$(dir $$<)")

vcf/$1_$2.%.sufam.tmp3.vcf : vcf/$1_$2.%.sufam.tmp2.vcf
	$$(call CHECK_VCF,$$<,$$@,\
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(JAVA8_MODULE),"\
	$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_LOW_MEM_JAVA)) -R $$(REF_FASTA) -V $$<  -o $$@  \
	--filterName interrogation \
	--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") > 0 || vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") > 0' \
	--filterName interrogation_Absent \
	--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"FAO\") == 0 && vc.getGenotype(\"$1\").getAnyAttribute(\"AO\") == 0' \
	--maskName HOTSPOT --mask $(CANCER_HOTSPOT_VCF) \
	&& $$(RM) $$< $$<.idx"))

vcf/$1_$2.%.sufam.tmp4.vcf : vcf/$3.%.sufam.tmp1.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(JAVA8_MODULE),"\
	$$(call GATK,SelectVariants,$$(RESOURCE_REQ_LOW_MEM_JAVA)) -R $$(REF_FASTA) -V $$< -o $$@ -sn $2")

vcf/$1_$2.%.sufam.vcf : vcf/$1_$2.%.sufam.tmp3.vcf vcf/$1_$2.%.sufam.tmp4.vcf vcf/$1_$2.%.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(JAVA8_MODULE),"\
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_LOW_MEM_JAVA)) -R $$(REF_FASTA) \
	$$(foreach vcf,$$^,-V $$(vcf)) -o $$@ --genotypemergeoption UNSORTED \
	-filteredRecordsMergeType KEEP_IF_ALL_UNFILTERED && $(RM) vcf/$1_$2.%.sufam.tmp*")

endif
endef

$(foreach set,$(SAMPLE_SETS),\
	$(foreach tumor,$(wordlist 1,$(shell expr $(words $(subst _,$(space),$(set))) - 1),$(subst _,$(space),$(set))),\
		$(eval $(call sufam,$(tumor),$(lastword $(subst _,$(space),$(set))),$(subst $(tumor)_,,$(set))))))

include usb-modules-v2/variant_callers/somatic/pon.mk
endif

VCFTOOLS_MK_TVC = true
endif
