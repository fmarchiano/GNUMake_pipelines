# vim: set ft=make :
# sub module containing vcf related tools

ifndef VCFTOOLS_MK

include usb-modules-v2/Makefile.inc

LOGDIR ?= log/vcf.$(NOW)
..DUMMY := $(shell mkdir -p version; echo "$(SNP_EFF) &> version/snp_eff.txt")

######### GZ & INDEX #####
%.vcf.idx : %.vcf
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	sleep 5 && $(RM) $@ && $(IGVTOOLS) index $< && sleep 5")

%.vcf.gz : %.vcf
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(TABIX_MODULE),"\
	sleep 5 && $(BGZIP) -c -f $< >$@ && sleep 5")

%.vcf : %.vcf.gz
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(TABIX_MODULE),"\
	sleep 5 && zcat $< > $@ && sleep 5")

%.vcf.gz.csi : %.vcf.gz
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BCFTOOLS_MODULE),"\
	sleep 5 && $(BCFTOOLS) index $<")

%.vcf.gz.tbi : %.vcf.gz
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(TABIX_MODULE),"\
	sleep 5 && $(TABIX) -f $< && sleep 5")

############ FILTERS #########

%.target_ft.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call GATK4241,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $(REF_FASTA) -V $< -O $@ \
	--mask $(TARGETS_FILE_INTERVALS) --mask-name targetInterval \
	--filter-not-in-mask && $(RM) $< $<.idx"))

%.het_ft.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call GATK4241,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $(REF_FASTA) -V $< -O $@ \
	--genotype-filter-expression 'isHet == 1' --genotype-filter-name 'Heterozygous positions'"))

%.biallelic_ft.vcf : %.vcf.gz
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BCFTOOLS_MODULE),"\
	$(BCFTOOLS) view -M2 $< | grep -v \"##contig\" > $@"))

%.multiallelic_ft.vcf : %.vcf.gz
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BCFTOOLS_MODULE),"\
	$(BCFTOOLS) view -m3 $< | grep -v \"##contig\" > $@"))

%.sdp_ft.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(SNP_EFF_MODULE),"\
		$(call SNP_SIFT,$(RESOURCE_REQ_LOW_MEM_JAVA)) filter $(SNP_SIFT_OPTS) -f $< '(exists GEN[?].DP) & (GEN[?].DP > 20)' > $@"))

%.strelka_ft.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call GATK4241,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $(REF_FASTA) -V $< -O $@ \
	--filter-expression 'QSI_NT < 30' --filter-name QSI_ref \
	--filter-expression 'IHP > 14' --filter-name iHpol \
	--filter-expression 'MQ0 > 1' --filter-name MQ0 && $(RM) $<"))

ifdef SAMPLE_PAIRS
define oxog-pair
vcf/$1_$2.%.oxog_ft.vcf : vcf/$1_$2.%.vcf metrics/$1.artifact_metrics.pre_adapter_detail_metrics
	$$(call CHECK_VCF,$$<,$$@,\
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
	$$(call GATK4141,FilterByOrientationBias,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $$(REF_FASTA) -V $$< -O $$@ \
	--artifact-modes 'G/T' -P $$(<<)"))
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call oxog-pair,$(tumor.$(pair)),$(normal.$(pair)))))
endif

define oxog-sample
vcf/$1.%.oxog_ft.vcf : vcf/$1.%.vcf metrics/$1.artifact_metrics.pre_adapter_detail_metrics
	$$(call CHECK_VCF,$$<,$$@,\
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
	$$(call GATK4141,FilterByOrientationBias,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $$(REF_FASTA) -V $$< -O $$@ \
	--artifact-modes 'G/T' -P $$(<<)"))
endef
$(foreach sample,$(SAMPLES),$(eval $(call oxog-sample,$(sample))))

ifdef SAMPLE_PAIRS
define ffpe-pair
vcf/$1_$2.%.ffpe_ft.vcf : vcf/$1_$2.%.vcf metrics/$1.artifact_metrics.pre_adapter_detail_metrics
	$$(call CHECK_VCF,$$<,$$@,\
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
	$$(call GATK4141,FilterByOrientationBias,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $$(REF_FASTA) -V $$< -O $$@ \
	--artifact-modes 'C/T' -P $$(<<)"))
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call ffpe-pair,$(tumor.$(pair)),$(normal.$(pair)))))
endif

define ffpe-sample
vcf/$1.%.ffpe_ft.vcf : vcf/$1.%.vcf metrics/$1.artifact_metrics.pre_adapter_detail_metrics
	$$(call CHECK_VCF,$$<,$$@,\
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
	$$(call GATK4141,FilterByOrientationBias,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $$(REF_FASTA) -V $$< -O $$@ \
	--artifact-modes 'C/T' -P $$(<<)"))
endef
$(foreach sample,$(SAMPLES),$(eval $(call ffpe-sample,$(sample))))
	
comma := ,
%.hotspot.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(SNP_EFF_MODULE) $(R4_MODULE),"\
	$(call SNP_EFF,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) ann $(SNP_EFF_OPTS) $(SNP_EFF_GENOME) $< | $(VCF_EFF_ONE_PER_LINE) > $<.tmp-hotspot.eff.OnePerLine.vcf &&\
	$(call SNP_SIFT,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) extractFields -s '|' -e '.' $<.tmp-hotspot.eff.OnePerLine.vcf CHROM POS ID REF ALT QUAL FILTER 'ANN[*].EFFECT' 'ANN[*].GENE' 'ANN[*].HGVS_C' 'ANN[*].AA' > $<.tmp-hotspot.eff.OnePerLine.tab &&\
	$(R4SCRIPT) $(HOTSPOTS_DIR)/prepareHotspotVCFs.R --mutations $<.tmp-hotspot.eff.OnePerLine.tab --tumor_ID $<.tmp-hotspot.eff.OnePerLine --hotspots $(HOTSPOTS_DIR)/hotspots_single_aa_change.$(REF).txt --noncoding_hotspots $(HOTSPOTS_DIR)/hotspots_non_coding.$(REF).txt --splicesite_hotspots $(HOTSPOTS_DIR)/hotspots_splice.$(REF).txt --indel_hotspots $(HOTSPOTS_DIR)/hotspots_indel.txt --skip_3D $(SKIP_3D_HOTSPOTS) --save_table no &&\
	echo -e '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO' > $<.tmp-hotspot.eff.OnePerLine_hotspot_merged.vcf &&\
	cat $<.tmp-hotspot.eff.OnePerLine_hotspot.vcf $<.tmp-hotspot.eff.OnePerLine_hotspot_indel.vcf $<.tmp-hotspot.eff.OnePerLine_hotspot_noncoding.vcf $<.tmp-hotspot.eff.OnePerLine_hotspot_sp.vcf | sed '/#/d' | sort -k1$(comma)1V -k2$(comma)2n >> $<.tmp-hotspot.eff.OnePerLine_hotspot_merged.vcf &&\
	$(call SNP_SIFT,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) annotate $<.tmp-hotspot.eff.OnePerLine_hotspot_merged.vcf $< > $<.tmp-hotspot.eff.OnePerLine_hotspot_merged.annotate.vcf &&\
	line=\`grep -n '##INFO=<ID=' $<.tmp-hotspot.eff.OnePerLine_hotspot_merged.annotate.vcf | tail -n1 | cut -f1 -d:\` &&\
	head -n \$$line $<.tmp-hotspot.eff.OnePerLine_hotspot_merged.annotate.vcf > $@ &&\
	cat $(HOTSPOTS_DIR)/vcf_info_header.txt >> $@ &&\
	tail -n+\`echo \$$line + 1 | bc -l\` $<.tmp-hotspot.eff.OnePerLine_hotspot_merged.annotate.vcf >> $@ &&\
	$(RM) $< $<.idx $<.tmp-hotspot*"))

%.nft.vcf : %.vcf $(PON_VCF)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
		$(call GATK4241,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
		-R $(REF_FASTA) -V $< -O $@ --mask-name 'PoN' --mask $(word 2,$^) && $(RM) $< $<.idx")

%.pass.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
		$(VCF_PASS) -n $(VCF_PASS_MAX_FILTERS) $< $@ $(subst pass,fail,$@) &&\
		$(RM) $< $<.idx"))

%.SVpass.vcf : %.vcf
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BCFTOOLS_MODULE),"\
	$(BCFTOOLS) view -f PASS -o $@ $<")

## This is definitely broken
# somatic filter for structural variants
#vcf/$1_$2.%.sv_som_ft.vcf : vcf/$1_$2.%.vcf
#	$(call CHECK_VCF,$<,$@,\
#	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
#		$$(call GATK,VariantFiltration,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $$(REF_FASTA) -V $$< -o $$@ \
#		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"SU\") <= $$(DEPTH_FILTER)' \
#		--filterName svSupport \
#		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"SU\") < 5 * vc.getGenotype(\"$2\").getAnyAttribute(\"SU\")' \
#		--filterName somaticSvSupport \
#		&& sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@ && $$(RM) $$< $$<.idx"))

#%.altad_ft.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),"$(LOAD_JAVA8_MODULE); \
#		$(call GATK,VariantFiltration,$(RESOURCE_REQ_MEDIUM_MEM)) \
#		-R $(REF_FASTA) -V $< -o $@ --filterExpression 'vc.getGenotype(\"$1\").getAD().1 < 0' --filterName nonZeroAD && $(RM) $< $<.idx")

#%.mapq_ft.vcf : %.vcf.gz %.vcf.gz.tbi $(foreach sample,$(SAMPLES),bam/$(sample).bam)
#	$(call 

#%.encode_ft.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,8G,01:59:59,"$(LOAD_JAVA8_MODULE); $(call VARIANT_FILTRATION,7G) \
#		-R $(REF_FASTA) -V $< -o $@ --maskName 'encode' --mask $(ENCODE_BED) && $(RM) $< $<.idx")



##### PROCESSING #######

%.norm.vcf : %.vcf.gz %.vcf.gz.tbi
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BCFTOOLS_MODULE),"\
	$(BCFTOOLS) norm -m -both $< | grep -v \"##contig\" > $@"))

%.norm2.vcf.gz : %.vcf
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BCFTOOLS_MODULE),"\
	$(BCFTOOLS) norm -f $(REF_FASTA) -m -both $< -Oz -o $@")

%.left_align.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(JAVA8_MODULE),"\
	$(call GATK4241,LeftAlignAndTrimVariants,$(RESOURCE_REQ_LOW_MEM_JAVA)) -R $(REF_FASTA) -V $< -o $@"))

%.post_bcftools.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	grep -v '##contig' $< | $(VCF_SORT) $(REF_DICT) - > $@"))

%.sorted.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	$(VCF_SORT) $(REF_DICT) $< > $@ && $(RM) $<"))

############ ANNOTATION #########

%.eff.vcf : %.vcf %.vcf.idx
	$(call CHECK_VCF,$<,$@,\
		$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(SNP_EFF_MODULE),"\
		$(call SNP_EFF,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) ann $(SNP_EFF_OPTS) -noStats $(SNP_EFF_GENOME) $< > $@ && $(RM) $^ $<*"))

%.nsfp.vcf : %.vcf %.vcf.idx
	$(call CHECK_VCF,$<,$@,\
		$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(SNP_EFF_MODULE),"\
		$(call SNP_SIFT,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) dbnsfp $(SNP_SIFT_OPTS) -db $(DB_NSFP) \
		$(if $(findstring hg38,$(REF)),-f 'rs_dbSNP$(,)hg19_chr$(,)hg19_pos(1-based)$(,)HGVSc_VEP$(,)HGVSp_VEP$(,)VEP_canonical$(,)clinvar_id$(,)clinvar_clnsig$(,)clinvar_trait$(,)MutPred2_score$(,)MutPred2_rankscore$(,)MutPred2_pred$(,)MutPred2_top5_mechanisms$(,)MVP_score$(,)MVP_rankscore$(,)gMVP_score$(,)gMVP_rankscore$(,)MisFit_D_score$(,)AlphaMissense_score$(,)AlphaMissense_rankscore$(,)AlphaMissense_pred$(,)PHACTboost_score$(,)PHACTboost_rankscore$(,)Aloft_Fraction_transcripts_affected$(,)Aloft_prob_Tolerant$(,)Aloft_prob_Recessive$(,)Aloft_prob_Dominant$(,)Aloft_pred$(,)Aloft_Confidence$(,)CADD_raw$(,)CADD_raw_rankscore$(,)CADD_phred$(,)fathmm-XF_coding_score$(,)fathmm-XF_coding_rankscore$(,)fathmm-XF_coding_pred$(,)AllofUs_ALL_AC$(,)AllofUs_ALL_AN$(,)AllofUs_ALL_AF$(,)AllofUs_POPMAX_AF$(,)AllofUs_POPMAX_AC$(,)AllofUs_POPMAX_AN$(,)AllofUs_POPMAX_POP$(,)gnomAD4.1_joint_flag$(,)gnomAD4.1_joint_AC$(,)gnomAD4.1_joint_AN$(,)gnomAD4.1_joint_AF$(,)gnomAD4.1_joint_nhomalt$(,)gnomAD4.1_joint_POPMAX_AC$(,)gnomAD4.1_joint_POPMAX_AN$(,)gnomAD4.1_joint_POPMAX_AF$(,)gnomAD4.1_joint_POPMAX_nhomalt$(,)ALFA_Total_AC$(,)ALFA_Total_AN$(,)ALFA_Total_AF$(,)dbNSFP_POPMAX_AF$(,)dbNSFP_POPMAX_AC$(,)dbNSFP_POPMAX_POP',,) \
		$< | sed '/^##INFO=<ID=dbNSFP/ s/Character/String/' \
		> $@ && $(RM) $^"))

#%.gatk_eff.vcf : %.vcf %.vcf.idx
#	$(call CHECK_VCF,$<,$@,\
#		$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(SNP_EFF_MODULE),"\
#		$(call SNP_EFF,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) eff -i vcf -o gatk $(SNP_EFF_GENOME) $< > $@"))

%.annotated.vcf : %.vcf %.gatk_eff.vcf %.gatk_eff.vcf.idx %.vcf.idx 
	$(call RUN,4,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(JAVA8_MODULE),"\ 
		$(call GATK4241,VariantFiltration,$(RESOURCE_REQ_LOWMEM)) \
		-R $(REF_FASTA) -nt 4 -A SnpEff --variant $< --snpEffFile $(word 2,$^) -o $@ &> $(LOGDIR)/$@.log")

%.amplicons.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,\
	$(INIT) module load $(TABIX_MODULE) $(BEDTOOLS_MODULE) $(BCFTOOLS_MODULE); \
		$(BEDTOOLS) annotate -counts -i $< -files $(TARGETS_FILE_INTERVALS) | grep -v "#" | \
		awk 'BEGIN {OFS="\t"}{print $$1$(,)$$2-1$(,)$$2$(,)$$NF}'> $<.bed && \
		$(BGZIP) -c $<.bed > $<.bed.gz && $(TABIX) -p bed $<.bed.gz && \
		echo "##INFO=<ID=AMPLICON_NUM$(,)Number=1$(,)Type=Integer$(,)Description=\"Number of the amplicons covering this position\">" > $<.header && \
		$(BCFTOOLS) annotate $< -a $<.bed.gz -c CHROM$(,)FROM$(,)TO$(,)AMPLICON_NUM -h $<.header > $@ && \
		$(RM) $<.bed $<.bed.gz $<.bed.gz.tbi $<.header)

define annotate-sample
vcf/$1.%.ann.vcf : vcf/$1.%.vcf bam/$1.bam bam/$1.bai
	$$(call RUN,4,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(JAVA8_MODULE),"\
	$$(call GATK4241,VariantAnnotator,$$(RESOURCE_REQ_LOW_MEM_JAVA)) -nt 4 -R $$(REF_FASTA) \
	$$(foreach ann,$$(VCF_ANNOTATIONS),-A $$(ann) ) --dbsnp $$(DBSNP) \
	$$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -L $$< -V $$< -o $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach sample,$(SAMPLES),$(eval $(call annotate-sample,$(sample))))

define hrun-sample
vcf/$1.%.hrun.vcf : vcf/$1.%.vcf bam/$1.bam bam/$1.bai
	$$(call RUN,2,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(JAVA8_MODULE),"\
	$$(call GATK4241,VariantAnnotator,$$(RESOURCE_REQ_LOW_MEM_JAVA)) -nt 2 -R $$(REF_FASTA) \
	-L $$< -A HomopolymerRun --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) \
	-V $$< -o $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach sample,$(SAMPLES),$(eval $(call hrun-sample,$(sample))))

%.dbsnp.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(SNP_EFF_MODULE),"\
	$(call SNP_SIFT,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) annotate $(SNP_SIFT_OPTS) $(DBSNP) $< > $@ && $(RM) $^"))

%.cosmic.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(SNP_EFF_MODULE),"\
	$(call SNP_SIFT,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) annotate $(SNP_SIFT_OPTS) $(COSMIC) $< > $@ && $(RM) $^"))

%.clinvar.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(SNP_EFF_MODULE),"\
	$(call SNP_SIFT,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) annotate $(SNP_SIFT_OPTS) $(CLINVAR) $< > $@ && $(RM) $^"))

%.exac_nontcga.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_VSHORT),$(SNP_EFF_MODULE),"\
	$(call SNP_SIFT,$(RESOURCE_REQ_HIGH_MEM_JAVA)) annotate $(SNP_SIFT_OPTS) $(EXAC_NONTCGA) \
	$< > $@ && $(RM) $^"))

%.exac_nonpsych.vcf : %.vcf %.vcf.idx
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_VSHORT),$(SNP_EFF_MODULE),"\
	$(call SNP_SIFT,$(RESOURCE_REQ_HIGH_MEM_JAVA)) annotate $(SNP_SIFT_OPTS) $(EXAC_NONPSYCH) $< > $@ && $(RM) $^"))

%.cadd.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_SHORT),$(SNP_EFF_MODULE),"\
	$(call SNP_SIFT,$(RESOURCE_REQ_HIGH_MEM_JAVA)) annotate $(SNP_SIFT_OPTS) \
	$(if $(findstring indel,$<),$(CADD_INDEL),$(CADD_SNV)) $< > $@ && $(RM) $^"))

#%.gene_ann.vcf : %.vcf
#	$(call CHECK_VCF,$<,$@,\
#	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
#	$(ADD_GENE_LIST_ANNOTATION) --genome $(REF) \
#	--geneBed $(HAPLOTYPE_INSUF_BED)$(,)$(CANCER_GENE_CENSUS_BED)$(,)$(KANDOTH_BED)$(,)$(LAWRENCE_BED)$(,)$(DGD_BED) \
#	--name hap_insuf$(,)cancer_gene_census$(,)kandoth$(,)lawrence$(,)duplicatedGenesDB --outFile $@ $< && $(RM) $< $<.idx"))

%.gene_ann.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,\
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(SNP_EFF_MODULE),"\
	$(call SNP_SIFT,$(RESOURCE_REQ_LOW_MEM_JAVA)) geneSets -v $(GENE_SETS_GMT) $< | \
	sed 's/MSigDb/CancerGeneSets/' > $@ && $(RM) $<")) 

ifdef SAMPLE_PAIRS
define annotate-facets-pair
vcf/$1_$2.%.facets.vcf : vcf/$1_$2.%.vcf facets/cncf/$1_$2.cncf.txt facets/cncf/$1_$2.out 
	$$(call CHECK_VCF,$$<,$$@,\
	purity=`grep Purity $$(<<<) | cut -f2 -d'=' | sed 's/NA/0.1/; s/ //g;'` && \
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
		$$(ANNOTATE_FACETS_VCF) --genome \"$$(REF)\" --tumor \"$1\" \
		--facetsSegTxt $$(<<) --purity $$$$purity --outFile $$@ $$< && \
		$$(RM) $$<"))
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call annotate-facets-pair,$(tumor.$(pair)),$(normal.$(pair)))))
endif

#%.mut_taste.vcf : %.vcf
#	$(INIT) $(call CHECK_VCF,$<,$@,$(MUTATION_TASTER) $< > $@ 2> $(LOG))

#%.chasm.vcf : %.vcf
#	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,8G,00:59:59,"unset PYTHONPATH && source $(CHASM_PYTHON_ENV)/bin/activate $(CHASM_PYTHON_ENV) && \
#		$(CHASM) --genome $(REF) --classifier $(subst $( ),$(,),$(CHASM_CLASSIFIER)) --chasmDir $(CHASM_DIR) --python $(shell which python) --outFile $@ $< && $(RM) $< $<.idx"))

#%.fathmm.vcf : %.vcf
#	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,8G,00:59:59,"PYTHONPATH=$(FATHMM_PYTHONPATH) $(FATHMM) $(FATHMM_OPTS) --outFile $@ $< && $(RM) $< $<.idx"))

#%.mutass.vcf : %.vcf
#	$(call LSCRIPT_MEM,12G,00:59:59,$(MUT_ASS) --outFile $@ --maData $(MUT_ASS_RDATA) $<)

#%.transfic.vcf : %.vcf
#	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_MEM,9G,00:59:59,"$(TRANSFIC) --genome $(REF) --transfic $(TRANSFIC_PERL_SCRIPT) --outFile $@ $< && $(RM) $< $<.idx"))

#%.provean.vcf : %.vcf
#	$(call LSCRIPT_MEM,8G,00:29:29,"$(PROVEAN) $(PROVEAN_OPTS) --outFile $@ $<")

#%.pathogenic.vcf : %.vcf
#	$(INIT) $(call CHECK_VCF,$<,$@,$(CLASSIFY_PATHOGENICITY) $< > $@ 2> $(LOG))

#%.$(ANNOVAR_REF)_multianno.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,7G,9G,"$(ANNOVAR) -out $* $(ANNOVAR_OPTS) $< $(ANNOVAR_DB) && $(RM) $< $<.idx")


define rename-samples-tumor-normal
vcf/$1_$2.%.rn.vcf : vcf/$1_$2.%.vcf
	$$(INIT) module load $(PERL_MODULE); \
	perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' $$< > $$@ && $$(RM) $$< $$<.idx
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call rename-samples-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


#define ad-tumor-normal
#vcf/$1_$2.%.ad.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
#	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,00:29:29,"$$(LOAD_JAVA8_MODULE); $$(call VARIANT_ANNOTATOR,1.5G) -nt 4 -R $$(REF_FASTA) \
#		-A DepthPerAlleleBySample --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@ -L $$<")
#endef
#$(foreach pair,$(SAMPLE_PAIRS),\
#	$(eval $(call ad-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

#define annotate-tumor-normal
#vcf/$1_$2.%.ann.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
#	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,00:29:29,"$$(LOAD_JAVA8_MODULE); $$(call VARIANT_ANNOTATOR,1.5G) -nt 4 -R $$(REF_FASTA) \
#		$$(foreach ann,$$(VCF_ANNOTATIONS),-A $$(ann) ) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@ -L $$< && $$(RM) $$< $$<.idx")
#endef
#$(foreach pair,$(SAMPLE_PAIRS),\
#		$(eval $(call annotate-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

#define hrun-tumor-normal
#vcf/$1_$2.%.hrun.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
#	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,00:29:29,"$$(LOAD_JAVA8_MODULE); $$(call VARIANT_ANNOTATOR,1.5G) -nt 4 -R $$(REF_FASTA) \
#		-A HomopolymerRun --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -L $$< -o $$@ && $$(RM) $$< $$<.idx")
#endef
#$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call hrun-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
#endif

##### extract vcf to table

VCF_FIELDS = CHROM POS ID REF ALT QUAL FILTER
ANN_FIELDS = $(addprefix ANN[*].,ALLELE EFFECT IMPACT GENE GENEID FEATURE FEATUREID BIOTYPE RANK \
	HGVS_C HGVS_P CDNA_POS CDNA_LEN CDS_POS CDS_LEN AA_POS AA_LEN DISTANCE ERRORS)

tables/%.opl_tab.txt : vcf/%.vcf
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(SNP_EFF_MODULE),"\
	format_fields=\$$(grep '^##FORMAT=<ID=' $< | sed 's/dbNSFP_GERP++/dbNSFP_GERP/g' | sed 's/.*ID=//; s/,.*//;' | tr '\n' ' '); \
	N=\$$(expr \$$(grep '^#CHROM' $< | wc -w) - 10); \
	fields='$(VCF_FIELDS)'; \
	for f in \$$format_fields; do \
		for i in \$$(seq 0 \$$N); do \
			fields+=' 'GEN[\$$i].\$$f; \
		done; \
	done; \
	fields+=' '\$$(grep '^##INFO=<ID=' $< | grep -v '=REF,' | sed 's/dbNSFP_GERP++/dbNSFP_GERP/g' | \
		sed 's/.*ID=//; s/,.*//; s/\bANN\b/$(ANN_FIELDS)/; ' | tr '\n' ' '); \
	module load $(PERL_MODULE); $(VCF_EFF_ONE_PER_LINE) < $< | sed 's/dbNSFP_GERP++/dbNSFP_GERP/g' | \
		$(call SNP_SIFT,$(RESOURCE_REQ_LOW_MEM_JAVA)) extractFields - \$$fields > $@; \
	for i in \`seq 0 \$$N\`; do \
	S=\$$(grep '^#CHROM' $< | cut -f \$$((\$$i + 10))); \
	sed -i \"1s/GEN\[\$$i\]/\$$S/g;\" $@; \
	done")

hotspots/%.opl_tab.txt : hotspots/%.vcf
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(SNP_EFF_MODULE),"\
	format_fields=\$$(grep '^##FORMAT=<ID=' $< | sed 's/dbNSFP_GERP++/dbNSFP_GERP/g' | sed 's/.*ID=//; s/,.*//;' | tr '\n' ' '); \
	N=\$$(expr \$$(grep '^#CHROM' $< | wc -w) - 10); \
	fields='$(VCF_FIELDS)'; \
	for f in \$$format_fields; do \
		for i in \$$(seq 0 \$$N); do \
			fields+=' 'GEN[\$$i].\$$f; \
		done; \
	done; \
	fields+=' '\$$(grep '^##INFO=<ID=' $< | grep -v '=REF,' | sed 's/dbNSFP_GERP++/dbNSFP_GERP/g' | \
		sed 's/.*ID=//; s/,.*//; s/\bANN\b/$(ANN_FIELDS)/; ' | tr '\n' ' '); \
	module load $(PERL_MODULE); $(VCF_EFF_ONE_PER_LINE) < $< | sed 's/dbNSFP_GERP++/dbNSFP_GERP/g' | \
		$(call SNP_SIFT,$(RESOURCE_REQ_LOW_MEM_JAVA)) extractFields - \$$fields > $@; \
	for i in \`seq 0 \$$N\`; do \
	S=\$$(grep '^#CHROM' $< | cut -f \$$((\$$i + 10))); \
	sed -i \"1s/GEN\[\$$i\]/\$$S/g;\" $@; \
	done")

sufamscreen/%.opl_tab.txt : sufamscreen/%.vcf
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(SNP_EFF_MODULE),"\
	format_fields=\$$(grep '^##FORMAT=<ID=' $< | sed 's/dbNSFP_GERP++/dbNSFP_GERP/g' | sed 's/.*ID=//; s/,.*//;' | tr '\n' ' '); \
	N=\$$(expr \$$(grep '^#CHROM' $< | wc -w) - 10); \
	fields='$(VCF_FIELDS)'; \
	for f in \$$format_fields; do \
		for i in \$$(seq 0 \$$N); do \
			fields+=' 'GEN[\$$i].\$$f; \
		done; \
	done; \
	fields+=' '\$$(grep '^##INFO=<ID=' $< | grep -v '=REF,' | sed 's/dbNSFP_GERP++/dbNSFP_GERP/g' | \
		sed 's/.*ID=//; s/,.*//; s/\bANN\b/$(ANN_FIELDS)/; ' | tr '\n' ' '); \
	module load $(PERL_MODULE); $(VCF_EFF_ONE_PER_LINE) < $< | sed 's/dbNSFP_GERP++/dbNSFP_GERP/g' | \
		$(call SNP_SIFT,$(RESOURCE_REQ_LOW_MEM_JAVA)) extractFields - \$$fields > $@; \
	for i in \`seq 0 \$$N\`; do \
	S=\$$(grep '^#CHROM' $< | cut -f \$$((\$$i + 10))); \
	sed -i \"1s/GEN\[\$$i\]/\$$S/g;\" $@; \
	done")

%.tab.txt : %.opl_tab.txt
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(PERL_MODULE),"\
	$(VCF_JOIN_EFF) < $< > $@ && $(RM) $<")
	

# merge tables
alltables/all$(PROJECT_PREFIX).%.txt : $(foreach sample,$(SAMPLES),tables/$(sample).%.txt)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(RBIND) --sampleName $< $^ > $@")
ifdef SAMPLE_SETS
alltables/allSS$(PROJECT_PREFIX).%.txt : $(foreach set,$(SAMPLE_SETS),tables/$(set).%.txt)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(RBIND) --normalLast $^ > $@")
endif
ifdef SAMPLE_PAIRS
alltables/allTN$(PROJECT_PREFIX).%.txt : $(foreach pair,$(SAMPLE_PAIRS),tables/$(pair).%.txt)
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(RBIND) --tumorNormal $^ > $@")
endif

%.all_eff.txt : %.txt
	$(INIT) ln -f $< $@

%.synonymous.txt : %.txt
	col_imp=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	col_eff=$$(head -1 $< | tr '\t' '\n' | grep -n "EFFECT" | sed 's/:.*//'); \
	$(INIT) head -1 $< > $@ && awk -v col_imp=$$col_imp -v col_eff=$$col_eff '! (match($$col_imp, /MODERATE/) || match($$col_imp, /HIGH/)) && \
		(match($$col_imp, /LOW/) && (match($$col_eff, /synonymous_variant/)))' $< >> $@

%.nonsynonymous.txt : %.txt
	col_imp=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	$(INIT) head -1 $< > $@ && awk -v col_imp=$$col_imp 'match($$col_imp, /MODERATE/) || match($$col_imp, /HIGH/)' $< >> $@

%.nonsynonymous_synonymous.txt : %.txt
	col_imp=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	col_eff=$$(head -1 $< | tr '\t' '\n' | grep -n "EFFECT" | sed 's/:.*//'); \
	$(INIT) head -1 $< > $@ && awk -v col_imp=$$col_imp -v col_eff=$$col_eff \
	'(match($$col_imp, /MODERATE/) || match($$col_imp, /HIGH/) || \
	(match($$col_imp, /LOW/) && match($$col_eff, /synonymous_variant/)))' $< >> $@

%.nonsynonymous_synonymous_lincRNA.txt : %.txt
	col_imp=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	col_eff=$$(head -1 $< | tr '\t' '\n' | grep -n "EFFECT" | sed 's/:.*//'); \
	col_biotype=$$(head -1 $< | tr '\t' '\n' | grep -n "BIOTYPE" | sed 's/:.*//'); \
	$(INIT) head -1 $< > $@ && awk -v col_imp=$$col_imp -v col_eff=$$col_eff -v col_biotype=$$col_biotype \
	'(match($$col_imp, /MODERATE/) || match($$col_imp, /HIGH/) || \
	(match($$col_imp, /LOW/) && match($$col_eff, /synonymous_variant/)) || \
	(match($$col_biotype, /lincRNA/) && match($$col_eff, /non_coding_exon/))' $< >> $@

%.nonsynonymous_synonymous_lincRNA_upstream.txt : %.txt
	col_imp=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	col_eff=$$(head -1 $< | tr '\t' '\n' | grep -n "EFFECT" | sed 's/:.*//'); \
	col_biotype=$$(head -1 $< | tr '\t' '\n' | grep -n "BIOTYPE" | sed 's/:.*//'); \
	$(INIT) head -1 $< > $@ && awk -v col_imp=$$col_imp -v col_eff=$$col_eff -v col_biotype=$$col_biotype \
	'(match($$col_imp, /MODERATE/) || match($$col_imp, /HIGH/) || \
	(match($$col_imp, /LOW/) && match($$col_eff, /synonymous_variant/)) || \
	(match($$col_biotype, /lincRNA/) && match($$col_eff, /non_coding_exon/)) || \
	(match($$col_eff, /upstream_gene_variant/) && match($$col_imp, /MODIFIER/)))' $< >> $@

#%.nonsynonymous_hotspot.txt : %.txt
#	col_imp=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
#	col_info=$$(head -1 $< | tr '\t' '\n' | grep -n "INFO"); \
#	$(INIT) head -1 $< > $@ && awk -v col_imp=$$col_imp -v col_filter=$$col_info \
#	'(match($$col_imp, /MODERATE/) || match($$col_imp, /HIGH/) || match($$col_info, /HOTSPOT/))' $< >> $@

#%.nonsynonymous_synonymous_hotspot.txt : %.txt
#	col_imp=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
#	col_eff=$$(head -1 $< | tr '\t' '\n' | grep -n "EFFECT" | sed 's/:.*//'); \
#	col_info=$$(head -1 $< | tr '\t' '\n' | grep -n "INFO"); \
#	$(INIT) head -1 $< > $@ && awk -v col_imp=$$col_imp -v col_eff=$$col_eff -v col_info=$$col_info \
#	'(match($$col_imp, /MODERATE/) || match($$col_imp, /HIGH/) || \
#	 (match($$col_imp, /LOW/) && match($$col_eff, /synonymous_variant/)) || \
#	match($$col_info, /HOTSPOT/))' $< >> $@

#%.nonsynonymous_synonymous_hotspot_lincRNA.txt : %.txt
#	col_imp=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
#	col_eff=$$(head -1 $< | tr '\t' '\n' | grep -n "EFFECT" | sed 's/:.*//'); \
#	col_biotype=$$(head -1 $< | tr '\t' '\n' | grep -n "BIOTYPE" | sed 's/:.*//'); \
#	col_info=$$(head -1 $< | tr '\t' '\n' | grep -n "INFO"); \
#	$(INIT) head -1 $< > $@ && awk -v col_imp=$$col_imp -v col_eff=$$col_eff -v col_biotype=$$col_biotype -v col_info=$$col_info \
#	'(match($$col_imp, /MODERATE/) || match($$col_imp, /HIGH/) || (match($$col_imp, /LOW/) && match($$col_eff, /synonymous_variant/)) || \
#	(match($$col_biotype, /lincRNA/) && match($$col_eff, /non_coding_exon/)) || match($$col_info, /HOTSPOT/))' $< >> $@

#%.nonsynonymous_synonymous_hotspot_lincRNA_upstream.txt : %.txt
#	col_imp=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
#	col_eff=$$(head -1 $< | tr '\t' '\n' | grep -n "EFFECT" | sed 's/:.*//'); \
#	col_biotype=$$(head -1 $< | tr '\t' '\n' | grep -n "BIOTYPE" | sed 's/:.*//'); \
#	col_info=$$(head -1 $< | tr '\t' '\n' | grep -n "INFO"); \
#	$(INIT) head -1 $< > $@ && awk -v col_imp=$$col_imp -v col_eff=$$col_eff -v col_biotype=$$col_biotype -v col_filter=$$col_filter \
#	'(match($$col_imp, /MODERATE/) || match($$col_imp, /HIGH/) || (match($$col_imp, /LOW/) && match($$col_eff, /synonymous_variant/)) || \
#	(match($$col_biotype, /lincRNA/) && match($$col_eff, /non_coding_exon/)) || match($$col_info, /HOTSPOT/) || \
#	(match($$col_eff, /upstream_gene_variant/) && match($$col_imp, /MODIFIER/)))' $< >> $@

#### vcf stats

%.vcf.stats : %.vcf
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BCFTOOLS_MODULE),"\
	$(BCFTOOLS) stats $< > $@")

alltables/all$(PROJECT_PREFIX).%.vcf.stats : $(foreach sample,$(SAMPLES),vcf/$(sample).%.vcf.stats)
	$(INIT) \
	echo "SAMPLE" > $@.tmp; \
	egrep "^SN" $< | cut -f 3 >> $@.tmp; \
	for metrics in $^; do \
		echo $${metrics} | sed 's/vcf\///; s/\..*$$//;' > $@.tmp2;\
		egrep "^SN" $${metrics} | cut -f 4 >> $@.tmp2; \
		paste $@.tmp $@.tmp2 > $@;\
		cp $@ $@.tmp;\
	done;
	rm -rf $@.tmp $@.tmp2


#################### MAF ###################
ifdef SAMPLE_PAIRS
define vcf2maf-tumor-normal
maf/$1_$2.%.maf : vcf/$1_$2.%.vcf
	$$(MKDIR) $$(@D)
	$$(call RUN,4,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(VEP_ONCOKB_MODULE)," 	sleep 5 && $$(VCF2MAF) 	--input-vcf $$< 	--tumor-id $1 	--normal-id $2 	--ref-fasta $$(REF_FASTA) 	--vep-path $$(VEP_PATH) 	--vep-data $$(VEP_DATA)  $$(if $$(findstring hg38,$(REF)), --ncbi-build GRCh38) --output-maf $$@ && $$(RM) vcf/$1_$2.*.vep.vcf")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call vcf2maf-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

#define vcf2maf-sample
#maf/$1.%.maf : vcf/$1.%.vcf
#	$$(call LSCRIPT_MEM,9G,12G,"$$(VCF2MAF) --input-vcf $$< --tumor-id $1 --ref-fasta $$(REF_FASTA) --vep-path $$(VEP_PATH) --vep-data $$(VEP_DATA) --output-maf $$@")
#endef
#$(foreach sample,$(SAMPLES),$(eval $(call vcf2maf-sample,$(sample))))

#allmaf/all.%.maf : $(foreach sample,$(SAMPLES),maf/$(sample).%.maf)
#	$(INIT) \
#	{ \
#	sed -n 2p $<; \
#	sed 1,2d $^; \
#	} > $@
#endif


# Copy number regulated genes annotated per subtype
# FYI Endometrioid_MSI-L has no copy number regulated genes
#CN_ENDOMETRIAL_SUBTYPES = CN_high CN_low Endometrioid_MSI_H Endometrioid_MSS Endometrioid MSI POLE Serous
#CN_BREAST_SUBTYPES = ER_negative ER_positive HER2_postitive Pam50_Basal Pam50_Her2 Pam50_LumA Pam50_LumB Pam50_Normal Triple_negative
#CN_ENDOMETRIAL_BED = $(foreach set,$(CN_ENDOMETRIAL_SUBTYPES), $(HOME)/share/reference/annotation_gene_lists/cn_reg/endometrial/copy_number_regulated_genes_subtype_$(set)_spearmanrsquare0.4_fdrbh_adjp_lt0.05.HUGO.bed)
#CN_BREAST_BED = $(foreach set,$(CN_BREAST_SUBTYPES), $(HOME)/share/reference/annotation_gene_lists/cn_reg/breast/metabric_subtype_$(set)_copy_number_regulated_genes_std0.5_spearmanrsquare0.4_fdrbh_adjp_lt0.05.HUGO.bed)
#%.cn_reg.vcf : %.vcf
#	$(call LSCRIPT_MEM,8G,12G,"$(ADD_GENE_LIST_ANNOTATION) --genome $(REF) --geneBed \
#        $(subst $(space),$(,),$(strip $(CN_ENDOMETRIAL_BED)) $(strip $(CN_BREAST_BED))) --name $(subst $(space),$(,),$(foreach set,$(strip $(CN_ENDOMETRIAL_SUBTYPES)),endometrial_#$(set)) $(foreach set,$(strip $(CN_BREAST_SUBTYPES)),breast_$(set))) --outFile $@ $< && $(RM) $< $<.idx")

#-cancer does nothing
#%.som_eff.vcf : %.vcf %.vcf.pair
#	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,14G,"$(call SNP_EFF_MEM,8G) ann -cancer -cancerSamples $(<<) $(SNP_EFF_OPTS) $(SNP_EFF_GENOME) $< > $@ && $(RM) $^"))

# VariantEval: generate vcf report
#reports/%.grp : $(foreach sample,$(SAMPLES),vcf/$(sample).%.vcf) $(foreach sample,$(SAMPLES),vcf/$(sample).%.vcf.idx)
#	$(call LSCRIPT_MEM,2G,5G,"$(call GATK_MEM,2G) -T VariantEval $(foreach sm,$(REPORT_STRATIFICATION), --stratificationModule $(sm)) -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach eval,$(filter %.vcf,$^), --eval:$(call strip-suffix,$(notdir $(eval))) $(eval)) -o $@")
#ifdef SAMPLE_PAIRS
#reports/%.grp : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).%.vcf vcf/$(pair).%.vcf.idx)
#	$(call LSCRIPT_MEM,2G,5G,"$(call GATK_MEM,2G) -T VariantEval $(foreach sm,$(REPORT_STRATIFICATION), --stratificationModule $(sm)) -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach eval,$(filter %.vcf,$^), --eval:$(call strip-suffix,$(notdir $(eval))) $(eval)) -o $@")
#endif

# apply dp filter for somatic sniper
#%.ss_dp_ft.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --filterExpression 'vc.getGenotype(\"TUMOR\").getDP() < $(DEPTH_FILTER) || vc.getGenotype(\"NORMAL\").getDP() < $(DEPTH_FILTER)' --filterName depthFilter && $(RM) $< $<.idx")

# varscan TN variant allele frequency: min tumor freq > 5% ; max normal freq < 5%
#%.freq_ft.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,2G,5G,"sed '/##FORMAT=<ID=FREQ/ s/String/Float/; /^#/! s/%//g' $< | $(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) '(exists GEN[*].FREQ) & (GEN[0].FREQ < 5) & (GEN[1].FREQ[0] > 5)' > $@ && $(RM) $< $<.idx")

#%.vaf_ft.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,2G,5G,"$(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) '(exists GEN[*].VAF) & (GEN[0].VAF > 0.05) & (GEN[1].VAF < 0.05)' < $< > $@ && $(RM) $< $<.idx")


# varscan depth filter (b/c varscan is dumb and only gives variant depth)
#%.vdp_ft.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,2G,5G,"cat $< | $(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) '(exists GEN[*].AD) & (GEN[*].AD > $(DEPTH_FILTER))' > $@ && $(RM) $< $<.idx")

# add exon distance
#%.exondist.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,2G,3G,"$(INTRON_POSN_LOOKUP) $< > $@")

#%.common_ft.vcf : %.vcf
#	$(call LSCRIPT_MEM,4G,5G,"$(COMMON_FILTER_VCF) $< > $@")


#%.fp_ft.vcf : %.vcf
#	$(call LSCRIPT_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration \
#	-R $(REF_FASTA) -V $< -o $@ --maskName 'FuentesFalsePositive' --mask $(FALSE_POSITIVE_BED) && $(RM) $< $<.idx")


#%.hrun_ft.vcf : %.vcf
#	$(call LSCRIPT_CHECK_MEM,8G,01:59:59,"$(LOAD_JAVA8_MODULE); $(VARIANT_FILTRATION) \
#		-R $(REF_FASTA) -V $< -o $@ --filterExpression 'HRun > $(HRUN_FILTER)' --filterName HRun && $(RM) $< $<.idx")

# workaround to do a double pass
#%.pass2.vcf : %.vcf
#	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,2G,00:29:29,"$(LOAD_SNP_EFF_MODULE); $(SNP_SIFT) filter \
#		$(SNP_SIFT_OPTS) -f $< \"( na FILTER ) | (FILTER = 'PASS')\" > $@"))

ifeq ($(findstring tvc,$(MUT_CALLER)),tvc)
include usb-modules-v2/vcf_tools/vcftools_tvc.mk
else
include usb-modules-v2/vcf_tools/vcftools_nontvc.mk
endif

include usb-modules-v2/qc/bamMetrics.mk
VCFTOOLS_MK = true
endif

