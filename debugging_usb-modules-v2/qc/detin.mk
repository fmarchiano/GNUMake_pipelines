
##### DEFAULTS ######
LOGDIR ?= log/detin.$(NOW)

##### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc


.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : detin

detin : $(foreach pair,$(SAMPLE_PAIRS), detin/$(tumor.$(pair))_$(normal.$(pair))/$(tumor.$(pair)).TiN_estimate.txt) 

define detin-tumor-normal
detin/$1_$2/$1_$2.cncf4detin.txt : facets/cncf/$1_$2.cncf.txt
	$$(INIT) awk '{FS="\t";OFS="\t"}{print $$$$1,$$$$10,$$$$11,$$$$3,$$$$11-$$$$10,$$$$3,$$$$6,2^($$$$5+1)}' $$^ | sed 's/^chrom.*/Chromosome\tStart.bp\tEnd.bp\tn_probes\tlength\tn_hets\tf\ttau/' > $$@

detin/$1_$2/$1.hets4detin.tsv : gatk/dbsnp/$1.gatk_snps.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE) $$(SNP_EFF_MODULE),"\
	$$(BCFTOOLS) view -g het $$^ | $$(call SNP_SIFT,$$(RESOURCE_REQ_LOW_MEM_JAVA)) extractFields - CHROM POS GEN[$1].AD[0] GEN[$1].AD[1] REF ALT | sed -e 's/^CHROM.*/CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE/g' -e 's/^chr//' > $$@")

detin/$1_$2/$2.hets4detin.tsv : gatk/dbsnp/$2.gatk_snps.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE) $$(SNP_EFF_MODULE),"\
	$$(BCFTOOLS) view -g het $$^ | $$(call SNP_SIFT,$$(RESOURCE_REQ_LOW_MEM_JAVA)) extractFields - CHROM POS GEN[$2].AD[0] GEN[$2].AD[1] REF ALT | sed -e 's/^CHROM.*/CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE/g' -e 's/^chr//' > $$@")

detin/$1_$2/$1.muts4detin.tsv : $(wildcard vcf/$1_$2.mutect2.*.hotspot.vcf)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE) $$(SNP_EFF_MODULE),"\
	$$(BCFTOOLS) view --max-alleles 2 $$^ | $$(call SNP_SIFT,$$(RESOURCE_REQ_LOW_MEM_JAVA)) extractFields - CHROM POS REF ALT GEN[$1].AD[0] GEN[$1].AD[1] GEN[$2].AD[0] GEN[$2].AD[1] FILTER | sed -E -e '/\tPASS/! s/$$$$/\tREJECT/' -e 's/\tPASS/\t\tKEEP/' -e 's/^chr//' -e 's/^/$1\t$2\t/' -e 's/.+\tCHROM\tPOS.+/tumor_name\tnormal_name\tcontig\tposition\tref_allele\talt_allele\tt_ref_count\tt_alt_count\tn_ref_count\tn_alt_count\tfailure_reasons\tjudgement/g' -e 's/normal_artifact/alt_allele_in_normal/' > $$@")

detin/$1_$2/$1.TiN_estimate.txt : detin/$1_$2/$1_$2.cncf4detin.txt detin/$1_$2/$1.hets4detin.tsv detin/$1_$2/$2.hets4detin.tsv detin/$1_$2/$1.muts4detin.tsv
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT), $$(MATPLOTLIB_MODULE),"\
	$$(DETIN) --mutation_data_path detin/$1_$2/$1.muts4detin.tsv --cn_data_path detin/$1_$2/$1_$2.cncf4detin.txt \
	--tumor_het_data detin/$1_$2/$1.hets4detin.tsv --normal_het_data detin/$1_$2/$2.hets4detin.tsv \
	--exac_data_path $$(DETIN_PICKLE) \
	--output_name $1 --output_dir detin/$1_$2")

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call detin-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


