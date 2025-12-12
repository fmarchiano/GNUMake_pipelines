ifndef ALIGN_INC

BAM_SUFFIX1 = $(subst $(space),.,$(strip \
	$(if $(findstring true,$(BAM_REORDER)),reordered)\
	$(if $(findstring true,$(BAM_SORT)),sorted)\
	$(if $(findstring true,$(BAM_FIX_RG)),rg)\
	$(if $(findstring true,$(BAM_FILTER)),filtered)\
	$(if $(findstring true,$(BAM_INTRACHR)),intrachr)\
	$(if $(findstring true,$(BAM_SPLIT_N_TRIM)),splitntrim)))

BAM_SUFFIX2 = $(subst $(space),.,$(strip \
	$(if $(findstring true,$(SPLIT_CHR)),splitchr)\
	$(if $(findstring rmdup,$(BAM_DUP_TYPE)),rmdup)\
	$(if $(findstring markdup,$(BAM_DUP_TYPE)),markdup)\
	$(if $(findstring true,$(BAM_REALN)),realn)\
	$(if $(findstring true,$(BAM_RECAL)),recal)))
	
BAM_SUFFIX = $(subst $(space),.,$(strip \
	$(BAM_SUFFIX1) \
	$(BAM_SUFFIX2) \
	$(if $(findstring true,$(SPLIT_CHR)),merged)\
	bam))
$(info BAM_SUFFIX-ALIGN is $(BAM_SUFFIX))

endif
ALIGN_INC = true

