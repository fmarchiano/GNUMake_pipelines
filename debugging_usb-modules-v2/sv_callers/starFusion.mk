include usb-modules-v2/Makefile.inc

LOGDIR ?= log/star_fusion.$(NOW)

PHONY += star_fusion
.PHONY : $(PHONY)
.DELETE_ON_ERROR:

star_fusion : star_fusion/all.star-fusion.STAR-Fusion.final

define star-fusion
star_fusion/$1/star-fusion.STAR-Fusion.filter.ok : star/$1.Chimeric.out.junction.gz
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(STAR_FUSION_MODULE),"\
	gunzip -c $$< | \
	awk '\$$$$1!~\"GL0\" && \$$$$1!~\"GRCm\" && \$$$$1!~\"hs\" && \$$$$4!~\"GL0\" && \$$$$4!~\"GRCm\" && \$$$$4!~\"hs\"' | \
	awk 'BEGIN {OFS = \"\t\"} {print \"chr\"$$(,) \$$$$1$$(,) \$$$$2$$(,) \$$$$3$$(,) \"chr\"$$(,) \$$$$4$$(,) \$$$$5$$(,) \
	\$$$$6$$(,) \$$$$7$$(,) \$$$$8$$(,) \$$$$9$$(,) \$$$$10$$(,) \$$$$11$$(,) \$$$$12$$(,) \$$$$13$$(,) \$$$$14}' | \
	sed 's/chr\t/chr/g' | sed 's/chrMT/chrM/g' > $$<.tmp && \
	$$(STAR_FUSION) --genome_lib_dir $$(STAR_CTAT_DIR) \
	-J $$<.tmp --output_dir star_fusion/$1 \
	--min_junction_reads $2 --min_sum_frags $3 --max_promiscuity $4 \
	--min_novel_junction_support $5 --min_alt_pct_junction $6 || exit 0; \
	$$(RM) $$<.tmp; $$(RMR) $$(@D)/*intermediates_dir; ")

star_fusion/$1/star-fusion.fusion_candidates.final : star_fusion/$1/star-fusion.STAR-Fusion.filter.ok
	

endef

ifdef NORMAL_SAMPLES
$(foreach sample,$(NORMAL_SAMPLES),$(eval $(call star-fusion,$(sample),$(STAR_FUSION_MIN_JUNCTION_READS_NORMAL),\
	$(STAR_FUSION_MIN_SUM_FRAGS_NORMAL),$(STAR_FUSION_MAX_PROMISCUITY_NORMAL),\
	$(STAR_FUSION_MIN_NOVEL_JUNCTION_SUPPORT_NORMAL),$(STAR_FUSION_MIN_ALT_PCT_JUNC_NORMAL))))
$(foreach sample,$(TUMOR_SAMPLES),$(eval $(call star-fusion,$(sample),$(STAR_FUSION_MIN_JUNCTION_READS),\
	$(STAR_FUSION_MIN_SUM_FRAGS),$(STAR_FUSION_MAX_PROMISCUITY),\
	$(STAR_FUSION_MIN_NOVEL_JUNCTION_SUPPORT),$(STAR_FUSION_MIN_ALT_PCT_JUNC))))
else
$(foreach sample,$(SAMPLES),$(eval $(call star-fusion,$(sample),$(STAR_FUSION_MIN_JUNCTION_READS),\
	$(STAR_FUSION_MIN_SUM_FRAGS),$(STAR_FUSION_MAX_PROMISCUITY),\
	$(STAR_FUSION_MIN_NOVEL_JUNCTION_SUPPORT),$(STAR_FUSION_MIN_ALT_PCT_JUNC))))
endif

ifdef NORMAL_SAMPLES
star_fusion/all.star-fusion.STAR-Fusion.final : $(foreach sample,$(TUMOR_SAMPLES) $(NORMAL_SAMPLES),star_fusion/$(sample)/star-fusion.fusion_candidates.final)
	$(INIT) \
	{ \
	sed "s/^#fusion/SAMPLE\tfusion/" $< | head -1; \
	for fusion in $^; do \
		samplename=`dirname $${fusion}`; \
		samplename=`basename $${samplename}`; \
		sed "/^#/d; s/^/$$samplename\t/" $$fusion; \
	done; \
	} > $@
else
star_fusion/all.star-fusion.STAR-Fusion.final : $(foreach sample,$(SAMPLES),star_fusion/$(sample)/star-fusion.fusion_candidates.final)
	$(INIT) \
	{ \
	sed "s/^#fusion/SAMPLE\tfusion/" $< | head -1; \
	for fusion in $^; do \
		samplename=`dirname $${fusion}`; \
		samplename=`basename $${samplename}`; \
		sed "/^#/d; s/^/$$samplename\t/" $$fusion; \
	done; \
	} > $@
endif
include usb-modules-v2/aligners/starAligner.mk
