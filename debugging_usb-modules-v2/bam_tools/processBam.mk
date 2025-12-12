# various bam processing steps

# can be used to reprocess bam files, merge them, or merge and reprocess bam files
# possible post-processing steps are defined in modules/aligners/align.inc
##### MAKE INCLUDES #####

ifndef PROCESS_BAM_MK

include usb-modules-v2/Makefile.inc

LOGDIR ?= log/process_bam.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY: 

BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)

ifeq ($(BAM_REPROCESS),true)
processed_bams : $(BAMS) $(addsuffix .bai,$(BAMS))
bam/%.bam : unprocessed_bam/%.$(BAM_SUFFIX)
	$(INIT) ln -f $< $@
else
ifeq ($(MERGE_SPLIT_BAMS),true)
merged_bams : $(BAMS) $(addsuffix .bai,$(BAMS))
bam/%.bam : unprocessed_bam/%$(if $(findstring true,$(BAM_FIX_RG)),.rg).bam
	$(INIT) ln -f $< $@
endif
endif
$(info CHROMOSOMES $(CHROMOSOMES))
ifeq ($(MERGE_SPLIT_BAMS),true)
define bam-header
unprocessed_bam/$1.header.sam : $$(foreach split,$2,unprocessed_bam/$$(split).bam)
	$$(INIT) module load $$(SAMTOOLS_MODULE); $$(SAMTOOLS) view -H $$< | grep -v '^@RG' > $$@.tmp; \
	for bam in $$(^M); do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $$(RM) $$@.tmp
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call bam-header,$(sample),$(split.$(sample)))))

define merged-bam
unprocessed_bam/$1.bam : unprocessed_bam/$1.header.sam $$(foreach split,$2,unprocessed_bam/$$(split).bam)
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(SAMTOOLS_MODULE),"\
	$$(SAMTOOLS) merge -f -h $$< $$@ $$(filter %.bam,$$^)")
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-bam,$(sample),$(split.$(sample)))))
endif

# indices
# if bam file is a symlink, need to create a symlink to index
index : $(addsuffix .bai,$(BAMS))

%.bam.bai : %.bam
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(SAMTOOLS_MODULE),"\
	sleep 5 && $(SAMTOOLS) index $< && ln -f $@ $*.bai")

%.bai : %.bam.bai
	$(INIT) sleep 5; ln -f $@ $<

%.bam : %.sam
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(SAMTOOLS_MODULE),"\
	$(SAMTOOLS) view -bSh $< > $@" && $(RM) $<)

# limit coverage
#%.dcov.bam : %.bam
#	$(call LSCRIPT_MEM,18G,00:59:59,"$(call GATK_MEM,18G) -T PrintReads -R $(REF_FASTA) -I $< -dcov 50 -o $@")

%.downsampled.bam : %.bam
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(SAMTOOLS_MODULE),"\
	$(SAMTOOLS) view -bh -s $(SAMTOOLS_DOWNSAMPLE_FACTOR) $< > $@")

%.fixmate.bam : %.bam
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call PICARD,FixMateInformation,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) I=$< O=$@ && $(RM) $<")

# needs LENIENT when working on bams from bwa aln ("MAPQ should be 0 for unmapped read" error)
%.reordered.bam : %.bam $(REF_DICT)
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call PICARD,ReorderSam,$(RESOURCE_REQ_HIGH_MEM_JAVA)) VALIDATION_STRINGENCY=LENIENT I=$< O=$@ SEQUENCE_DICTIONARY=$(REF_DICT) && $(RM) $<")

# needs LENIENT when working on bams from bwa aln ("MAPQ should be 0 for unmapped read" error)
%.sorted.bam : %.bam
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_LONG),$(JAVA8_MODULE),"\
	$(call PICARD,SortSam,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) VALIDATION_STRINGENCY=LENIENT I=$< O=$@ SO=coordinate VERBOSITY=ERROR && $(RM) $^")

%.nsorted.bam : %.bam
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_LONG),$(JAVA8_MODULE),"\
	$(call PICARD,SortSam,$(RESOURCE_REQ_HIGH_MEM_JAVA)) I=$< O=$@ SO=queryname")

%.clean.bam : %.bam
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call PICARD,CleanSam,$(RESOURCE_REQ_LOW_MEM_JAVA)) I=$< O=$@")

%.filtered.bam : %.bam %.bai
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_MEDIUM),$(SAMTOOLS_MODULE),"\
	$(SAMTOOLS) view -bF $(BAM_FILTER_FLAGS) $< > $@ && $(RM) $^")

%.markdup.bam : %.bam %.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_VVHIGH_MEM),$(RESOURCE_REQ_LONG),$(JAVA8_MODULE),"\
	$(MKDIR) metrics; \
	$(call PICARD,MarkDuplicates,$(RESOURCE_REQ_VVHIGH_MEM_JAVA)) I=$< O=$@ \
	METRICS_FILE=metrics/$(call strip-suffix,$(@F)).dup_metrics.txt && $(RM) $^")

%.rmdup.bam : %.bam
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(SAMTOOLS_MODULE),"\
	$(SAMTOOLS) rmdup $< $@ && $(RM) $^")

%.splitntrim.bam : %.bam %.bai
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_MEDIUM),$(JAVA8_MODULE),"\
	$(call GATK,SplitNCigarReads,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -I $< -o $@ \
	-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -R $(REF_FASTA)")

%.intrachr.bam : %.bam %.bai
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),$(SAMTOOLS_MODULE),"\
	$(SAMTOOLS) view -h $< | awk '\$$1~\"@\" || \$$7==\"=\"' | $(SAMTOOLS) view -Sb >> $@ && $(RM) $^")


ifneq ($(KNOWN_INDELS),)
BAM_REALN_OPTS = --knownAlleles $(KNOWN_INDELS)
BAM_REALN_TARGET_OPTS = --known $(KNOWN_INDELS)
endif

%.intervals : %.bam %.bam.bai
	$(call RUN,4,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(JAVA8_MODULE),"\
	$(call GATK,RealignerTargetCreator,$(RESOURCE_REQ_LOW_MEM_JAVA)) \
	-I $< -nt 4 -R $(REF_FASTA) -o $@ $(BAM_REALN_TARGET_OPTS)")

%.realn.bam : %.bam %.intervals %.bai %.bam.bai
	if [[ -s $(word 2,$^) ]]; then \
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call GATK,IndelRealigner,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-I $< -R $(REF_FASTA) -targetIntervals $(<<) -o $@ $(BAM_REALN_OPTS) && $(RM) $^") ; \
	else mv $< $@ ; fi

ifeq ($(findstring true,$(CHOOSE_CHR_FOR_RECAL)),true)
%.$(BAM_SUFFIX1).$(RECAL_CHR).$(BAM_SUFFIX2)_report.grp : %.$(BAM_SUFFIX1).$(RECAL_CHR).$(subst .recal,$(),$(BAM_SUFFIX2)).bam %.$(BAM_SUFFIX1).$(RECAL_CHR).$(subst .recal,$(),$(BAM_SUFFIX2)).bam.bai
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_LONG),$(JAVA8_MODULE),"\
	$(call GATK,BaseRecalibrator,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $(REF_FASTA) -knownSites $(DBSNP) $(if $(TARGETS_FILE_INTERVALS),-L $(TARGETS_FILE_INTERVALS)) \
	-L $(RECAL_CHR) -I $< -o $@ \
	$(foreach chr,$(filter-out $(RECAL_CHR),$(CHROMOSOMES)), && ln -f $@ $*.$(BAM_SUFFIX1).$(chr).$(BAM_SUFFIX2)_report.grp)")

define chr-recal-report
%.$(BAM_SUFFIX1).$1.$(BAM_SUFFIX2)_report.grp : %.$(BAM_SUFFIX1).$(RECAL_CHR).$(BAM_SUFFIX2)_report.grp
	$$(INIT) ln -f $$< $$@
endef
$(foreach chr,$(filter-out $(RECAL_CHR),$(CHROMOSOMES)),$(eval $(call chr-recal-report,$(chr))))
else
%.recal_report.grp : %.bam %.bai
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_LONG),$(JAVA8_MODULE),"\
	$(call GATK,BaseRecalibrator,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) \
	-R $(REF_FASTA) -knownSites $(DBSNP) \
	$(if $(TARGETS_FILE_INTERVALS),-L $(TARGETS_FILE_INTERVALS)) -I $< -o $@")
endif

%.recal.bam : %.bam %.recal_report.grp %.bai %.bam.bai
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call GATK,PrintReads,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -R $(REF_FASTA) -I $< \
	-BQSR $(word 2,$^) -o $@ && $(RM) $^")


ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
%.rg.bam : %.bam
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call PICARD,AddOrReplaceReadGroups,$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) I=$< O=$@ \
	RGLB=$(call strip-suffix,$(@F)) RGPL=$(SEQ_PLATFORM) RGPU=00000000 \
	RGSM=$(call strip-suffix,$(@F)) RGID=$(call strip-suffix,$(@F)) \
	VERBOSITY=ERROR && $(RM) $<")
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
%.rg.bam : %.bam
	$(INIT) module load $(SAMTOOLS_MODULE); \
		samplename=`basename $< .bam` && \
		$(SAMTOOLS) view -H $< | sed "s/SM:[a-zA-Z0-9 _\.-]*/SM:$${samplename}/" > $<.header
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(SAMTOOLS_MODULE),"\
	$(SAMTOOLS) reheader $<.header $< > $@")
endif

ifeq ($(SPLIT_CHR),true)
define chr-splitchr
%.$1.splitchr.bam : %.bam %.bam.bai
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_MEDIUM),$$(JAVA8_MODULE),"\
	$$(call GATK,PrintReads,$$(RESOURCE_REQ_MEDIUM_MEM_JAVA)) -L $1 -R $$(REF_FASTA) -I $$< -o $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-splitchr,$(chr))))

%.$(BAM_SUFFIX2).merged.bam : $(foreach chr,$(CHROMOSOMES),%.$(chr).$(BAM_SUFFIX2).bam) $(foreach chr,$(CHROMOSOMES),%.$(chr).$(BAM_SUFFIX2).bam.bai) $(foreach chr,$(CHROMOSOMES),%.$(chr).$(BAM_SUFFIX2).bai)
	$(call RUN,1,$(RESOURCE_REQ_VHIGH_MEM),$(RESOURCE_REQ_VLONG),$(JAVA8_MODULE),"\
	$(call PICARD,MergeSamFiles,$(RESOURCE_REQ_VHIGH_MEM_JAVA)) \
	$(foreach i,$(filter %.bam,$^), I=$(i)) SORT_ORDER=coordinate O=$@ USE_THREADING=true && $(RM) $^")

endif # end chrom splitting



endif
PROCESS_BAM_MK = true


