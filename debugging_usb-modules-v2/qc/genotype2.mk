##### DEFAULTS ######
LOGDIR ?= log/genotype2.$(NOW)

##### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

ifeq ($(shell test $(words $(SAMPLES)) -gt 2; echo $$?),0)
all : genotype/BAMixChecker/Total_result.txt
endif

ifneq ($(findstring RNA,$(CAPTURE_METHOD)),RNA)
genotype/BAMixChecker/Total_result.txt : $(foreach sample,$(SAMPLES),bam/$(sample).bam)
	$(call RUN,6,$(RESOURCE_REQ_VHIGH_MEM),$(RESOURCE_REQ_MEDIUM),$(SINGULARITY_MODULE),"\
	echo $^ | tr ' ' '\n' > genotype/BAMixChecker/bamlist &&\
	$(SINGULARITY_EXEC) $(BAMIXCHECKER_IMG) python /BAMixChecker-1.0.1/BAMixChecker.py -l genotype/BAMixChecker/bamlist \
	$(if $(findstring BAITS, $(CAPTURE_METHOD)), -b $(TARGETS_FILE_INTERVALS), $(if $(findstring PCR,$(CAPTURE_METHOD)), -b $(TARGETS_FILE_INTERVALS))) \
	-r $(REF_FASTA) \
	-o genotype --OFFFileNameMatching -p 6")
endif

ifeq ($(findstring RNA,$(CAPTURE_METHOD)),RNA)
bam_splitncigar/%.splitncigar.bam : bam/%.bam
	$(call RUN,1,$(RESOURCE_REQ_VHIGH_MEM),$(RESOURCE_REQ_LONG),$(JAVA8_MODULE),"\
	$(call GATK4241,SplitNCigarReads,$(RESOURCE_REQ_VHIGH_MEM_JAVA)) \
	--tmp-dir $(TMPDIR) \
	-R $(REF_FASTA) -I $< -O $@")

genotype/BAMixChecker/Total_result.txt : $(foreach sample,$(SAMPLES),bam_splitncigar/$(sample).splitncigar.bam)
	$(call RUN,6,$(RESOURCE_REQ_VHIGH_MEM),$(RESOURCE_REQ_MEDIUM),$(SINGULARITY_MODULE),"\
	echo $^ | tr ' ' '\n' > genotype/BAMixChecker/bamlist &&\
	$(SINGULARITY_EXEC) $(BAMIXCHECKER_IMG) python /BAMixChecker-1.0.1/BAMixChecker.py -l genotype/BAMixChecker/bamlist \
	-r $(REF_FASTA) \
	-o genotype --OFFFileNameMatching -p 6")
endif
