HPC ?= scicore
include usb-modules-v2/env/$(HPC)/Makefile.inc


# GLOBAL CONF:


############################################
######## DEFINE SAMPLES/SETS ###############
############################################


ifndef MAKEFILE_INC

PROJECT_NAME = $(shell pwd | sed 's:.*/projects/::; s:.*/data/::; s:.*kinglab/::; s:/:_:g')

ifneq ($(wildcard $(SAMPLE_FILE)),)
  SAMPLES ?= $(shell sed '/^\#/d' $(SAMPLE_FILE))
endif
ifneq ($(wildcard $(SAMPLE_PON_FILE)),)
  PANEL_OF_NORMAL_SAMPLES ?= $(shell sed '/^\#/d' $(SAMPLE_PON_FILE))
endif
ifneq ($(wildcard $(SAMPLE_POOLEDNORM_FILE)),)
  POOLED_NORM_SAMPLES ?= $(shell sed '/^\#/d' $(SAMPLE_POOLEDNORM_FILE))
endif
ifneq ($(wildcard $(MERGE_SAMPLE_FILE)),)
  MERGE_SAMPLES ?= $(shell sed '/^\#/d' $(MERGE_SAMPLE_FILE))
endif


get_tumors = $(wordlist 1,$(shell expr $(words $1) - 1),$1)
get_normal = $(lastword $1)
get_space_sets = $(shell sed '/^\#/d; s/\s\+/ /g; s/\s\+$$//;' $(SAMPLE_SET_FILE) | sed -n '$(1)p')
get_underscore_sets = $(shell sed '/^\#/d; s/\s\+/_/g; s/\s\+$$//;' $(SAMPLE_SET_FILE) | sed -n '$(1)p')

ifneq ($(wildcard $(SAMPLE_SET_FILE)),)
  NUM_SETS := $(shell sed '/^\#/d' $(SAMPLE_SET_FILE) | wc -l | cut -f 1 -d' ')
  SETS_SEQ := $(shell seq 1 $(NUM_SETS))
  $(foreach i,$(SETS_SEQ), \
	  $(eval set.$i := $(call get_space_sets,$i)))
  # define set_lookup.$(sample)
  $(foreach i,$(SETS_SEQ), \
	  $(foreach sample,$(set.$i), \
		  $(eval set_lookup.$(sample) := $i)))
  # define SAMPLE_SETS to contain sample_sample_normal
  $(foreach i,$(SETS_SEQ), \
	  $(eval SAMPLE_SETS += $(call get_underscore_sets,$i)))
  # tumor.SET => tumors
  $(foreach i,$(SETS_SEQ), \
	$(eval tumor.$(call get_underscore_sets,$i) := $(call get_tumors,$(set.$i))))
  # normal.SET => normal
  $(foreach i,$(SETS_SEQ), \
	$(eval normal.$(call get_underscore_sets,$i) := $(call get_normal,$(set.$i))))
  NORMAL_SAMPLES = $(filter-out poolednorm,$(foreach i,$(SETS_SEQ),$(call get_normal,$(set.$i))))
  TUMOR_SAMPLES = $(foreach i,$(SETS_SEQ),$(call get_tumors,$(set.$i)))
  SAMPLE_PAIRS = $(foreach set,$(SAMPLE_SETS),$(foreach tumor,$(tumor.$(set)),$(tumor)_$(normal.$(set))))
  # define tumor.pair and normal.pair to retrieve tumor/normal from pairs
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor.$(tumor)_$(normal.$(set)) := $(tumor))))
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval normal.$(tumor)_$(normal.$(set)) := $(normal.$(set)))))
  # tumor_normal.$(tumor) to retrieve tumor_normal
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor_normal.$(tumor) := $(tumor)_$(normal.$(set)))))
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor_normal.$(normal.$(set)) := $(tumor)_$(normal.$(set)))))
  # set.$(sample) to retrieve sample sets (underscores)
  $(foreach i,$(SETS_SEQ), \
	  $(foreach sample,$(set.$i), \
		  $(eval set.$(sample) := $(word $i,$(SAMPLE_SETS)))))
  UNMATCHED_SAMPLES = $(shell sed '/^\#/d' $(SAMPLE_FILE) $(SAMPLE_SET_FILE) | tr ' ' '\n' | sort | uniq -u)

  SAMPLE_SET_PAIRS = $(shell echo "$(SAMPLE_PAIRS) $(SAMPLE_SETS)" | tr ' ' '\n' | sort | uniq)
  $(foreach set,$(SAMPLE_SET_PAIRS), \
	$(eval samples.$(set) := $(shell echo "$(set)" | sed 's/_/ /g')))
endif

ifneq ($(wildcard $(SAMPLE_SPLIT_FILE)),)
  A = $(shell sed '/^\#/d' $(SAMPLE_SPLIT_FILE) | cut -f1)
  B = $(shell sed '/^\#/d' $(SAMPLE_SPLIT_FILE) | cut -f2)
  $(foreach i,$(shell seq 1 $(words $(A))),$(eval split.$(word $i,$(A)) += $(word $i,$(B))))
  UNSPLIT_SAMPLES = $(B)
  SPLIT_SAMPLES = $(shell sed '/^\#/d' $(SAMPLE_SPLIT_FILE) | cut -f1 | sort | uniq)
endif

ALL_SAMPLES = $(SAMPLE_PAIRS) $(SAMPLES)
COUNT_SAMPLES = $(shell expr `sed 's/\..*//; s:.*/::' <<< $1 | grep -o "_" | wc -l` + 1)
PANEL_OF_NORMAL_SAMPLES ?= $(NORMAL_SAMPLES)
POOLED_NORM_SAMPLES ?= $(NORMAL_SAMPLES)


#######################################################
########### DEFINE GENOME AND RESOURCES ###############
#######################################################

RESOURCE_REQ ?= medium
include usb-modules-v2/genome_inc/$(REF).inc
include usb-modules-v2/resources.inc

#####################################
### DEFINE BASIC FUNCTIONS ##########
#####################################

NO_RM ?= false
ifeq ($(NO_RM),true)
    RM := touch
    RMR = touch
else
    RM := rm -f
    RMR := rm -fr
endif

# define $(,) and $( ) for escaping commas and spaces
, := ,
# $( ) variable not working in make v4.3 (empty variable name error for '$(space) :=':
# space :=
# space +=
# $(space) :=
# $(space) +=
E :=
space := $E $E


# $(call strip-suffix,filename)
strip-suffix = $(firstword $(subst ., ,$(1)))

LINK = ln -svf $(notdir $1) $(2) && touch $1
NOW := $(shell date +"%F")
MKDIR = mkdir -p -m 775
MKDIRS = $(MKDIR) $(LOGDIR)/$(@D) $(@D)
LOG = $(LOGDIR)/$(@).log
UMASK = 002

INIT = module purge; umask $(UMASK); $(MKDIRS); set -o pipefail;
PURGE_AND_LOAD = module purge && module load $1
#PURGE_AND_LOAD = export LMOD_DISABLE_SAME_NAME_AUTOSWAP='no' && module purge && module load $1

MD5 = md5sum $(@:.md5=) > $@
CHECK_MD5 = for i in {0..20}; do if md5sum -c $(filter %.md5,$^); then break; fi; sleep 5; done;
<M = $(<:.md5=)
^M = $(^:.md5=)
@M = $(@:.md5=)
<<M = $(word 2,$(^M))
<<<M = $(word 3,$(^M))

# $(call CHECK_VCF,vcf.file,target,command)
# check for valid header, ie non-empty vcf file
CHECK_VCF = if [ `grep -v '^\#' $1 | wc -l` -eq 0 ] && [ `grep '^\#CHROM' $1 | wc -l` -eq 1 ]; then cp $1 $2; else $3; fi
CHECK_MAF = if [ `grep -v '^\#' $1 | wc -l` -eq 1 ] && [ `grep '^Hugo_Symbol' $1 | wc -l` -eq 1 ]; then cp $1 $2; else $3; fi
CHECK_VCF_CMD = if [ `grep -v '^\#' $1 | wc -l` -eq 0 ] && [ `grep '^\#CHROM' $1 | wc -l` -eq 1 ]; then $2; else $3; fi
CHECK_VCF_HIGH_MODERATE_CMD = if [ `grep -v '^\#' $1 | egrep "HIGH|MODERATE" | wc -l` -eq 0 ] && [ `grep '^\#CHROM' $1 | wc -l` -eq 1 ]; then $2; else $3; fi

<< = $(word 2,$^)
<<< = $(word 3,$^)
<<<< = $(word 4,$^)
<<<<< = $(word 5,$^)
<<<<<< = $(word 6,$^)

define R_INIT
dir.create('$(@D)', showWarnings = F, recursive = T)
qw <- function(x) unlist(strsplit(x, "[[:space:]]+"))
endef

endif
MAKEFILE_INC = true
