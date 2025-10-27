# work_standard_makefile

# change this as necessary; see utils/makefile_header for definitions.

include $(MESA_DIR)/utils/makefile_header

LOAD_MESA = $(LOAD_EXTRAS1) $(LOAD_MESA_STAR) $(LOAD_EXTRAS)

ifndef STAR
STAR = star
endif

# Set path
WORK_DIR = ..
WORK_SRC_DIR = $(WORK_DIR)/src
STAR_JOB_DIR = $(MESA_DIR)/star/job

# Set object files
RUN_STAR_EXTRA_OBJS  =
RUN_STAR_EXTRA_OBJS += custom_eos.o
RUN_STAR_EXTRA_OBJS += custom_kap.o

OBJS_STAR = $(RUN_STAR_EXTRA_OBJS) $(STAR_WORK_OBJS) run.o
OBJS_PLANET = $(RUN_STAR_EXTRA_OBJS) run_star_extras.o run_star.o run.o

LOAD_MESA = $(LOAD_MESA_STAR)

WORK_COMPILE = \
   $(LOADER) $(FCbasic) $(FCopenmp) $(FCchecks) $(FCopt) $(FCdebug) $(FCfree) $(FCstandard) \
   -I$(MESA_INCLUDE_DIR) -I$(WORK_SRC_DIR) -I$(MESA_DIR)/star/make -c

# Target files

STAR = star
PLANET = planet

TARGET =
ifdef TAG
	TARGET+=$(TAG)
else
	TARGET+=$(STAR)
	TARGET+=$(PLANET)
endif

all: $(TARGET)

# Compile targets

$(STAR) : $(OBJS_STAR)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(STAR) $(OBJS_STAR) $(LOAD_MESA)

$(PLANET) : $(OBJS_PLANET)
	$(LOADER) $(FCopenmp) -o $(WORK_DIR)/$(PLANET) $(OBJS_PLANET) $(LOAD_MESA)

# Compile object files
%.o: $(WORK_SRC_DIR)/%.f90
	$(WORK_COMPILE) $<

%.o: $(STAR_JOB_DIR)/%.f90
	$(WORK_COMPILE) $<

clean:
	-@rm -f *.o *.mod *.smod *.so *.a $(WORK_DIR)/$(STAR)
	-@rm -f *.o *.mod *.smod *.so *.a $(WORK_DIR)/$(PLANET)

remk:
	-@rm -f run.o $(WORK_DIR)/$(STAR)
	-@rm -f run.o $(WORK_DIR)/$(PLANET)
