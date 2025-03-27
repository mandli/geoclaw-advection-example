
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = geoclaw                  # Clawpack package to use
EXE = xgeoclaw                 # Executable to create
SETRUN_FILE = setrun.py        # File containing function to make data
OUTDIR = _output               # Directory for output
SETPLOT_FILE = setplot.py      # File containing function to set plots
PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran
# NETCDF4_FORT_PREFIX = /opt/homebrew/Cellar/netcdf-fortran/4.6.1
# FFLAGS += -DNETCDF -I$(NETCDF4_FORT_PREFIX)/include -L$(NETCDF4_FORT_PREFIX)/LIB
# FFLAGS += -DNETCDF $(shell nf-config --fflags)
# LFLAGS += $(FFLAGS) $(shell nf-config --flibs) $(shell nc-config --libs)

# ---------------------------------
# package sources for this program:
# ---------------------------------

GEOLIB = $(CLAW)/geoclaw/src/2d/shallow
include $(GEOLIB)/Makefile.geoclaw

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

RIEMANN = $(CLAW)/riemann/src

MODULES = \

SOURCES = \
  ./rpt2_geoclaw.f \
  ./valout.f90 \
  ./src2.f90 \
  ./qinit.f90

RP ?= simple
ifeq ($(RP), simple)
  SOURCES += ./rpn2_shallow_fwave.f90
else ifeq ($(RP), geoclaw)
  SOURCES += $(RIEMANN)/rpn2_geoclaw.f
  SOURCES += $(RIEMANN)/geoclaw_riemann_utils.f
else
  $(error Invalid Riemann solver $(RP))
endif

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)
