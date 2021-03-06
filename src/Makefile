# |================|
# | SES3D-Makefile |
# |================|

# FIXME: Directories and sub-directories are not created if needed
# binaries will stored here
BINDIR=../bin/
# object files will stored here
OBJDIR=../obj/
# Fortran module files will stored here
MODDIR=../mod/


## The Fortran MPI wrapper script if not given as environment variable
ifndef MPI_FC
    MPI_FC=mpif90
endif


### GCC ###
# Fortran module files will be stored here
MDIR=-J $(MODDIR)
## compiler flags, e.g. -I<include dir> if you have headers in a nonstandard directory
#FFLAGS+=-I~/.local/include
# flags for maximum performance
FFLAGS+=-O3
FFLAGS+=-mtune=native
FFLAGS+=-ffast-math -funroll-loops -ftree-vectorize
# flags for debugging
#FFLAGS+= -g -Wall -std=f2008
#FFLAGS+= -fcheck=all -pedantic


#### NAG ###
#FFLAGS=-g -C=all -mismatch
## Fortran module files will stored here
#MDIR=-mdir $(MODDIR)


#### Intel ###
#FC=ifort
#FFLAGS=-g -warn all -g
## Important flag to enable F03's automatic (re-)allocation
#FFLAGS+=-assume realloc_lhs
## Important flag to direct large arrays to the heap
## An alternative is to increase stack memory size ulimit -s unlimited
#FFLAGS+=-heap-arrays
## Fortran module files will stored here
#MDIR=-module $(MODDIR)


## preprocessor flags
# If set INCLUDE 'mpif.h' will be use instead of USE mpi
#CPPFLAGS+=-DMPI_INCLUDE
# Activates netCDF output
#CPPFLAGS+=-DNetCDF


# Linker flags, e.g. -L<lib dir> if you have libraries in a nonstandard directory
#LDFLAGS+=-L~/.local/lib


## libraries to be passed to the linker, e.g. -l<library>
# Compiles with netCDF libraries
#LIBS+=-lnetcdff -lnetcdf



# ============================================================================

# Implicit make rules
${OBJDIR}%.o : %.f90
	$(MPI_FC) $(FFLAGS) -I $(MODDIR) $(MDIR)  -c $< -o $@ $(LIBS)
# with preprocessing
${OBJDIR}%.o : %.F90
	$(MPI_FC) $(FFLAGS) -I $(MODDIR) $(MDIR) $(CPPFLAGS) -c $< -o $@ $(LIBS)


# Object files
# Unfortunately, the order matters
OBJ = $(OBJDIR)misc/error_mod.o \
      $(OBJDIR)misc/string_utilities_mod.o \
      $(OBJDIR)parameters.o \
      $(OBJDIR)misc/checks_mod.o \
      $(OBJDIR)misc/coordinate_utilities_mod.o \
      $(OBJDIR)misc/source_signature_mod.o \
      $(OBJDIR)misc/auto_partitionize_mod.o \
      $(OBJDIR)models/homogeneous.o \
      $(OBJDIR)models/prem_mod.o \
      $(OBJDIR)models/ak135_f_mod.o \
      $(OBJDIR)sac_io/sac_header_mod.o \
      $(OBJDIR)sac_io/sac_io_mod.o \
      $(OBJDIR)gll_mod.o \
      $(OBJDIR)coordinates_mod.o \
      $(OBJDIR)interpolation_mod.o \
      $(OBJDIR)parser/parser.o \
      $(OBJDIR)parser/bcast.o \
      $(OBJDIR)config/grid.o \
      $(OBJDIR)config/general.o \
      $(OBJDIR)config/model.o \
      $(OBJDIR)config/time.o \
      $(OBJDIR)config/config.o \
      $(OBJDIR)mpi_parallel_io.o \
      $(OBJDIR)netcdf_parallel_io.o \
      $(OBJDIR)communicate_fields.o \
      $(OBJDIR)geometric_paras.o \
      $(OBJDIR)model_paras.o \
      $(OBJDIR)elastic_vars.o \
      $(OBJDIR)SOURCE/ses3d_grad.o \
      $(OBJDIR)volume_snapshot.o \
      $(OBJDIR)config/point.o \
      $(OBJDIR)config/channel_mod.o \
      $(OBJDIR)output_raw.o \
      $(OBJDIR)output_netcdf.o \
      $(OBJDIR)evolution.o \
      $(OBJDIR)config/receiver.o \
      $(OBJDIR)config/source.o \


# Makefile targets

.PHONY: ses3d-nt
ses3d-nt: $(BINDIR)ses3d-nt


$(BINDIR)ses3d-nt: $(OBJ) ses3d-nt.F90
	$(MPI_FC) $(FFLAGS) $(CPPFLAGS) -I $(MODDIR) $^ -o $@ $(LDFLAGS) $(LIBS)


.PHONY: clean
clean:
	rm -f $(MODDIR)*.mod
	rm -f $(OBJ)
