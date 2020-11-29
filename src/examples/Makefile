# GCC
FC=gfortran
FFLAGS+=-g -Wall -std=f2008
FFLAGS+=-fcheck=all -pedantic
#FFLAGS+=-fdefault-integer-8

# Intel
#FC=ifort
#FFLAGS=-warn all -g
## Important flag to enable F03's automatic (re-)allocation
#FFLAGS+=-assume realloc_lhs
## Important flag to direct large arrays to the heap
## An alternative is to increase stack memory size ulimit -s unlimited
#FFLAGS+=-heap-arrays

# NAG
#FC=nagfor
#FFLAGS=-C=all


PROGRAMS = source_signature_example \
           coordinate_utilities_test \
           sac_io_example\
           make_isotropic_homogeneous_model_example \
           make_random_model_example \
           make_PREM_model_example \
           checks_test \
           string_utilities_test \
           auto_partitionize_example \
           sac_header_test \
           gll_test \
           interpolation_test \
           coordinates_test \
           ak135_f_test \
           prem_test


# Default target
all: $(PROGRAMS)


## Modules / Dependencies

coordinate_utilities_mod.o: ../misc/coordinate_utilities_mod.f90 \
                            dummy_parameters.o
	$(FC) $(FFLAGS) -c $<

source_signature_mod.o: ../misc/source_signature_mod.f90 \
                        dummy_parameters.o \
                        coordinate_utilities_mod.o
	$(FC) $(FFLAGS) -c $<

string_utilities_mod.o: ../misc/string_utilities_mod.f90
	$(FC) $(FFLAGS) -c $<

sac_header_mod.o: ../sac_io/sac_header_mod.f90 \
                  dummy_parameters.o \
                  string_utilities_mod.o
	$(FC) $(FFLAGS) -c $<

sac_io_mod.o: ../sac_io/sac_io_mod.f90 \
           sac_header_mod.o
	$(FC) $(FFLAGS) -c $<

parser.o: ../parser/parser.f90 \
          dummy_error.o \
          dummy_parameters.o
	$(FC) $(FFLAGS) -c $<

coordinates_mod.o: ../coordinates_mod.f90 \
                   dummy_parameters.o \
                   gll_mod.o
	$(FC) $(FFLAGS) -c $<

gll_mod.o: ../gll_mod.f90 \
       dummy_parameters.o
	$(FC) $(FFLAGS) -c $<

interpolation_mod.o: ../interpolation_mod.f90 \
       dummy_parameters.o
	$(FC) $(FFLAGS) -c $<

prem_mod.o: ../models/prem_mod.F90 \
        dummy_parameters.o
	$(FC) $(FFLAGS) -c $<

checks_mod.o: ../misc/checks_mod.f90 \
          dummy_parameters.o
	$(FC) $(FFLAGS) -c $<

auto_partitionize_mod.o: ../misc/auto_partitionize_mod.f90 \
                         dummy_error.o
	$(FC) $(FFLAGS) -c $<

ak135_f_mod.o: ../models/ak135_f_mod.F90 \
               dummy_parameters.o
	$(FC) $(FFLAGS) -c $<

## Actual programs

source_signature_example: source_signature_example.f90 \
                          dummy_parameters.o \
                          coordinate_utilities_mod.o \
                          source_signature_mod.o
	$(FC) $(FFLAGS) $^ -o $@

coordinate_utilities_test: coordinate_utilities_test.f90 \
                              dummy_parameters.o \
                              coordinate_utilities_mod.o
	$(FC) $(FFLAGS) $^ -o $@

sac_io_example: sac_io_example.f90 \
                dummy_parameters.o \
                string_utilities_mod.o \
                sac_header_mod.o \
                sac_io_mod.o
	$(FC) $(FFLAGS) $^ -o $@

make_isotropic_homogeneous_model_example: make_isotropic_homogeneous_model_example.f90 \
                                dummy_parameters.o \
                                dummy_error.o \
                                parser.o
	$(FC) $(FFLAGS) $^ -o $@

make_random_model_example: make_random_model_example.f90 \
                           dummy_parameters.o \
                           dummy_error.o \
                           parser.o
	$(FC) $(FFLAGS) $^ -o $@

make_PREM_model_example: make_PREM_model_example.f90 \
                         dummy_parameters.o \
                         dummy_error.o \
                         coordinates_mod.o \
                         gll_mod.o \
                         parser.o \
                         prem_mod.o
	$(FC) $(FFLAGS) $^ -o $@

checks_test: checks_test.f90 \
                 dummy_parameters.o \
                 checks_mod.o
	$(FC) $(FFLAGS) $^ -o $@

string_utilities_test: string_utilities_test.f90 \
                       string_utilities_mod.o
	$(FC) $(FFLAGS) $^ -o $@ 

auto_partitionize_example: auto_partitionize_example.f90 \
                           dummy_error.o \
                           auto_partitionize_mod.o
	$(FC) $(FFLAGS) $^ -o $@ 

sac_header_test: sac_header_test.f90 \
                 string_utilities_mod.o \
                 sac_header_mod.o
	$(FC) $(FFLAGS) $^ -o $@ 

gll_test: gll_test.f90 \
          gll_mod.o
	$(FC) $(FFLAGS) $^ -o $@ 

interpolation_test: interpolation_test.f90 \
                    gll_mod.o \
                    interpolation_mod.o
	$(FC) $(FFLAGS) $^ -o $@ 

coordinates_test: coordinates_test.f90 \
                  dummy_parameters.o \
                  gll_mod.o \
                  coordinates_mod.o
	$(FC) $(FFLAGS) $^ -o $@ 

ak135_f_test: ak135_f_test.f90 \
              dummy_parameters.o \
              ak135_f_mod.o
	$(FC) $(FFLAGS) $^ -o $@ 

prem_test: prem_test.f90 \
              dummy_parameters.o \
              prem_mod.o
	$(FC) $(FFLAGS) $^ -o $@ 



.PHONY: clean
clean:
	rm -f *.o
	rm -f *.mod

veryclean: clean
	rm -f $(PROGRAMS)


## Implicit rule to build object files from Fortran sources
%.o %.mod: %.f90
	$(FC) $(FFLAGS) -c $<

