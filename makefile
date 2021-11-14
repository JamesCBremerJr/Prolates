.phony:  all PROGRAM clean

# Choose the compiler and set compiler options
OMP_NUM_THREADS      := $(shell nproc)
export OMP_NUM_THREADS

OMP_STACKSIZE        := 1024M
export OMP_STACKSIZE

FCOMP     =  gfortran
FOPTS     =  -fopenmp -Ofast -march=native -w -fexternal-blas -fallow-argument-mismatch 
#FOPTS    +=  -freal-8-real-10 
#FOPTS    +=  -fdefault-real-8

# specify BLAS and LAPACK library
LDOPT     =  -lblas -llapack
#LDOPT     = OpenBLAS/libopenblas.a
#LDOPT     =  FLAME/libflame.a BLIS/libblis.a

FCOMPQUAD = $(FCOMP)
FC        = $(FCOMP)
FFLAGS    = $(FOPTS)

CC        = gcc
COPTS     = -O3  -I./include

export FC FFLAGS

# Set the list of programs to compile

PROGRAMS = test_prolates test_chebyshev test_prolates_phase test_prolates_exps            \
           test_spheroidal test_spheroidal_phase experiment1 experiment2 

# Compile all of the test programs and the library
all	      	             : clean $(PROGRAMS) 


# List the dependencies for each module's test program

EXPERIMENTS1_FILES           = utils.o                                                    \
                               chebyshev.o                                                \
                               prolates10.o                                               \
                               prolates.o                                                 \
                               prolates_phase.o                                           \
                               prolates_exps.o

EXPERIMENTS2_FILES           = utils.o                                                    \
                               chebyshev.o                                                \
                               prolates.o                                                 \
                               prolates_phase.o                                           \
                               prolates_exps.o

SPHEROIDAL_PHASE_FILES       = utils.o                                                    \
                               chebyshev.o                                                \
                               spheroidal.o                                               \
                               spheroidal_phase.o

SPHEROIDAL_FILES             = utils.o                                                    \
                               spheroidal.o                                               \

PROLATES_EXPS_FILES          = utils.o                                                    \
                               chebyshev.o                                                \
                               prolates.o                                                 \
                               prolates_phase.o                                           \
                               prolates_exps.o


PROLATES_PHASE_FILES         = utils.o                                                    \
                               chebyshev.o                                                \
                               prolates.o                                                 \
                               prolates_phase.o

CHEBYSHEV_FILES              = utils.o                                                    \
                               chebyshev.o


PROLATES_FILES               = utils.o                                                    \
                               prolates.o                                                 \



test_prolates_phase2.o       : $(PROLATES_PHASE2_FILES) test_prolates_phase2.f90
test_prolates_phase2         : $(PROLATES_PHASE2_FILES) test_prolates_phase2.o


experiments1                 : $(EXPERIMENTS1_FILES) 

experiments2                 : $(EXPERIMENTS2_FILES) 

test_spheroidal_phase.o      : $(SPHEROIDAL_PHASE_FILES) test_spheroidal_phase.f90
test_spheroidal_phase        : $(SPHEROIDAL_PHASE_FILES) test_spheroidal_phase.o

test_spheroidal.o            : $(SPHEROIDAL_FILES) test_spheroidal.f90
test_spheroidal              : $(SPHEROIDAL_FILES) test_spheroidal.o

test_prolates_exps.o         : $(PROLATES_EXPS_FILES) test_prolates_exps.f90
test_prolates_exps           : $(PROLATES_EXPS_FILES) test_prolates_exps.o

test_prolates_phase.o        : $(PROLATES_PHASE_FILES) test_prolates_phase.f90
test_prolates_phase          : $(PROLATES_PHASE_FILES) test_prolates_phase.o

test_chebyshev.o             : $(CHEBYSHEV_FILES) test_chebyshev.f90
test_chebyshev               : $(CHEBYSHEV_FILES) test_chebyshev.o

test_prolates.o              : $(PROLATES_FILES) test_prolates.f90
test_prolates                : $(PROLATES_FILES) test_prolates.o



# Setup the general compilation rules

%		: %.o
	$(FCOMP) $(FOPTS) -o $@ $^ $(LDOPT)
	@echo  
	@echo 
	@echo "---------[ $@     ]--------------------------------------------"
	@echo 
	@./$@
	@echo 
	@echo "--------------------------------------------------------------------------"
	@echo 


%.o		: %.f90
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.f
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.cpp
	$(CPPCOMP) -c $(CPPOPTS)  $<

%.o		: %.c
	$(CC) -c $(COPTS)  $<

.PHONY: clean

clean:
	rm -f *.o *.mod *~ fort.* a.out *.a
	rm -f $(PROGRAMS)
	rm -f *.dat
	rm -f *.py
	rm -f *.pdf
	rm -f table*.tex

