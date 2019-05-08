# Choose the compiler and set compiler options

FCOMP    = gfortran-7
FOPTS    = -w -fopenmp -O3 -march=native -fexternal-blas -fdefault-integer-8 
LDOPT    = -llapack -lblas -lfftw3

PYTHON   = python

OMP_STACKSIZE        := 4096M
export OMP_STACKSIZE

OMP_NUM_THREADS      := $(shell nproc)
export OMP_NUM_THREADS


OMP_NESTED           := true
export OMP_NESTED

FC     = $(FCOMP)
FFLAGS = $(FOPTS)
export FC FFLAGS

# Set the list of programs to compile

PROGRAMS               = test_chebyshev test_tensor test_adapquad create_boundchi          \
                         test_linearode  test_riccati test_prolates test_prolates_window   \
                         test_prolates_exps create_prolexps test_prolates_fast             \
                         plot_chi plot_A plot_phases test_chiexp test_pseval



###################################################################################

TEST_PSEVAL_FILES      = utils.o                                 \
                         chebyshev.o                             \
                         tensor.o                                \
                         riccati.o                               \
                         linearode.o                             \
                         prolates.o                              \
                         prolates_window.o                       \
                         prolates_exps.o                         \
                         prolexps.o                              \
                         prolates_fast.o

TEST_CHIEXP_FILES      = utils.o                                 \
                         chebyshev.o                             \
                         tensor.o                                \
                         riccati.o                               \
                         linearode.o                             \
                         prolates.o                              \
                         prolates_window.o                       \
                         prolexps.o

PLOTPHASES_FILES       = utils.o                                 \
                         chebyshev.o                             \
                         tensor.o                                \
                         riccati.o                               \
                         linearode.o                             \
                         prolates.o                              \
                         prolates_window.o                       \
                         prolates_exps.o


PLOTA_FILES            = utils.o                                 \
                         chebyshev.o                             \
                         tensor.o                                \
                         riccati.o                               \
                         linearode.o                             \
                         prolates.o                              \
                         prolates_window.o                       \
                         prolates_exps.o

PLOTCHI_FILES          = utils.o                                 \
                         chebyshev.o                             \
                         tensor.o                                \
                         riccati.o                               \
                         linearode.o                             \
                         prolates.o                              \
                         prolates_window.o                       \
                         prolates_exps.o

TEST_PROLATES_FAST_FILES = utils.o                               \
                         chebyshev.o                             \
                         tensor.o                                \
                         riccati.o                               \
                         linearode.o                             \
                         prolates.o                              \
                         prolates_window.o                       \
                         prolexps.o                              \
                         prolates_fast.o

CREATE_PROLEXPS_FILES  = utils.o                                 \
                         chebyshev.o                             \
                         tensor.o                                \
                         riccati.o                               \
                         linearode.o                             \
                         prolates.o                              \
                         prolates_window.o                       \
                         prolates_exps.o                         \
                         create_prolexps.o

PROLATES_EXPS_FILES    = utils.o                                 \
                         chebyshev.o                             \
                         tensor.o                                \
                         riccati.o                               \
                         linearode.o                             \
                         prolates.o                              \
                         prolates_window.o                       \
                         prolates_exps.o

PROLATES_WINDOW_FILES  = utils.o                                 \
                         chebyshev.o                             \
                         tensor.o                                \
                         riccati.o                               \
                         linearode.o                             \
                         prolates.o                              \
                         prolates_window.o


PROLATES_FILES         = utils.o                                 \
                         prolates.o

RICCATI_FILES          = utils.o                                 \
                         chebyshev.o                             \
                         riccati.o

LINEARODE_FILES        = utils.o                                  \
                         chebyshev.o                              \
                         linearode.o

CREATE_BOUNDCHI_FILES  = utils.o                                 \
			 chebyshev.o                             \
			 adapquad.o

ADAPQUAD_FILES         = utils.o                                  \
                         adapquad.o

TENSOR_FILES           = utils.o                                  \
                         chebyshev.o                              \
                         tensor.o

CHEBYSHEV_FILES        = utils.o                                  \
                         chebyshev.o


###################################################################################

all	                : clean $(PROGRAMS) 

test_pseval.o           : $(TEST_PSEVAL_FILES) test_pseval.f90
test_pseval             : $(TEST_PSEVAL_FILES) test_pseval

test_chiexp.o           : $(TEST_CHIEXP_FILES) test_chiexp.f90
test_chiexp             : $(TEST_CHIEXP_FILES) test_chiexp.o

plot_phases.o           : $(PLOTPHASES_FILES) plot_phases.f90
plot_phases             : $(PLOTPHASES_FILES) plot_phases.o

plot_A.o                : $(PLOTA_FILES) plot_A.f90
plot_A                  : $(PLOTA_FILES) plot_A.o

plot_chi.o              : $(PLOTCHI_FILES) plot_chi.f90
plot_chi                : $(PLOTCHI_FILES) plot_chi.o

test_prolates_fast.o    : $(TEST_PROLATES_FAST_FILES) test_prolates_fast.f90
test_prolates_fast      : $(TEST_PROLATES_FAST_FILES) test_prolates_fast.o

create_prolexps.o       : $(CREATE_PROLEXPS_FILES) create_prolexps.f90
create_prolexps         : $(CREATE_PROLEXPS_FILES) create_prolexps.o

test_prolates_exps.o    : $(PROLATES_EXPS_FILES) test_prolates_exps.f90
test_prolates_exps      : $(PROLATES_EXPS_FILES) test_prolates_exps.o

test_prolates_window.o  : $(PROLATES_WINDOW_FILES) test_prolates_window.f90
test_prolates_window    : $(PROLATES_WINDOW_FILES) test_prolates_window.o

test_prolates.o         : $(PROLATES_FILES) test_prolates.f90
test_prolates           : $(PROLATES_FILES) test_prolates.o

test_riccati.o          : $(RICCATI_FILES) test_riccati.f90
test_riccati            : $(RICCATI_FILES) test_riccati.o

test_linearode.o        : $(LINEARODE_FILES) test_linearode.f90
test_linearode          : $(LINEARODE_FILES) test_linearode.o

create_boundchi.o       : $(CREATE_BOUNDCHI_FILES) create_boundchi.f90
create_boundchi         : $(CREATE_BOUNDCHI_FILES) create_boundchi.o

test_adapquad.o         : $(ADAPQUAD_FILES) test_adapquad.f90
test_adapquad           : $(ADAPQUAD_FILES) test_adapquad.o

test_tensor.o           : $(TENSOR_FILES) test_tensor.f90
test_tensor             : $(TENSOR_FILES) test_tensor.o

test_chebyshev.o        : $(CHEBYSHEV_FILES) test_chebyshev.f90
test_chebyshev          : $(CHEBYSHEV_FILES) test_chebyshev.o



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

.PHONY: clean

clean:
	rm -f *.o *.mod *~ fort.* a.out *.a
	rm -f $(PROGRAMS)
	rm -f *.py
	rm -f *.pdf
	rm -f *.tex

