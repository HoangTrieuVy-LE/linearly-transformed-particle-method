############################ -*- Mode: Makefile -*- ###########################
## Makefile --- MPI Hands-On session: Poisson's equation
##
## Author  : Isabelle Dupays (CNRS/IDRIS) <Isabelle.Dupays@idris.fr>
###############################################################################
SHELL := /bin/bash
# Include Machine-specific Parameters
include ../../arch/make_inc

OBJS = module_params.o module_parallel_mpi.o module_compute.o poisson.o
OBJS1 = read.o
OBJS2 = exact_solution.o
OBJS3 = module_params.o module_parallel_mpi.o module_compute.o poisson_nonblocking.o

# Implicit Rules
.SUFFIXES: .o .f90
.f90.o:
	$(CF95) -c  $(FFLAGS_TP8) $<

# Makefile Rules
default: poisson

all: poisson read exact_solution

poisson: $(OBJS)
	$(CF95) $(LDFLAGS_TP8) -o $@ $(OBJS)

poisson_nonblocking: $(OBJS3)
	$(CF95) $(LDFLAGS_TP8) -o $@ $(OBJS3)

read:  $(OBJS1)
	$(CF95) $(LDFLAGS_TP8) -o $@ $(OBJS1)

exact_solution: $(OBJS2)
	$(CF95) $(LDFLAGS_TP8) -o $@ $(OBJS2)

verification: fort.11 exact_solution
	$(EXEC_TP8_VERIF) ./exact_solution

fort.11: read data.dat
	$(EXEC_TP8_VERIF) ./read

exe data.dat: poisson poisson.data
	rm -f data.dat
	$(MPIEXEC_TP8) ./poisson

clean:
	rm -f $(OBJS) $(OBJS1) $(OBJS2) $(OBJS3) poisson poisson_nonblocking exact_solution read *.mod core data.dat fort.10 fort.11
