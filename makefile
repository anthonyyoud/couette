OUTDIR		= ./
OBJECTS		= parameters.o ic_bc.o variables.o derivs.o pressure.o \
                  variables.o magnetic.o current.o matrices.o io.o linear.o \
                  nonlinear.o solve.o couette_mod.o
COMPFLAGS	= -O3 -w95 -tpp7 -xW -unroll #-parallel
#COMPFLAGS	= -d0 -CA -CB -CS -CU -CV
LINKFLAGS	= -i_dynamic
TYPE            = implicit
LIBS            = -L $(HOME)/lib/ -lscalapack -lblacsF77init_MPI-LINUX-0 \
                  -lblacs_MPI-LINUX-0 -lblacsCinit_MPI-LINUX-0 -lblas_ref
COMPILER	= mpif77

#-----------------------------------------------------------------------
all:	$(OBJECTS)
	$(COMPILER) $(COMPFLAGS) $(LINKFLAGS) -o\
        $(OUTDIR)couette_mod.out \
        $(OBJECTS) $(LIBS)

#-----------------------------------------------------------------------
parameters.o : parameters.f90
	$(COMPILER) $(COMPFLAGS) -c parameters.f90

#-----------------------------------------------------------------------
ic_bc.o : ic_bc.f90
	$(COMPILER) $(COMPFLAGS) -c ic_bc.f90

#-----------------------------------------------------------------------
variables.o : variables.f90
	$(COMPILER) $(COMPFLAGS) -c variables.f90

#-----------------------------------------------------------------------
io.o : io.f90
	$(COMPILER) $(COMPFLAGS) -c variables.f90 io.f90

#-----------------------------------------------------------------------
derivs.o : derivs.f90
	$(COMPILER) $(COMPFLAGS) -c derivs.f90

#-----------------------------------------------------------------------
pressure.o : pressure.f90
	$(COMPILER) $(COMPFLAGS) -c variables.f90 io.f90 pressure.f90

#-----------------------------------------------------------------------
magnetic.o : magnetic.f90
	$(COMPILER) $(COMPFLAGS) -c io.f90 magnetic.f90

#-----------------------------------------------------------------------
current.o : current.f90
	$(COMPILER) $(COMPFLAGS) -c derivs.f90 io.f90 current.f90

#-----------------------------------------------------------------------
matrices.o : matrices.f90
	$(COMPILER) $(COMPFLAGS) -c matrices.f90

#-----------------------------------------------------------------------
linear.o : linear.f90
	$(COMPILER) $(COMPFLAGS) -c linear.f90

#-----------------------------------------------------------------------
nonlinear.o : nonlinear.f90
	$(COMPILER) $(COMPFLAGS) -c nonlinear.f90

#-----------------------------------------------------------------------
solve.o : solve.f90
	$(COMPILER) $(COMPFLAGS) -c solve.f90

#-----------------------------------------------------------------------
couette_mod.o : couette_mod_$(TYPE).f90
	$(COMPILER) $(COMPFLAGS) -c -o couette_mod.o \
                                    couette_mod_$(TYPE).f90
#-----------------------------------------------------------------------
clean :
	rm -f *.o *.mod *.d *.out *.pc *.pcl *.il
