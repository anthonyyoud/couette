OUTDIR		= ./
OBJECTS		= parameters.o ic_bc.o variables.o derivs.o stream.o \
                  magnetic.o current.o matrices.o io.o linear.o \
                  nonlinear.o solve.o couette_mod.o
FFLAGS	        = -O2 -w95 -tpp7 -xW -unroll -vec_report0
#FFLAGS	        = -pg -d0 -CA -CB -CS -CU -CV
LINKFLAGS	= -static
COMPILER	= mpif77
LDBLAS          = -L$(BLASHOME)/lib -lblas
LDSCALA         = -L$(SCALAPACKHOME)/lib -lscalapack
#LDSCALA         = -L/work/n8049290/mpi/SCALAPACK -lscalapack_ifc8_gcc
LDBLACS         = -L$(BLACSHOME)/lib -lblacsF77init -lblacs -lblacsCinit

LIBS            = $(LDSCALA) $(LDBLACS) $(LDBLAS)
COMPFLAGS       = $(FFLAGS)
#-----------------------------------------------------------------------
all:	$(OBJECTS)
	$(COMPILER) $(COMPFLAGS) $(LINKFLAGS) -o \
        $(OUTDIR)couette_mod.out \
        $(OBJECTS) $(LIBS)

#-----------------------------------------------------------------------
parameters.o : parameters.f90
	$(COMPILER) $(COMPFLAGS) -c parameters.f90

#-----------------------------------------------------------------------
ic_bc.o : ic_bc.f90
	$(COMPILER) $(COMPFLAGS) -c ic_bc.f90

#-----------------------------------------------------------------------
io.o : io.f90
	$(COMPILER) $(COMPFLAGS) -c io.f90

#-----------------------------------------------------------------------
variables.o : variables.f90
	$(COMPILER) $(COMPFLAGS) -c variables.f90

#-----------------------------------------------------------------------
derivs.o : derivs.f90
	$(COMPILER) $(COMPFLAGS) -c derivs.f90

#-----------------------------------------------------------------------
stream.o : stream.f90
	$(COMPILER) $(COMPFLAGS) -c stream.f90

#-----------------------------------------------------------------------
magnetic.o : magnetic.f90
	$(COMPILER) $(COMPFLAGS) -c magnetic.f90

#-----------------------------------------------------------------------
current.o : current.f90
	$(COMPILER) $(COMPFLAGS) -c current.f90

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
couette_mod.o : couette_mod.f90
	$(COMPILER) $(COMPFLAGS) -c couette_mod.f90
#-----------------------------------------------------------------------
clean :
	rm -f *.o *.mod *.d *.out *.pc *.pcl *.il
