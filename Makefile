# $Id$
#----------------------------------------------------------------------------

OBJECT = couette_mod
OBJS = couette_mod.o current.o derivs.o ic_bc.o io.o linear.o magnetic.o \
       matrices.o nonlinear.o parameters.o solve.o stream.o variables.o
FC = sunf95
FFLAGS = -fast
LDFLAGS = -llapack -lblas
#-----------------------------------------------------------------------
%.o : %.f90
	$(FC) $(FFLAGS) -c $*.f90

all:	$(OBJECT)

clean :
	rm -f $(OBJECT) *.o *.mod
#-----------------------------------------------------------------------
$(OBJECT): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) $(LDFLAGS) -o $@

ic_bc.o: parameters.o
variables.o: parameters.o ic_bc.o derivs.o
io.o: parameters.o variables.o ic_bc.o derivs.o
derivs.o: parameters.o 
stream.o: parameters.o ic_bc.o
magnetic.o: parameters.o ic_bc.o derivs.o
current.o: parameters.o variables.o ic_bc.o derivs.o
matrices.o: parameters.o variables.o ic_bc.o
linear.o: parameters.o variables.o derivs.o ic_bc.o
nonlinear.o: parameters.o variables.o derivs.o ic_bc.o
solve.o: parameters.o variables.o ic_bc.o
couette_mod.o: parameters.o stream.o matrices.o io.o ic_bc.o variables.o \
  linear.o nonlinear.o solve.o magnetic.o current.o
