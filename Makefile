# Compiler

#FC = epcf90
#FC = gfortran
FC = ifort

# Flags
FCFLAGS = -g -pg -llapack -lblas
#FCFLAGS = -O4 -llapack -lblas

# Objects
OBJECTS = main.o bsimple.o subroutines.o typedefs.o parameters.o mc_toolbox.o spincalc.o

# Program Name
PROGRAM = kmc_spin

# Recipes
$(PROGRAM): $(OBJECTS)
	$(FC) -o $(PROGRAM) *.o $(FCFLAGS)

main.o: main.f90 subroutines.o parameters.o typedefs.o spincalc.o bsimple.o
	$(FC) -c main.f90 $(FCFLAGS)

typedefs.o: typedefs.f90
	$(FC) -c typedefs.f90 $(FCFLAGS)

subroutines.o: subroutines.f90 parameters.o bsimple.o typedefs.o mc_toolbox.o
	$(FC) -c subroutines.f90 $(FCFLAGS)

bsimple.o: bsimple.f90 typedefs.o parameters.o
	$(FC) -c bsimple.f90 $(FCFLAGS)

parameters.o: parameters.f90 typedefs.o
	$(FC) -c parameters.f90 $(FCFLAGS)

mc_toolbox.o: mc_toolbox.f90
	$(FC) -c mc_toolbox.f90 $(FCFLAGS)

spincalc.o: spincalc.f90 subroutines.o parameters.o typedefs.o bsimple.o
	$(FC) -c spincalc.f90 $(FCFLAGS)

clean:
	 rm -f *.o *.mod *.csv *.txt
