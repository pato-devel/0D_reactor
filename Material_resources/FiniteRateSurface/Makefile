#define fortran compiler here
F90C = gfortran

#define any flags that the compiler needs here
FOPTS = -O3

#define bounds checking for diagnostic programs
FDIAG = -fbounds-check

#build various tool aplications
all : 
	@echo Compiling checkders....
	make checkders
	@echo Compiling checkinputs....
	make checkinputs
	@echo Compiling checksystem....
	make checksystem
	@echo Compiling srmodel....
	make srmodel

checkders: surfacemod.F FiniteRateSurface.F main-check-ders.F
	$(F90C) $(FDIAG) -o checkders surfacemod.F FiniteRateSurface.F main-check-ders.F

checkinputs: surfacemod.F FiniteRateSurface.F main-check-inputs.F
	$(F90C) $(FDIAG) -o checkinputs surfacemod.F FiniteRateSurface.F main-check-inputs.F

checksystem: surfacemod.F FiniteRateSurface.F readlewis.F main-check-system.F
	$(F90C) $(FDIAG) -o checksystem surfacemod.F FiniteRateSurface.F readlewis.F main-check-system.F

srmodel: surfacemod.F FiniteRateSurface.F readlewis.F lumatrix.F main-srmodel.F equilibrium.F
	$(F90C) $(FOPTS) -o srmodel surfacemod.F FiniteRateSurface.F readlewis.F lumatrix.F equilibrium.F main-srmodel.F

clean:
	-rm -f checkders checkinputs srmodel checksystem *.o *.exe *.mod
