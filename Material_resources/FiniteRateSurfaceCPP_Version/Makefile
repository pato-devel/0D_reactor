#define fortran compiler here
CXX = g++

#define any flags that the compiler needs here
CXXOPTS = -O3

#define bounds checking for diagnostic programs
CXXDIAG = -fbounds-check

#define lists of source and headers
#define all the source files 
SRC =	FiniteRateSurface.cpp \
	lumatrix.cpp \
	FRMArray.cpp

HEADER = FiniteRateSurface.hxx \
	lumatrix.hxx \
	FRMArray.hxx

#build various tool aplications
#all : 
#	@echo Compiling checkders....
#	make checkders
#	@echo Compiling checkinputs....
#	make checkinputs
#	@echo Compiling checksystem....
#	make checksystem
#	@echo Compiling srmodel....
#	make srmodel

#checkders: $(SRC) $(HEADER) main-check-ders.cpp
#	$(CXX) $(CXXDIAG) -o checkders $(SRC) main-check-ders.cpp

#checkinputs: $(SRC) $(HEADER) main-check-inputs.cpp
#	$(CXX) $(CXXDIAG) -o checkinputs $(SRC) main-check-inputs.cpp

#checksystem: $(SRC) $(HEADER) main-check-system.cpp
#	$(CXX) $(CXXDIAG) -o checksystem $(SRC) main-check-system.cpp

srmodel: $(SRC) $(HEADER) main-srmodel.cpp
	$(CXX) $(CXXOPTS) -o srmodel $(SRC) main-srmodel.cpp

gsi:  $(SRC) $(HEADER) gsi-wrapper-class.cpp
	$(CXX) $(CXXOPTS) -o gsiwrapper $(SRC) gsi-wrapper-class.cpp

clean:
	-rm -f -r checkders checkinputs srmodel checksystem gsiwrapper *.o *.exe *.dSYM
