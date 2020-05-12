##################################
# makefile: 
##################################

FC 	= gfortran
FFLAGS= -O3 -g -fbounds-check -fbacktrace

LOADLIBES = 
#LOADLIBES = libarpack_mac.a fftlib.a -llapack

LOADLIBES = -llapack -lblas

# to link statically we might need something like: 
#FFLAGS += -static
#LOADLIBES += -lblas -lgfortran -lm

SRC90 = esvar.f90
SRC77 = modmain.f 


SRC90 = estatic.f90 mkstrxd.f90 

SRC77 += scg.f sylm.f makcg9.f hstra.f tbshfl.f rcnsl0.f

SRC90 += test.f90

OBJ = $(SRC90:.f90=.o)    $(SRC77:.f=.o) 

####.SUFFIXES: .o .f


test: $(OBJ)
	$(FC) -o test $(FFLAGS) $(OBJ) $(LOADLIBES)


# Here are the compile steps
# IMPORTANT: remove unused dependencies after testing all modules.

esvar.o: esvar.f90
	$(FC) -c $(FFLAGS) esvar.f90


estatic.o: estatic.f90 modmain.o esvar.o mkstrxd.o scg.o sylm.o makcg9.o
	$(FC) -c $(FFLAGS) estatic.f90

mkstrxd.o: mkstrxd.f90 hstra.o esvar.o
	$(FC) -c $(FFLAGS) mkstrxd.f90


rcnsl0.o: rcnsl0.f esvar.o
	$(FC) -c $(FFLAGS) rcnsl0.f


tbshfl.o: tbshfl.f esvar.o
	 $(FC) -c $(FFLAGS) tbshfl.f


scg.o: scg.f
	 $(FC) -c $(FFLAGS) scg.f

sylm.o: sylm.f
	 $(FC) -c $(FFLAGS) sylm.f

makcg9.o: makcg9.f
	 $(FC) -c $(FFLAGS) makcg9.f

hstra.o: hstra.f esvar.o tbshfl.o
	 $(FC) -c $(FFLAGS) hstra.f



modmain.o: modmain.f  
	 $(FC) -c $(FFLAGS) modmain.f 


test.o: test.f90 modmain.o estatic.o esvar.o
	 $(FC) -c $(FFLAGS) test.f90 

clean:
	rm -rf *.o *.mod test
  
