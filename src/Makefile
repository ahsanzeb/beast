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


SRC77 = main.f modmain.f skoster.f hamiltonian.f readinput.f \
fermi.f 

SRC90 = modgaunt.f90 hubbard.f90 scf.f90 pdos.f90 mixing.f90 \
fixmom.f90

SRC90 += esvar.f90 estatic.f90 mkstrxd.f90 tbeseld.f90

SRC77 += scg.f sylm.f makcg9.f hstra.f tbshfl.f rcnsl0.f
SRC77 += tbmpol.f


OBJ = $(SRC77:.f=.o) $(SRC90:.f90=.o)

####.SUFFIXES: .o .f


beast: $(OBJ)
	$(FC) -o beast $(FFLAGS) $(OBJ) $(LOADLIBES)
	cp beast ~/bin/.


# Here are the compile steps
# IMPORTANT: remove unused dependencies after testing all modules.

modgaunt.o: modgaunt.f90  
	 $(FC) -c $(FFLAGS) modgaunt.f90 

modmain.o: modmain.f  
	 $(FC) -c $(FFLAGS) modmain.f 

hubbard.o: hubbard.f90 modmain.o
	 $(FC) -c $(FFLAGS) hubbard.f90

fixmom.o: fixmom.f90 modmain.o
	 $(FC) -c $(FFLAGS) fixmom.f90


pdos.o: pdos.f90 modmain.o
	 $(FC) -c $(FFLAGS) pdos.f90

mixing.o: mixing.f90 modmain.o esvar.o
	 $(FC) -c $(FFLAGS) mixing.f90

skoster.o: skoster.f modmain.o
	 $(FC) -c $(FFLAGS) skoster.f	


hamiltonian.o: hamiltonian.f modmain.o esvar.o
	 $(FC) -c $(FFLAGS) hamiltonian.f

readinput.o: readinput.f modmain.o
	 $(FC) -c $(FFLAGS) readinput.f

fermi.o: fermi.f modmain.o
	 $(FC) -c $(FFLAGS) fermi.f	

scf.o: scf.f90 modmain.o hamiltonian.o fermi.o modgaunt.o hubbard.o \
       mixing.o fixmom.o estatic.o tbmpol.o tbeseld.o
	 $(FC) -c $(FFLAGS) scf.f90	

main.o: main.f modmain.o skoster.o hamiltonian.o readinput.o fermi.o \
        modgaunt.o hubbard.o scf.o pdos.o estatic.o tbeseld.o
	 $(FC) -c $(FFLAGS) main.f	



####.SUFFIXES: .o .f

# Here are the compile steps
# IMPORTANT: remove unused dependencies after testing all modules.

esvar.o: esvar.f90
	$(FC) -c $(FFLAGS) esvar.f90


estatic.o: estatic.f90 modmain.o esvar.o mkstrxd.o scg.o sylm.o makcg9.o
	$(FC) -c $(FFLAGS) estatic.f90

mkstrxd.o: mkstrxd.f90 hstra.o esvar.o rcnsl0.o
	$(FC) -c $(FFLAGS) mkstrxd.f90


rcnsl0.o: rcnsl0.f esvar.o sylm.o
	$(FC) -c $(FFLAGS) rcnsl0.f


tbshfl.o: tbshfl.f esvar.o
	 $(FC) -c $(FFLAGS) tbshfl.f

tbeseld.o: tbeseld.f90 esvar.o
	 $(FC) -c $(FFLAGS) tbeseld.f90

tbmpol.o: tbmpol.f esvar.o
	 $(FC) -c $(FFLAGS) tbmpol.f


scg.o: scg.f
	 $(FC) -c $(FFLAGS) scg.f

sylm.o: sylm.f
	 $(FC) -c $(FFLAGS) sylm.f

makcg9.o: makcg9.f
	 $(FC) -c $(FFLAGS) makcg9.f

hstra.o: hstra.f esvar.o tbshfl.o
	 $(FC) -c $(FFLAGS) hstra.f





clean:
	rm -rf *.o *.mod beast
  
