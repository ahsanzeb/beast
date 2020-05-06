##################################
# makefile: 
##################################

FC 	= gfortran
FFLAGS= -O3 -g -fbounds-check -fbacktrace

LOADLIBES = 
#LOADLIBES = libarpack_mac.a fftlib.a -llapack

LOADLIBES = -llapack

# to link statically we might need something like: 
#FFLAGS += -static
#LOADLIBES += -lblas -lgfortran -lm



SRC = main.f modmain.f skoster.f hamiltonian.f readinput.f \
fermi.f modgaunt.f90 hubbard.f90 scf.f90 pdos.f90 mixing.f90 \
fixmom.f90 

OBJ = $(SRC:.f=.o)

.SUFFIXES: .o .f

perovskite: $(OBJ)
	$(FC) -o perovskite $(FFLAGS) $(OBJ) $(LOADLIBES)
	cp perovskite ~/bin/.


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

mixing.o: mixing.f90 modmain.o
	 $(FC) -c $(FFLAGS) mixing.f90

skoster.o: skoster.f modmain.o
	 $(FC) -c $(FFLAGS) skoster.f	


hamiltonian.o: hamiltonian.f modmain.o
	 $(FC) -c $(FFLAGS) hamiltonian.f

readinput.o: readinput.f modmain.o
	 $(FC) -c $(FFLAGS) readinput.f

fermi.o: fermi.f modmain.o
	 $(FC) -c $(FFLAGS) fermi.f	

scf.o: scf.f90 modmain.o hamiltonian.o fermi.o modgaunt.o hubbard.o mixing.o fixmom.o
	 $(FC) -c $(FFLAGS) scf.f90	

main.o: main.f modmain.o skoster.o hamiltonian.o readinput.o fermi.o modgaunt.o hubbard.o scf.o pdos.o
	 $(FC) -c $(FFLAGS) main.f	

clean:
	rm -rf *.o *.mod perovskite
  
