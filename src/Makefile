##################################
# makefile: 
##################################

FC 	= gfortran
FFLAGS= -O3 -g -fbounds-check -fbacktrace

LOADLIBES = 
#LOADLIBES = libarpack_mac.a fftlib.a -llapack

# to link statically we might need something like: 
#FFLAGS += -static
#LOADLIBES += -lblas -lgfortran -lm



SRC = main.f modmain.f skoster.f

OBJ = $(SRC:.f=.o)

.SUFFIXES: .o .f

perovskite: $(OBJ)
	$(FC) -o perovskite $(FFLAGS) $(OBJ) $(LOADLIBES)
	cp perovskite ~/bin/.


# Here are the compile steps
# IMPORTANT: remove unused dependencies after testing all modules.
 
modmain.o: modmain.f  
	 $(FC) -c $(FFLAGS) modmain.f 

skoster.o: skoster.f modmain.o
	 $(FC) -c $(FFLAGS) skoster.f	


main.o: main.f modmain.o skoster.o
	 $(FC) -c $(FFLAGS) main.f	

clean:
	rm -rf *.o *.mod perovskite
  
