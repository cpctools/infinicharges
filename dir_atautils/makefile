ifeq ("$(wildcard .f2py_f2cmap)","")
$(error .f2py_f2cmap not found)
endif
SOURCES=cubecruncher.f90 parameters.f90 ewald.f90
OBJECTS=$(SOURCES:.f90=.o)
FFLAGS=-c -fPIC
#gfortran fast
GFFFLAGS= -m64 -Ofast -flto -march=native -funroll-loops
#gfortran debug
GFDFLAGS=-g -Wall -Wextra -fbounds-check
#ifort fast
IFFLAGS=-xHost -fast
#ifort debug
IDFLAGS=-O0 -CB

fast: 
		gfortran $(FFLAGS) $(GFFFLAGS) $(SOURCES) 
		f2py --opt="-Ofast -flto -funroll-loops" --arch="-m64 -march=native" -c -m atautils atautils.f90 $(OBJECTS)

fast-sig: 
		gfortran $(FFLAGS) $(GFFFLAGS) $(SOURCES) 
		f2py -m atautils -h atautils.pyf atautils.f90 $(OBJECTS)
		f2py --opt="-Ofast -flto -funroll-loops" --arch="-m64 -march=native" -c atautils.pyf atautils.f90 $(OBJECTS)

fast-gibbs:
		gfortran $(FFLAGS) $(GFFFLAGS) $(SOURCES) 
		scl enable python27 "f2py --opt="-Ofast -flto -funroll-loops" --arch="-m64 -march=native" -c -m atautils atautils.f90 $(OBJECTS)"

debug:
		gfortran $(FFLAGS) $(GFDFLAGS) $(SOURCES) 
		f2py -c -m atautils atautils.f90 $(OBJECTS)

ifort:
		ifort $(FFLAGS) $(IFFLAGS) $(SOURCES) 
		f2py --opt="-fast" --arch="-xHost" --fcompiler=intelem -c -m atautils atautils.f90 $(OBJECTS)

ifort-debug:
		ifort $(FFLAGS) $(IDFLAGS) $(SOURCES) 
		f2py  --fcompiler=intelem -c -m atautils atautils.f90 $(OBJECTS)

profile:
		gfortran $(FFLAGS) $(GFDFLAGS) -p $(SOURCES) 
		gfortran $(GFDFLAGS) -p atautils.f90  maintest.f90 $(OBJECTS) -o maintest.x

testfi:
		ifort $(FFLAGS) $(IFFLAGS) $(SOURCES) 
		ifort $(IFFLAGS) atautils.f90 maintest.f90 $(OBJECTS) -o maintest.x

testfi-debug:
		ifort $(FFLAGS) $(IDFLAGS) $(SOURCES) 
		ifort $(IDFLAGS) atautils.f90 maintest.f90 $(OBJECTS) -o maintest.x

clean:
	    rm -f *.o *.so *.mod *.pyf *.x
