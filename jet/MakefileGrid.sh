FFLAGS=  -O3 -r8 -w
#FFLAGS=  -fdefault-double-8

EXENAME= grid.run
FC= ifort 
#FC= gfortran 

OBJS=malla_genera.o

$(EXENAME):  $(OBJS)
	$(FC)  $(OBJS) $(FFLAGS) -o $(EXENAME)

malla_genera.o: malla_genera.f90
	$(FC) $(FFLAGS) -c $*.f90

clean :
	rm *.o
	rm *.mod

