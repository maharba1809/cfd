
FFLAGS=  -O3 -r8 -w
#FFLAGS=  -fdefault-real-16
EXENAME= tp.run
#FC= gfortran 
FC= ifort 

OBJS=tp_maker.o

$(EXENAME):  $(OBJS)
	$(FC)  $(OBJS) $(FFLAGS) -o $(EXENAME)

tp_maker.o: tp_maker.f90
	$(FC) $(FFLAGS) -c $*.f90

clean :
	rm *.o
	rm *.mod

