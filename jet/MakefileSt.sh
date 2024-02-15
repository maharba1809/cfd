
FFLAGS=  -O3 -r8 -w
EXENAME= st.run
FC= ifort 

OBJS=estadia.o

$(EXENAME):  $(OBJS)
	$(FC)  $(OBJS) $(FFLAGS) -o $(EXENAME)

estadia.o: estadia.f90
	$(FC) $(FFLAGS) -c $*.f90

clean :
	rm *.o
	rm *.mod

