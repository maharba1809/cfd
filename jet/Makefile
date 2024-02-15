FFLAGS=  -O3 -r8 -w
#FFLAGS=  -fdefault-double-8
 
EXENAME= jetrun
FC= ifort 
#FC= gfortran 

OBJS=\
flow_mod.o\
aleatorio.o\
derivada.o\
deriv_vel.o\
divergencia.o\
ehecatl.o\
esponja.o\
estadisticas.o\
filtro2b.o\
flow_alloc.o\
flujos.o\
frontera_z.o\
frontera_xss.o\
frontera_y.o\
gravar.o\
herramientas2.o\
iniconst.o\
jacob_malla.o\
leer.o\
rungek2.o\
sgdm.o\
source.o\
tstep.o\
viscosidad.o\

$(EXENAME):  $(OBJS)
	$(FC)  $(OBJS) $(FFLAGS) -o $(EXENAME)
flow_mod.o: flow_mod.f90
	$(FC) $(FFLAGS) -c $*.f90
aleatorio.o: aleatorio.f90
	$(FC) $(FFLAGS) -c $*.f90
derivada.o: derivada.f90
	$(FC) $(FFLAGS) -c $*.f90
deriv_vel.o: deriv_vel.f90
	$(FC) $(FFLAGS) -c $*.f90
divergencia.o: divergencia.f90
	$(FC) $(FFLAGS) -c $*.f90
ehecatl.o: ehecatl.f90
	$(FC) $(FFLAGS) -c $*.f90
esponja.o: esponja.f90
	$(FC) $(FFLAGS) -c $*.f90
estadisticas.o:estadisticas.f90
	$(FC) $(FFLAGS) -c $*.f90
filtro2b.o: filtro2b.f90
	$(FC) $(FFLAGS) -c $*.f90
flow_alloc.o: flow_alloc.f90
	$(FC) $(FFLAGS) -c $*.f90
flujos.o: flujos.f90
	$(FC) $(FFLAGS) -c $*.f90
frontera_z.o: frontera_z.f90
	$(FC) $(FFLAGS) -c $*.f90
frontera_xss.o: frontera_xss.f90
	$(FC) $(FFLAGS) -c $*.f90
frontera_y.o: frontera_y.f90
	$(FC) $(FFLAGS) -c $*.f90
gravar.o: gravar.f90
	$(FC) $(FFLAGS) -c $*.f90
herramientas2.o: herramientas2.f90
	$(FC) $(FFLAGS) -c $*.f90
iniconst.o: iniconst.f90
	$(FC) $(FFLAGS) -c $*.f90
jacob_malla.o: jacob_malla.f90
	$(FC) $(FFLAGS) -c $*.f90
leer.o: leer.f90
	$(FC) $(FFLAGS) -c $*.f90
rungek2.o: rungek2.f90
	$(FC) $(FFLAGS) -c $*.f90
sgdm.o: sgdm.f90
	$(FC) $(FFLAGS) -c $*.f90
source.o: source.f90
	$(FC) $(FFLAGS) -c $*.f90
tstep.o: tstep.f90
	$(FC) $(FFLAGS) -c $*.f90
viscosidad.o: viscosidad.f90
	$(FC) $(FFLAGS) -c $*.f90
clean :
	rm *.mod
	/bin/rm *.o
	
