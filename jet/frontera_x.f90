!__________________________________________________
       SUBROUTINE frontera_x ()
!__________________________________________________

      use dimensiones
      use velocidades
      use variables
      use deltas
      use consadim
      use flow
      use fronterapl
      use jacobtools
      use derivtools
      use consderper
      use consdernper

      IMPLICIT NONE
      integer :: i,j,k,l,m
      real rii,pii,uii,vii,wii,tii,uci
      real rii_1,pii_1,uii_1,vii_1,wii_1
      real rii_2,pii_2,uii_2,vii_2,wii_2
      real rii_3,pii_3,uii_3,vii_3,wii_3
      real rii_4,pii_4,uii_4,vii_4,wii_4
      real rif,pif,uif,vif,wif,tif,ucf
      real rif_1,pif_1,uif_1,vif_1,wif_1
      real rif_2,pif_2,uif_2,vif_2,wif_2
      real rif_3,pif_3,uif_3,vif_3,wif_3
      real rif_4,pif_4,uif_4,vif_4,wif_4
      real dudxi,drdxi,dvdxi,dpdxi,dwdxi
      real dudxf,drdxf,dvdxf,dpdxf,dwdxf
      real p_inf,xmachi,xmachf
      real pl_l1,pl_l2,pl_l3,pl_l4,pl_l5
      real pl_d1,pl_d2,pl_d3,pl_d4,pl_d5

!______________________________________________
!    z=1
!______________________________________________

      DO k=1,nz
       DO j=1,ny

        uii=u(1,j,k,1)
        vii=u(1,j,k,2)
        wii=u(1,j,k,3)
        pii=pres(1,j,k)
        rii=um(1,j,k,1)
        tii=pii/rii

        uii_1=u(2,j,k,1)
        vii_1=u(2,j,k,2)
        wii_1=u(2,j,k,3)
        pii_1=pres(2,j,k)
        rii_1=um(2,j,k,1)

        uii_2=u(3,j,k,1)
        vii_2=u(3,j,k,2)
        wii_2=u(3,j,k,3)
        pii_2=pres(3,j,k)
        rii_2=um(3,j,k,1)

        uii_3=u(4,j,k,1)
        vii_3=u(4,j,k,2)
        wii_3=u(4,j,k,3)
        pii_3=pres(4,j,k)
        rii_3=um(4,j,k,1)

        uii_4=u(5,j,k,1)
        vii_4=u(5,j,k,2)
        wii_4=u(5,j,k,3)
        pii_4=pres(5,j,k)
        rii_4=um(5,j,k,1)

        dudxi=jbn(1,j,k,1)*deltax*((-25./12.)*uii+4.*uii_1-   &
              3.*uii_2+(4./3.)*uii_3-(1./4.)*uii_4)
        dvdxi=jbn(1,j,k,1)*deltax*((-25./12.)*vii+4.*vii_1-   &
              3.*vii_2+(4./3.)*vii_3-(1./4.)*vii_4)
        dwdxi=jbn(1,j,k,1)*deltax*((-25./12.)*wii+4.*wii_1-   &
              3.*wii_2+(4./3.)*wii_3-(1./4.)*wii_4)
        dpdxi=jbn(1,j,k,1)*deltax*((-25./12.)*pii+4.*pii_1-   &
              3.*pii_2+(4./3.)*pii_3-(1./4.)*pii_4)
        drdxi=jbn(1,j,k,1)*deltax*((-25./12.)*rii+4.*rii_1-   &
              3.*rii_2+(4./3.)*rii_3-(1./4.)*rii_4)


        uif=u(nx,j,k,1)
        vif=u(nx,j,k,2)
        wif=u(nx,j,k,3)
        pif=pres(nx,j,k)
        rif=um(nx,j,k,1)
        tif=pii/rii

        uif_1=u(nx-1,j,k,1)
        vif_1=u(nx-1,j,k,2)
        wif_1=u(nx-1,j,k,3)
        pif_1=pres(nx-1,j,k)
        rif_1=um(nx-1,j,k,1)

        uif_2=u(nx-2,j,k,1)
        vif_2=u(nx-2,j,k,2)
        wif_2=u(nx-2,j,k,3)
        pif_2=pres(nx-2,j,k)
        rif_2=um(nx-2,j,k,1)

        uif_3=u(nx-3,j,k,1)
        vif_3=u(nx-3,j,k,2)
        wif_3=u(nx-3,j,k,3)
        pif_3=pres(nx-3,j,k)
        rif_3=um(nx-3,j,k,1)

        uif_4=u(nx-4,j,k,1)
        vif_4=u(nx-4,j,k,2)
        wif_4=u(nx-4,j,k,3)
        pif_4=pres(nx-4,j,k)
        rif_4=um(nx-4,j,k,1)

        dudxf=jbn(nx,j,k,1)*deltax*((25./12.)*uif-4.*uif_1+   &
              3.*uif_2-(4./3.)*uif_3+(1./4.)*uif_4)
        dvdxf=jbn(nx,j,k,1)*deltax*((25./12.)*vif-4.*vif_1+   &
              3.*vif_2-(4./3.)*vif_3+(1./4.)*vif_4)
        dwdxf=jbn(nx,j,k,1)*deltax*((25./12.)*wif-4.*wif_1+   &
              3.*wif_2-(4./3.)*wif_3+(1./4.)*wif_4)
        dpdxf=jbn(nx,j,k,1)*deltax*((25./12.)*pif-4.*pif_1+   &
              3.*pif_2-(4./3.)*pif_3+(1./4.)*pif_4)
        drdxf=jbn(nx,j,k,1)*deltax*((25./12.)*rif-4.*rif_1+   &
              3.*rif_2-(4./3.)*rif_3+(1./4.)*rif_4)

        uci=SQRT(gamma*c3*tii)
        ucf=SQRT(gamma*c3*tif)
        xmachi=mach*uii/sqrt(tii)
        xmachf=mach*uif/sqrt(tif)

!__________________________________________________
!       ENTRADA 
!__________________________________________________
        pl_l1=(uii-uci)*(c3*dpdxi-rii*uci*dudxi)
        pl_l5=pl_l1
        pl_l2=0.5*c2*(pl_l1+pl_l5)
        pl_l3=0.0
        pl_l4=0.0

        pl_d1=(pl_l2+0.5*(pl_l1+pl_l5))/(uci*uci)
        pl_d2=0.0
        pl_d2=0.0
        pl_d3=0.0
        pl_d4=0.0

        frontx(1,1,j,k)=pl_d1
        frontx(1,2,j,k)=0.0         
        frontx(1,3,j,k)=0.0
        frontx(1,4,j,k)=0.0
        frontx(1,5,j,k)=0.0

!__________________________________________________
!       SALIDA 
!__________________________________________________

        p_inf=1.00
        pl_l5=(uif+ucf)*(c3*dpdxf+rif*ucf*dudxf)
        pl_l1=.1*(1.-xmachf*xmachf)*ucf*(pif-p_inf)*c3
        pl_l2=uif*(ucf*ucf*drdxf-c3*dpdxf)
        pl_l3=uif*dvdxf
        pl_l4=uif*dwdxf

        pl_d1=(pl_l2+0.5*(pl_l1+pl_l5))/(ucf*ucf)
        pl_d2=0.5*(pl_l1+pl_l5)
        pl_d3=(pl_l5-pl_l1)/(2.*rif*ucf)
        pl_d4=pl_l3
        pl_d5=pl_l4

        frontx(2,1,j,k)=pl_d1
        frontx(2,2,j,k)=uif*pl_d1+rif*pl_d3
        frontx(2,3,j,k)=vif*pl_d1+rif*pl_d4
        frontx(2,4,j,k)=wif*pl_d1+rif*pl_d5
        frontx(2,5,j,k)=0.5*(uif*uif+vif*vif+wif*wif)*pl_d1+pl_d2/c2+  &
                        rif*uif*pl_d3+rif*vif*pl_d4+rif*wif*pl_d5

       ENDDO
      ENDDO
      RETURN
      END SUBROUTINE FRONTERA_X

!____________________________________________________________________
      SUBROUTINE FXSOURCE (irk,neq)
!____________________________________________________________________

      use dimensiones
      use variables
      use fronterapl
      use tiempo
      implicit none
      integer neq,irk,i,j,k,l,m
      real deltat

       IF(irk.eq.1)THEN
        deltat=dt
       ELSEIF(irk.eq.2)THEN
        deltat=0.25*dt
       ELSEIF(irk.eq.3)THEN
        deltat=(2./3.)*dt
       ENDIF 

       Do k=1,nz
        Do j=1,ny
         um1(1,j,k,neq)=um1(1,j,k,neq)-deltat*frontx(1,neq,j,k)
         um1(nx,j,k,neq)=um1(nx,j,k,neq)-deltat*frontx(2,neq,j,k)
        Enddo
       Enddo

       RETURN
       END SUBROUTINE FXSOURCE

!__________________________________________________
       SUBROUTINE frontera_x_val (irk)
!__________________________________________________

       use dimensiones
       use velocidades
       use variables
       use consadim
       use fronterapl
       use frontera_uval
       use ranaleo
       IMPLICIT NONE
       integer :: i,j,k,l,m,irk
       real :: pe,te,re,eps,bruit

       eps=0.2
       do k=1,nz
       do j=1,ny
          ue=frontin(j,k,2)
          ve=frontin(j,k,3)
          we=frontin(j,k,4)
 
          if((ue.gt.0.159).and.(ue.le.1.097)) then
            ue=ue+eps*(ran(idum)-0.5)
            ve=ve+eps*(ran(idum)-0.5)
            we=we+eps*(ran(idum)-0.5)

          elseif((ue.gt.1.097))then
            ue=ue+eps*(ran(idum)-0.5)*0.2
            ve=ve+eps*(ran(idum)-0.5)*0.2
            we=we+eps*(ran(idum)-0.5)*0.2
          endif
          te=frontin(j,k,1)
          re=um1(1,j,k,1)
          pe=te*re
       
           um1(1,j,k,2)=ue*re
           um1(1,j,k,3)=ve*re
           um1(1,j,k,4)=we*re
           um1(1,j,k,5)=pe/c2+c1*re*(ue**2+ve**2+we**2)
!           u(1,j,k,1)=ue
!           u(1,j,k,2)=ve
!           u(1,j,k,3)=we
!           temp(1,j,k)=te
!           pres(1,j,k)=pe
        enddo
        enddo
        return
        end subroutine frontera_x_val

!__________________________________________________
       SUBROUTINE inifrontera_x_val ()
!__________________________________________________

       use dimensiones
       use velocidades
       use fronterapl
       IMPLICIT NONE
       integer :: i,j,k,l,m
       integer :: in_Nx,in_Ny,in_Nz,in_Nvar

      write(6,*)'LECTURA CAMPO DE ENTRADA X'
      open(11,file='frontx.in',form='unformatted')
      read(11,err=1002)in_Nvar
      if (Nd.ne.in_Nvar) goto 1002
      read(11,err=1003)in_Nx,in_Ny,in_Nz
      if (Nx.ne.in_Nx) goto 1003
      if (Ny.ne.in_Ny) goto 1003
      if (Nz.ne.in_Nz) goto 1003

      read(11,err=1004) frontin
      CLOSE (11)
      
      RETURN

 1002 Print*,in_Nvar,' leida ',Nd,' esperada'
      print*,'ERROR DETECTADO LEYENDO EN NUMERO DE VARIABLES'
      stop
 1003 print*,in_Nx,in_Ny,in_Nz,' leidas ',Nx,Ny,Nz,' eperadas'
      print*,'ERROR DETECTADO LEYENDO EN DIMENSIONES'
      stop
 1004 print*,'ERROR DETECTADO LEYENDO LA MALLA'
      stop

!       do k=1,nz
!       do j=1,ny
!        frontin(j,k,1)=temp(1,j,k)
!        frontin(j,k,2)=u(1,j,k,1)
!        frontin(j,k,3)=u(1,j,k,2)
!        frontin(j,k,4)=u(1,j,k,3)
!        frontin(j,k,5)=pres(1,j,k)
!       enddo
!       enddo
!       return
       end subroutine inifrontera_x_val

!__________________________________________________
       SUBROUTINE inialea ()
!__________________________________________________

       use ifport

       call seed(1995)
 
       return
       end subroutine inialea
