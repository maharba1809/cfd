!________________________________________________
       SUBROUTINE frontera_y ()
!________________________________________________
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
      use derivvel

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
      real dudyi,drdyi,dvdyi,dpdyi,dwdyi
      real dudyf,drdyf,dvdyf,dpdyf,dwdyf
      real pl_l1,pl_l2,pl_l3,pl_l4,pl_l5
      real pl_d1,pl_d2,pl_d3,pl_d4,pl_d5
      real p_inf,xmachi,xmachf

!______________________________________________
!    z=1
!______________________________________________

      DO i=1,nx
       DO k=1,nz
        uii=u(i,1,k,1)
        vii=u(i,1,k,2)
        wii=u(i,1,k,3)
        pii=pres(i,1,k)
        rii=um(i,1,k,1)
        tii=pii/rii

        uii_1=u(i,2,k,1)
        vii_1=u(i,2,k,2)
        wii_1=u(i,2,k,3)
        pii_1=pres(i,2,k)
        rii_1=um(i,2,k,1)

        uii_2=u(i,3,k,1)
        vii_2=u(i,3,k,2)
        wii_2=u(i,3,k,3)
        pii_2=pres(i,3,k)
        rii_2=um(i,3,k,1)

        uii_3=u(i,4,k,1)
        vii_3=u(i,4,k,2)
        wii_3=u(i,4,k,3)
        pii_3=pres(i,4,k)
        rii_3=um(i,4,k,1)

        uii_4=u(i,5,k,1)
        vii_4=u(i,5,k,2)
        wii_4=u(i,5,k,3)
        pii_4=pres(i,5,k)
        rii_4=um(i,5,k,1)

        dudyi=jbn(i,1,k,5)*deltay*((-25./12.)*uii+4.*uii_1-   &
              3.*uii_2+(4./3.)*uii_3-(1./4.)*uii_4)
        dvdyi=jbn(i,1,k,5)*deltay*((-25./12.)*vii+4.*vii_1-   &
              3.*vii_2+(4./3.)*vii_3-(1./4.)*vii_4)
        dwdyi=jbn(i,1,k,5)*deltay*((-25./12.)*wii+4.*wii_1-   &
              3.*wii_2+(4./3.)*wii_3-(1./4.)*wii_4)
        dpdyi=jbn(i,1,k,5)*deltay*((-25./12.)*pii+4.*pii_1-   &
              3.*pii_2+(4./3.)*pii_3-(1./4.)*pii_4)
        drdyi=jbn(i,1,k,5)*deltay*((-25./12.)*rii+4.*rii_1-   &
              3.*rii_2+(4./3.)*rii_3-(1./4.)*rii_4)
!        write(6,*)'i',dudyi,dvdyi,dwdyi,dpdyi,drdyi

        uif=u(i,ny,k,1)
        vif=u(i,ny,k,2)
        wif=u(i,ny,k,3)
        pif=pres(i,ny,k)
        rif=um(i,ny,k,1)
        tif=pif/rif

        uif_1=u(i,ny-1,k,1)
        vif_1=u(i,ny-1,k,2)
        wif_1=u(i,ny-1,k,3)
        pif_1=pres(i,ny-1,k)
        rif_1=um(i,ny-1,k,1)

        uif_2=u(i,ny-2,k,1)
        vif_2=u(i,ny-2,k,2)
        wif_2=u(i,ny-2,k,3)
        pif_2=pres(i,ny-2,k)
        rif_2=um(i,ny-2,k,1)

        uif_3=u(i,ny-3,k,1)
        vif_3=u(i,ny-3,k,2)
        wif_3=u(i,ny-3,k,3)
        pif_3=pres(i,ny-3,k)
        rif_3=um(i,ny-3,k,1)

        uif_4=u(i,ny-4,k,1)
        vif_4=u(i,ny-4,k,2)
        wif_4=u(i,ny-4,k,3)
        pif_4=pres(i,ny-4,k)
        rif_4=um(i,ny-4,k,1)

        dudyf=jbn(i,ny,k,5)*deltay*((25./12.)*uif-4.*uif_1+   &
              3.*uif_2-(4./3.)*uif_3+(1./4.)*uif_4)
        dvdyf=jbn(i,ny,k,5)*deltay*((25./12.)*vif-4.*vif_1+   &
              3.*vif_2-(4./3.)*vif_3+(1./4.)*vif_4)
        dwdyf=jbn(i,ny,k,5)*deltay*((25./12.)*wif-4.*wif_1+   &
              3.*wif_2-(4./3.)*wif_3+(1./4.)*wif_4)
        dpdyf=jbn(i,ny,k,5)*deltay*((25./12.)*pif-4.*pif_1+   &
              3.*pif_2-(4./3.)*pif_3+(1./4.)*pif_4)
        drdyf=jbn(i,ny,k,5)*deltay*((25./12.)*rif-4.*rif_1+   &
              3.*rif_2-(4./3.)*rif_3+(1./4.)*rif_4)

!        write(6,*)'f',dudyf,dvdyf,dwdyf,dpdyf,drdyf
        uci=SQRT(gamma*c3*tii)
        ucf=SQRT(gamma*c3*tif)
        xmachi=mach*vii/sqrt(tii)
        xmachf=mach*vif/sqrt(tif)
!__________________________________________________
!       PARED CON DESLIZAMIENTO (w=0)
!__________________________________________________
         pl_l1=(vii-uci)*(c3*dpdyi-rii*uci*dvdyi)
         pl_l5=pl_l1
         pl_l2=0.0
         pl_l3=0.0
         pl_l4=0.0

         pl_d1=(pl_l2+0.5*(pl_l1+pl_l5))/(uci*uci)
         pl_d2=0.5*(pl_l1+pl_l2)
         pl_d3=(pl_l5-pl_l1)/(2.*rii*uci)
         pl_d4=0.
         pl_d5=0.

!         write(6,*)'i',pl_l1,pl_d1,pl_d2,pl_d3

         fronty(1,1,i,k)=pl_d1
         fronty(1,2,i,k)=uif*pl_d1
         fronty(1,3,i,k)=0.0
         fronty(1,4,i,k)=wif*pl_d1
         fronty(1,5,i,k)=0.5*(uii*uii+wii*wii)*pl_d1+pl_d2/c2  

         pl_l5=(vif+ucf)*(c3*dpdyf+rif*ucf*dvdyf)
         pl_l1=pl_l5
         pl_l2=0.0
         pl_l3=0.0
         pl_l4=0.0

         pl_d1=(pl_l2+0.5*(pl_l1+pl_l5))/(ucf*ucf)
         pl_d2=0.5*(pl_l1+pl_l2)
         pl_d3=(pl_l5-pl_l1)/(2.*rif*ucf)
         pl_d4=0.
         pl_d5=0.
!         write(6,*)'f',pl_l1,pl_d1,pl_d2,pl_d3

         fronty(2,1,i,k)=pl_d1
         fronty(2,2,i,k)=uif*pl_d1
         fronty(2,3,i,k)=0.0
         fronty(2,4,i,k)=wif*pl_d1
         fronty(2,5,i,k)=0.5*(uii*uii+wii*wii)*pl_d1+pl_d2/c2

       ENDDO
      ENDDO
      RETURN
      END SUBROUTINE FRONTERA_Y

!____________________________________________________________________
      SUBROUTINE FYSOURCE (irk,neq)
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

       Do k=2,nz-1
        Do i=2,nx
         um1(i,1,k,neq)=um1(i,1,k,neq)-deltat*fronty(1,neq,i,k)
         um1(i,ny,k,neq)=um1(i,ny,k,neq)-deltat*fronty(2,neq,i,k)
        Enddo
       Enddo

       RETURN
       END SUBROUTINE FYSOURCE

!__________________________________________________
       SUBROUTINE frontera_y_val ()
!__________________________________________________

       use dimensiones
       use velocidades
       use variables
       use consadim

       IMPLICIT NONE
       integer :: i,j,k,l,m
       real :: ue,ve,we,re,pe

       DO k=1,nz
       DO i=1,nx
        um1(i,1,k,3)=0.0
        um1(i,ny,k,3)=0.0
        re=um1(i,1,k,1)
        ue=um1(i,1,k,2)/um1(i,1,k,1)
        ve=um1(i,1,k,3)/um1(i,1,k,1)
        we=um1(i,1,k,4)/um1(i,1,k,1)
        pe=1.0*um1(i,1,k,1)
        um1(i,1,k,5)=pe/c2+c1*re*(ue**2.+ve**2.+we**2.)
        re=um1(i,ny,k,1)
        ue=um1(i,ny,k,2)/um1(i,ny,k,1)
        ve=um1(i,ny,k,3)/um1(i,ny,k,1)
        we=um1(i,ny,k,4)/um1(i,ny,k,1)
        pe=1.0*um1(i,ny,k,1)
        um1(i,ny,k,5)=pe/c2+c1*re*(ue**2.+ve**2.+we**2.)
       ENDDO
       ENDDO
       RETURN
       END SUBROUTINE frontera_y_val
