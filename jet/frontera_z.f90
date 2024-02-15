!________________________________________________
       SUBROUTINE frontera_z ()
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
      real dudzi,drdzi,dvdzi,dpdzi,dwdzi
      real dudzf,drdzf,dvdzf,dpdzf,dwdzf
      real pl_l1,pl_l2,pl_l3,pl_l4,pl_l5
      real pl_d1,pl_d2,pl_d3,pl_d4,pl_d5
      real p_inf,xmachi,xmachf

!______________________________________________
!    z=1
!______________________________________________

      DO i=1,nx
       DO j=1,ny
        uii=u(i,j,1,1)
        vii=u(i,j,1,2)
        wii=u(i,j,1,3)
        pii=pres(i,j,1)
        rii=um(i,j,1,1)
        tii=pii/rii

        uii_1=u(i,j,2,1)
        vii_1=u(i,j,2,2)
        wii_1=u(i,j,2,3)
        pii_1=pres(i,j,2)
        rii_1=um(i,j,2,1)

        uii_2=u(i,j,3,1)
        vii_2=u(i,j,3,2)
        wii_2=u(i,j,3,3)
        pii_2=pres(i,j,3)
        rii_2=um(i,j,3,1)

        uii_3=u(i,j,4,1)
        vii_3=u(i,j,4,2)
        wii_3=u(i,j,4,3)
        pii_3=pres(i,j,4)
        rii_3=um(i,j,4,1)

        uii_4=u(i,j,5,1)
        vii_4=u(i,j,5,2)
        wii_4=u(i,j,5,3)
        pii_4=pres(i,j,5)
        rii_4=um(i,j,5,1)

        dudzi=jbn(i,j,1,9)*deltaz*((-25./12.)*uii+4.*uii_1-   &
              3.*uii_2+(4./3.)*uii_3-(1./4.)*uii_4)
        dvdzi=jbn(i,j,1,9)*deltaz*((-25./12.)*vii+4.*vii_1-   &
              3.*vii_2+(4./3.)*vii_3-(1./4.)*vii_4)
        dwdzi=jbn(i,j,1,9)*deltaz*((-25./12.)*wii+4.*wii_1-   &
              3.*wii_2+(4./3.)*wii_3-(1./4.)*wii_4)
        dpdzi=jbn(i,j,1,9)*deltaz*((-25./12.)*pii+4.*pii_1-   &
              3.*pii_2+(4./3.)*pii_3-(1./4.)*pii_4)
        drdzi=jbn(i,j,1,9)*deltaz*((-25./12.)*rii+4.*rii_1-   &
              3.*rii_2+(4./3.)*rii_3-(1./4.)*rii_4)
!        write(6,*)'i',dudzi,dvdzi,dwdzi,dpdzi,drdzi

        uif=u(i,j,nz,1)
        vif=u(i,j,nz,2)
        wif=u(i,j,nz,3)
        pif=pres(i,j,nz)
        rif=um(i,j,nz,1)
        tif=pif/rif

        uif_1=u(i,j,nz-1,1)
        vif_1=u(i,j,nz-1,2)
        wif_1=u(i,j,nz-1,3)
        pif_1=pres(i,j,nz-1)
        rif_1=um(i,j,nz-1,1)

        uif_2=u(i,j,nz-2,1)
        vif_2=u(i,j,nz-2,2)
        wif_2=u(i,j,nz-2,3)
        pif_2=pres(i,j,nz-2)
        rif_2=um(i,j,nz-2,1)

        uif_3=u(i,j,nz-3,1)
        vif_3=u(i,j,nz-3,2)
        wif_3=u(i,j,nz-3,3)
        pif_3=pres(i,j,nz-3)
        rif_3=um(i,j,nz-3,1)

        uif_4=u(i,j,nz-4,1)
        vif_4=u(i,j,nz-4,2)
        wif_4=u(i,j,nz-4,3)
        pif_4=pres(i,j,nz-4)
        rif_4=um(i,j,nz-4,1)

        dudzf=jbn(i,j,nz,9)*deltaz*((25./12.)*uif-4.*uif_1+   &
              3.*uif_2-(4./3.)*uif_3+(1./4.)*uif_4)
        dvdzf=jbn(i,j,nz,9)*deltaz*((25./12.)*vif-4.*vif_1+   &
              3.*vif_2-(4./3.)*vif_3+(1./4.)*vif_4)
        dwdzf=jbn(i,j,nz,9)*deltaz*((25./12.)*wif-4.*wif_1+   &
              3.*wif_2-(4./3.)*wif_3+(1./4.)*wif_4)
        dpdzf=jbn(i,j,nz,9)*deltaz*((25./12.)*pif-4.*pif_1+   &
              3.*pif_2-(4./3.)*pif_3+(1./4.)*pif_4)
        drdzf=jbn(i,j,nz,9)*deltaz*((25./12.)*rif-4.*rif_1+   &
              3.*rif_2-(4./3.)*rif_3+(1./4.)*rif_4)

!        write(6,*)'f',dudzf,dvdzf,dwdzf,dpdzf,drdzf
        uci=SQRT(gamma*c3*tii)
        ucf=SQRT(gamma*c3*tif)
        xmachi=mach*wii/sqrt(tii)
        xmachf=mach*wif/sqrt(tif)
!__________________________________________________
!       PARED CON DESLIZAMIENTO (w=0)
!__________________________________________________
         pl_l1=(wii-uci)*(c3*dpdzi-rii*uci*dwdzi)
         pl_l5=pl_l1
         pl_l2=0.0
         pl_l3=0.0
         pl_l4=0.0

         pl_d1=(pl_l2+0.5*(pl_l1+pl_l5))/(uci*uci)
         pl_d2=0.5*(pl_l1+pl_l2)
         pl_d3=(pl_l5-pl_l1)/(2.*rii*uci)
         pl_d4=0.
         pl_d5=0.

!         write(6,*)pl_l1,pl_d1,pl_d2,pl_d3

         frontz(1,1,i,j)=pl_d1
         frontz(1,2,i,j)=uif*pl_d1
         frontz(1,3,i,j)=vif*pl_d1
         frontz(1,4,i,j)=0.
         frontz(1,5,i,j)=0.5*(uii*uii+vii*vii)*pl_d1+pl_d2/c2  

         pl_l5=(wif+ucf)*(c3*dpdzf+rif*ucf*dwdzf)
         pl_l1=pl_l5
         pl_l2=0.0
         pl_l3=0.0
         pl_l4=0.0

         pl_d1=(pl_l2+0.5*(pl_l1+pl_l5))/(ucf*ucf)
         pl_d2=0.5*(pl_l1+pl_l2)
         pl_d3=(pl_l5-pl_l1)/(2.*rif*ucf)
         pl_d4=0.
         pl_d5=0.

         frontz(2,1,i,j)=pl_d1
         frontz(2,2,i,j)=uif*pl_d1
         frontz(2,3,i,j)=vif*pl_d1
         frontz(2,4,i,j)=0.
         frontz(2,5,i,j)=0.5*(uii*uii+vii*vii)*pl_d1+pl_d2/c2

       ENDDO
      ENDDO
      RETURN
      END SUBROUTINE FRONTERA_Z

!____________________________________________________________________
      SUBROUTINE FZSOURCE (irk,neq)
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

       Do j=1,ny
        Do i=2,nx
         um1(i,j,1,neq)=um1(i,j,1,neq)-deltat*frontz(1,neq,i,j)
         um1(i,j,nz,neq)=um1(i,j,nz,neq)-deltat*frontz(2,neq,i,j)
        Enddo
       Enddo

       RETURN
       END SUBROUTINE FZSOURCE

!__________________________________________________
       SUBROUTINE frontera_z_val ()
!__________________________________________________

       use dimensiones
       use velocidades
       use variables
       use consadim

       IMPLICIT NONE
       integer :: i,j,k,l,m
       real :: ue,ve,we,re,pe

       DO j=1,ny
       DO i=1,nx
        um1(i,j,1,4)=0.0
        um1(i,j,nz,4)=0.0
        re=um1(i,j,1,1)
        ue=um1(i,j,1,2)/um1(i,j,1,1)
        ve=um1(i,j,1,3)/um1(i,j,1,1)
        we=um1(i,j,1,4)/um1(i,j,1,1)
        pe=1.0*um1(i,j,1,1)
        um1(i,j,1,5)=pe/c2+c1*re*(ue**2+ve**2)
        re=um1(i,j,nz,1)
        ue=um1(i,j,nz,2)/um1(i,j,nz,1)
        ve=um1(i,j,nz,3)/um1(i,j,nz,1)
        we=um1(i,j,nz,4)/um1(i,j,nz,1)
        pe=1.0*um1(i,j,nz,1)
        um1(i,j,nz,5)=pe/c2+c1*re*(ue**2+ve**2)
       ENDDO
       ENDDO
       RETURN
       END SUBROUTINE frontera_z_val
