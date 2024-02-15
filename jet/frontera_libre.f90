!________________________________________________
       SUBROUTINE frontera_libre ()
!________________________________________________

      use dimensiones
      use velocidades
      use variables
      use deltas
      use consadim
      use jacobtools

      IMPLICIT NONE
      integer :: i,j,k,l,m,irk
      real ri,pi,ui,vi,wi
      real rt,pt,ut,vt,wt
      real re,pe,ue,ve,we
      real dlz

!______________________________________________
!    z=1
!______________________________________________
      dlz=deltaz

      DO j=1,ny
       DO i=1,nx
        ut=u(i,j,1,1)
        vt=u(i,j,1,2)
        wt=u(i,j,1,3)
        pt=pres(i,j,1)
        rt=um1(i,j,1,1)

        ui=u(i,j,2,1)
        vi=u(i,j,2,2)
        wi=u(i,j,2,3)
        pi=pres(i,j,2)
        ri=um1(i,j,2,1)

        ue=ut
        ve=vt
        we=wt
        pe=pt
        re=rt

!__________________________________________________
!       FRONTERA LBRE THOMPSON O(1)
!__________________________________________________
        CALL libre(rt,ut,vt,wt,pt,ri,ui,vi,wi,pi,re,ue,ve,we,pe, &
                       jbn(i,j,1,7),jbn(i,j,1,9),dlz,-1)

         um1(i,j,1,1)=rt
         um1(i,j,1,2)=rt*ut
         um1(i,j,1,3)=rt*vt
         um1(i,j,1,4)=rt*wt
         um1(i,j,1,5)=pt/c2+c1*rt*(ut**2+vt**2+wt**2)
!         pres(i,j,1)=pt
!         temp(i,j,1)=pt/rt
!         u(i,j,1,1)=ut
!         u(i,j,1,2)=vt
!         u(i,j,1,3)=wt
       ENDDO
      ENDDO
!______________________________________________________________________
!     z=nz
!______________________________________________________________________
      DO j=1,ny
       DO i=1,nx
        ut=u(i,j,nz,1)
        vt=u(i,j,nz,2)
        wt=u(i,j,nz,3)
        pt=pres(i,j,nz)
        rt=um1(i,j,nz,1)

        ui=u(i,j,nz-1,1)
        vi=u(i,j,nz-1,2)
        wi=u(i,j,nz-1,3)
        pi=pres(i,j,nz-1)
        ri=um1(i,j,nz-1,1)

        ue=ut
        ve=vt
        we=wt
        pe=pt
        re=rt

!__________________________________________________
!       FRONTERA LBRE THOMPSON O(1)
!__________________________________________________

!         write(6,*)'antes',wt,rt,pt
        CALL libre(rt,ut,vt,wt,pt,ri,ui,vi,wi,pi,re,ue,ve,we,pe, &
                       jbn(i,j,nz,7),jbn(i,j,nz,9),dlz,1)

         um1(i,j,nz,1)=rt
         um1(i,j,nz,2)=rt*ut
         um1(i,j,nz,3)=rt*vt
         um1(i,j,nz,4)=rt*wt
         um1(i,j,nz,5)=pt/c2+c1*rt*(ut**2+vt**2+wt**2)
!         pres(i,j,nz)=pt
!         temp(i,j,nz)=pt/rt
!         u(i,j,nz,1)=ut
!         u(i,j,nz,2)=vt
!         u(i,j,nz,3)=wt
!         write(6,*)'despues',wt,rt,pt
       ENDDO
      ENDDO
      RETURN
      END subroutine frontera_libre

!__________________________________________________________________
      SUBROUTINE libre(rt,ut,vt,wt,pt,ri,ui,vi,wi,pi,     &
                      re,ue,ve,we,pe,t1,t2,dl,isens)
!__________________________________________________________________
      use consadim
      use tiempo
      use flow
      IMPLICIT NONE 
!
! Conditions de sorties libres avec les caracteristiques misent a zero
! si les vitesses caracteristiquent sont negatives
!   isens donne le sens de sortie : 1 si positif (est, Nord)
!                                  -1 si negatif (ouest, sud)
! Formulation en ro, rou, rov, p

      integer isens,irk
      real rt,ut,vt,wt,pt,ri,ui,vi,wi,pi
      real re,ue,ve,we,pe,t1,t2,dl
      real dal, TJC,TJCQ,alpha,ct,cti,aa
      real wpt,wot,wptpa,wptma,wpi,woi,wpe,woe
      real w1t,w1i,w1e,w2t,w2i,w2e,w3t,w3i,w3e
      real w4t,w4i,w4e,w5t,w5i,w5e
      real uu,ww,vv,pp,rr,wo,wp
      real wptmad,wptpad,wptd

      dal=dt/dl
      TJC=t1**2+t2**2
      TJCQ=SQRT(TJC)
      alpha=SQRT(c3)
      ct=SQRT(gamma*pt/rt)
      cti=1.0/ct
      aa=alpha*ct

      wpt=t1*ut+t2*wt
      wot=t2*ut-t1*wt
      wptpa=wpt+aa*TJCQ
      wptma=wpt-aa*TJCQ
      wpi=t1*ui+t2*wi
      woi=t2*ui-t1*wi
      wpe=t1*ue+t2*we
      woe=t2*ue-t1*we

      w1t=-ct**2*rt+pt
      w1i=-ct**2*ri+pi
      w1e=-ct**2*re+pe

      w2t=rt*wot
      w2i=rt*woi
      w2e=rt*woe

      w3t=rt*vt
      w3i=rt*vi
      w3e=rt*ve

      w4t=ct*rt*wpt+pt*TJCQ*alpha
      w4i=ct*rt*wpi+pi*TJCQ*alpha
      w4e=ct*rt*wpe+pe*TJCQ*alpha

      w5t=-ct*rt*wpt+pt*TJCQ*alpha
      w5i=-ct*rt*wpi+pi*TJCQ*alpha
      w5e=-ct*rt*wpe+pe*TJCQ*alpha



!   Vitesse caracteristique et valeur ABSolue de ces vitesse sur-indicee "d"
!    Libre : si vitesse caracteristique negative on ne ramene rien (Thompson)
!            autre possibilite : voir Poinsot

      IF (wpt*isens.LT.0.) wpt=0.0
      IF (wptpa*isens.LT.0.) wptpa=0.0
      IF (wptma*isens.LT.0.) wptma=0.0
      wptd=ABS(wpt)
      wptpad=ABS(wptpa)
      wptmad=ABS(wptma)

! Schema implicite en temps

      w1t=1.0/(1.0+wptd*dal)*(w1t+      &
          (wptd-wpt)*0.50*dal*w1e+(wpt+wptd)*0.5*dal*w1i)
      w2t=1.0/(1.e0+wptd*dal)*(w2t+     &
          (wptd-wpt)*0.5*dal*w2e+(wpt+wptd)*0.5*dal*w2i)
      w3t=1.0/(1.e0+wptd*dal)*(w3t+     &
          (wptd-wpt)*0.5*dal*w3e+(wpt+wptd)*0.5*dal*w3i)
      w4t=1.0/(1.0+wptpad*dal)*(w4t+     &
          (wptpad-wptpa)*0.5*dal*w4e+(wptpa+wptpad)*0.5*dal*w4i)
      w5t=1.0/(1.0+wptmad*dal)*(w5t+     &
          (wptmad-wptma)*0.5*dal*w5e+(wptma+wptmad)*0.5*dal*w5i)

! recombinaison en variable r, u, v, p

      pp=(w4t+w5t)/(2.0*alpha*TJCQ)
      rr=(pp-w1t)*(cti**2)
      wo=w2t/rt
      vt=w3t/rt
      wp=(w4t-w5t)/(2.*rt*ct)
      uu=(wp*t1+wo*t2)/TJC
      ww=(wp*t2-wo*t1)/TJC
      rt=rr
      ut=uu
      wt=ww
      pt=pp

      RETURN
      END SUBROUTINE LIBRE

