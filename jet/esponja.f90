!__________________________________________________________________
         SUBROUTINE esponja()
!__________________________________________________________________

      use dimensiones
      use variables
      use bob
      IMPLICIT NONE
      INTEGER i,j,k,m,l,nzt,nyt
      

      nzt=int(float(nz)/2)
      nyt=int(float(ny)/2)

      do m=1,nd
      do k=1,nz
      do j=1,ny
       esplayer(j,k,m)=0.0
      enddo
      enddo
      enddo
!
      do m=1,nd
      do k=1,nz
      do j=1,ny
      do i=iei,ief
       esplayer(j,k,m)=esplayer(j,k,m)+um(i,j,k,m)
      enddo
      enddo
      enddo
      enddo
!
      do m=1,nd
      do j=1,ny
      do k=1,nzt
        l=nz-(k-1)
       IF(m.eq.4) then
        esplayer(j,k,m)=(esplayer(j,k,m)-esplayer(j,l,m))/2.
       ELSE
        esplayer(j,k,m)=(esplayer(j,k,m)+esplayer(j,l,m))/2.
       ENDIF
      enddo
      enddo
      enddo

      do m=1,nd
      do j=1,ny
      do k=1,nzt
        l=nz-(k-1)
       IF(m.eq.4) then
        esplayer(j,l,m)=-esplayer(j,k,m)
       ELSE
        esplayer(j,l,m)=esplayer(j,k,m)
       ENDIF
      enddo
      enddo
      enddo

      do m=1,nd
      do k=1,nz
      do j=1,nyt
        l=ny-(j-1)
       IF(m.eq.3) then
        esplayer(j,k,m)=(esplayer(j,k,m)-esplayer(l,k,m))/2.
       ELSE
        esplayer(j,k,m)=(esplayer(j,k,m)+esplayer(l,k,m))/2.
       ENDIF
      enddo
      enddo
      enddo

      do m=1,nd
      do k=1,nz
      do j=1,nyt
        l=ny-(j-1)
       IF(m.eq.3) then
        esplayer(j,l,m)=-esplayer(j,k,m)
       ELSE
        esplayer(j,l,m)=esplayer(j,k,m)
       ENDIF
      enddo
      enddo
      enddo
!

      do m=1,nd
      do j=1,ny
      do k=1,nz
       esplayer(j,k,m)=esplayer(j,k,m)/(float(ief-iei+1))
      enddo
      enddo
      enddo

      RETURN
      END SUBROUTINE ESPONJA


!________________________________________________________________
      SUBROUTINE INI_ESPONJA()
!________________________________________________________________

      use dimensiones
      use mallagrid
      use bob
      IMPLICIT NONE
      real axa
      integer i,j,k,l,m

      iei=90
      ief=120
      ieii=121
      iefi=150
      do i=1,nx
       if (i.lt.ieii) then
        esource(i)=0.0
       else
        axa=(x(iefi,1,1,1)-x(i,1,1,1))/(x(iefi,1,1,1)-x(ieii,1,1,1))
        esource(i)=0.8*exp(-7.0*(axa)**4.0)    
       endif
      enddo
       return
       end subroutine ini_esponja
