!   _______________________________________________________
      SUBROUTINE INI_SOURCE()
!   _______________________________________________________

      use dimensiones
      use variables
      use consource

      debit=4./3.
      xydim=float(nx*ny)
      velm=0.0
      deltau=0.0
      derho=0.0
      chight=2.0
      cener=c6
      RETURN
      END SUBROUTINE INI_SOURCE
!   _______________________________________________________
      SUBROUTINE SOURCE(neq)
!   _______________________________________________________

      use dimensiones
      use variables
      use consource
      use mallagrid

      IMPLICIT NONE
      integer :: neq,i,j,k,l,m
      real, dimension(nz) :: sum,sut,debiz,derhoz
      real :: debitu

        IF (neq.eq.2) THEN
        DO k=1,nz
         sum(k)=0.0
         sut(k)=0.0
        ENDDO
        DO k=1,nz
         DO j=1,ny
          DO i=1,nx
           sum(k)=sum(k)+um1(i,j,k,2)
           sut(k)=sut(k)+um1(i,j,k,1)
          ENDDO
         ENDDO
        ENDDO
        DO k=1,nz
         sum(k)=sum(k)/xydim
         sut(k)=sut(k)/xydim
        ENDDO

        DO k=2,nz
         debiz(k)=(sum(k)+sum(k-1))*0.5e0
         derhoz(k)=(sut(k)+sut(k-1))*0.5e0
        ENDDO
        debiz(1)=sum(1)
        derhoz(1)=sut(1)

        debitu=0.e0
        derho=0.e0
        DO k=2,nz
         debitu=debitu+debiz(k)*(x(1,1,k,3)-x(1,1,k-1,3))
         derho=derho+derhoz(k)*(x(1,1,k,3)-x(1,1,k-1,3))
        ENDDO
        debitu=debitu+debiz(1)*(x(1,1,2,3)-x(1,1,1,3))
        derho=derho+derhoz(1)*(x(1,1,2,3)-x(1,1,1,3))
!
        velm=debitu/derho
        deltau=(debit-debitu)/chight
        derho=derho/chight
        PRINT *,'Correcion del gasto',debitu,'  ','du=',deltau
!        DO k=2,nz-1
!         DO j=1,ny
!          DO i=1,nx
!           um1(i,j,k,2)=um1(i,j,k,2)+deltau
!          ENDDO
!         ENDDO
!        ENDDO
       ELSE IF (neq.eq.5) THEN
!        DO k=2,nz-1
!         DO j=1,ny
!          DO i=1,nx
!           um1(i,j,k,5)=um1(i,j,k,5)+cener*velm*deltau
!          ENDDO
!         ENDDO
!        ENDDO
       ENDIF
       RETURN
       END SUBROUTINE SOURCE
