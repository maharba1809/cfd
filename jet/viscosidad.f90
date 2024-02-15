      SUBROUTINE viscosidad()
      use dimensiones
      use velocidades
      use viscosidades
      use sgdmodel
      use flow
      use consadim
      IMPLICIT NONE
      integer :: i,j,k,l,m

      DO k=1,nz
       DO j=1,ny
        DO i=1,nx
          vis(i,j,k)=1.0
          vis6(i,j,k)=(c5*vis(i,j,k)+              &
                      gamma*amut(i,j,k)/prandtl_t)/c2
!          visl(i,j,k)=(c4*vis(i,j,k))
!          vist(i,j,k)=(gamma*amut(i,j,k)/prandtl_t)/c2
          vis2(i,j,k)=c4*c6*vis(i,j,k)
          vis(i,j,k)=c4*vis(i,j,k)+amut(i,j,k)
        END DO
       END DO
      END DO
      RETURN
      END SUBROUTINE VISCOSIDAD
