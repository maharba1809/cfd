!-----------------------------------------------------------------------
      module dimensiones
!-----------------------------------------------------------------------

      integer nx,ny,nz,nd,ne

      end  module dimensiones

!-----------------------------------------------------------------------
      module deltas
!-----------------------------------------------------------------------

      real :: deltax,deltay,deltaz

      end  module deltas

!----------------------------------------------------------------------
      module consderper
!-----------------------------------------------------------------------

      real, allocatable, dimension (:) :: axp,ayp,azp
      real, allocatable, dimension (:) :: bxp,byp,bzp
      real, allocatable, dimension (:) :: cxp,cyp,czp

      end  module consderper
!----------------------------------------------------------------------
      module consdernper
!-----------------------------------------------------------------------

      real, allocatable, dimension (:) :: axnp,aynp,aznp
      real, allocatable, dimension (:) :: bxnp,bynp,bznp
      real, allocatable, dimension (:) :: cxnp,cynp,cznp

      end  module consdernper

!----------------------------------------------------------------------
      module consrs
!-----------------------------------------------------------------------

      real ars,brs,ars1,brs1,crs1,ars2

      end module consrs


!----------------------------------------------------------------------
      module consderper_f
!-----------------------------------------------------------------------

      real, allocatable, dimension (:) :: axpf,aypf,azpf
      real, allocatable, dimension (:) :: bxpf,bypf,bzpf
      real, allocatable, dimension (:) :: cxpf,cypf,czpf

      end  module consderper_f
!----------------------------------------------------------------------
      module consdernper_f
!-----------------------------------------------------------------------

      real, allocatable, dimension (:) :: axnpf,aynpf,aznpf
      real, allocatable, dimension (:) :: bxnpf,bynpf,bznpf
      real, allocatable, dimension (:) :: cxnpf,cynpf,cznpf

      end  module consdernper_f

!----------------------------------------------------------------------
      module consrs_f
!-----------------------------------------------------------------------

      real :: arsf,brsf,crsf,drsf,ersf,frsf
      real :: arsf1,brsf1,crsf1,drsf1,ersf1
      real :: arsf2,brsf2,crsf2,drsf2
      real :: arsf3,brsf3,crsf3,drsf3,ersf3,frsf3,grsf3
      real :: arsf4,brsf4,crsf4,drsf4,ersf4,frsf4,grsf4

      end module consrs_f

!----------------------------------------------------------------------
      module consfiltro
!-----------------------------------------------------------------------

      real :: cfilt

      end module consfiltro

!----------------------------------------------------------------------
      module funciones
!-----------------------------------------------------------------------

      real, allocatable, dimension (:) :: u,u2

      end module funciones
