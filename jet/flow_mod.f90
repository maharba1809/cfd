!-----------------------------------------------------------------------
      module indices
!-----------------------------------------------------------------------

      integer i,j,k,l,m

      end  module indices
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
      module derivtools
!-----------------------------------------------------------------------

      real, allocatable, dimension (:) :: du,dv,dw

      end module derivtools

!----------------------------------------------------------------------
      module jacobtools
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: jbn
      real t1,t2,t3,t4,t5,t6,t7,t8,t9

      end module jacobtools

!----------------------------------------------------------------------
      module mallagrid
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: x

      end module mallagrid

!----------------------------------------------------------------------
      module derivvel
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: dvel
      real, allocatable, dimension (:,:,:,:) :: dconc
      real, allocatable, dimension (:,:,:,:) :: dtemp

      end module derivvel

!----------------------------------------------------------------------
      module mflujos
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: e,f,g
      real, allocatable, dimension (:,:,:) :: ev,fv,gv

      end module mflujos

!----------------------------------------------------------------------
      module dmflujos
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: dere,derf,derg
      real, allocatable, dimension (:,:,:) :: derev,derfv,dergv

      end module dmflujos
!----------------------------------------------------------------------
      module viscosidades
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: vis,vis2,vis6,visl,vist

      end module viscosidades

!----------------------------------------------------------------------
      module consadim
!-----------------------------------------------------------------------

      real c1,c2,c3,c4,c5,c6,c7,c8,c9

      end module consadim

!----------------------------------------------------------------------
      module velocidades
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: u
      real, allocatable, dimension (:,:,:) :: conc,temp,pres

      end module velocidades

!----------------------------------------------------------------------
      module variables
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: um,um1,um0

      end module variables

!----------------------------------------------------------------------
      module right
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: rs,rsv

      end module right

!----------------------------------------------------------------------
      module tiempo
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: dx,dy,dz 
      real :: dt,dto,olddt,cflm,cfl,diffdt,up,down,tmax,cvisc
      real :: cfla
      integer :: it,itf,itr

      end module tiempo

!----------------------------------------------------------------------
      module inputname
!-----------------------------------------------------------------------

      character(32) :: inputgrid,inputfield

      end module inputname
!----------------------------------------------------------------------
      module flow
!-----------------------------------------------------------------------

      real :: reynolds,froude,prandtl,prandtl_t,gamma,mach

      end module flow

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
      module sgdmodel
!-----------------------------------------------------------------------

      character(32) :: sgs_model
      real :: amu0,cutoff
      real, allocatable, dimension (:,:,:) :: amut,dxyzsgd
      real, allocatable, dimension (:,:,:) :: dxmsgd,dxpsgd
      real, allocatable, dimension (:,:,:) :: dymsgd,dypsgd
      real, allocatable, dimension (:,:,:) :: dzmsgd,dzpsgd
      integer :: isgm

      end module sgdmodel
!----------------------------------------------------------------------
      module vorticidad
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: wx,wy,wz

      end module vorticidad

!----------------------------------------------------------------------
      module consource
!-----------------------------------------------------------------------

      real :: debit,xydim,velm,deltau,derho,chight,cener

      end module consource

!----------------------------------------------------------------------
      module stat
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: st

      end module stat

!----------------------------------------------------------------------
      module fieldrec
!-----------------------------------------------------------------------

      integer :: n_serie,i_sample,igrava
      real :: dgrava,dgrava1
      character(2) serie
      character(3) sample

      end module fieldrec

!----------------------------------------------------------------------
      module vieu
!-----------------------------------------------------------------------

      real :: tvieu,tvieu1
      integer :: ivieu


      end module vieu
!----------------------------------------------------------------------
      module fronterapl
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:,:) :: frontx,frontz,fronty
      real, allocatable, dimension (:,:,:) :: frontin
      integer, allocatable, dimension (:,:,:) :: imask

      end module fronterapl
!----------------------------------------------------------------------
      module frontera_uval
!-----------------------------------------------------------------------
      real :: ue,ve,we

     end module frontera_uval

!----------------------------------------------------------------------
      module bob
!-----------------------------------------------------------------------

      real, allocatable, dimension (:,:,:) :: esplayer
      real, allocatable, dimension (:) :: esource
      integer :: ief,iei,iefi,ieii

      end module bob
!----------------------------------------------------------------------
      module ranaleo
!-----------------------------------------------------------------------

      integer :: idum

      end module ranaleo



