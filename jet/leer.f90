      SUBROUTINE leer_malla()
      use dimensiones
      use mallagrid
      use inputname
      IMPLICIT NONE
      integer :: i,j,k,l,m
      integer :: in_Nvar,in_Nx,in_Ny,in_Nz
      real, allocatable, dimension (:,:,:) :: buffer


      write(6,*)'LECTURA MALLA'
      open(11,file=inputgrid,form='unformatted')
      read(11,err=1002)in_Nvar
      write(*,*)'leido',in_nvar,'numero dimensiones',nd
      
      if (Nd.ne.in_Nvar) goto 1002
      read(11,err=1003)in_Nx,in_Ny,in_Nz
      if (Nx.ne.in_Nx) goto 1003
      if (Ny.ne.in_Ny) goto 1003
      if (Nz.ne.in_Nz) goto 1003

      allocate(buffer(nx,ny,nz))
      write(*,*)nx,ny,nz
      
      read(11,err=1004) buffer

      call scatter3d(buffer,x(:,:,:,1))
      read(11,err=1005) buffer
      call scatter3d(buffer,x(:,:,:,2))
      read(11,err=1006) buffer
      call scatter3d(buffer,x(:,:,:,3))

      deallocate(buffer)
      CLOSE (11)
      RETURN

 1002 Print*,in_Nvar,' leida ',Nd,' esperada'
      print*,'ERROR DETECTADO LEYENDO EN NUMERO DE VARIABLES'
      stop
 1003 print*,in_Nx,in_Ny,in_Nz,' leidas ',Nx,Ny,Nz,' esperadas'
      print*,'ERROR DETECTADO LEYENDO EN DIMENSIONES'
      stop
 1004 print*,'ERROR DETECTADO LEYENDO LA MALLA x',nx,in_Nx
      stop
 1005 print*,'ERROR DETECTADO LEYENDO LA MALLA y'
      stop
 1006 print*,'ERROR DETECTADO LEYENDO LA MALLA z'
      stop
      END SUBROUTINE LEER_MALLA







      SUBROUTINE leer_campos()
      use dimensiones
      use velocidades
      use variables
      use consadim
      use inputname
      IMPLICIT NONE
      integer :: i,j,k,l,m
      integer :: in_Nvar,in_Nx,in_Ny,in_Nz
      real, allocatable, dimension (:,:,:) :: buffer

      write(6,*)'LECTURA CAMPOS'
      open(11,file=inputfield,form='unformatted')
      read(11,err=1002)in_Nvar
      if (Nd.ne.in_Nvar) goto 1002
      read(11,err=1003)in_Nx,in_Ny,in_Nz
      if (Nx.ne.in_Nx) goto 1003
      if (Ny.ne.in_Ny) goto 1003
      if (Nz.ne.in_Nz) goto 1003

      allocate(buffer(nx,ny,nz))

      read(11,err=1004) buffer
      call scatter3d(buffer,U(:,:,:,1))
      read(11,err=1004) buffer
      call scatter3d(buffer,U(:,:,:,2))
      read(11,err=1004) buffer
      call scatter3d(buffer,U(:,:,:,3))
      read(11,err=1004) buffer
      call scatter3d(buffer,temp(:,:,:))
      read(11,err=1004) buffer
      call scatter3d(buffer,pres(:,:,:))

      CALL varnc_to_varc ()
      deallocate(buffer)
      CLOSE (11)


      RETURN

 1002 Print*,in_Nvar,' leida ',Nd,' esperada'
      print*,'ERROR DETECTADO LEYENDO EN NUMERO DE VARIABLES,CAMPOS'
      stop
 1003 print*,in_Nx,in_Ny,in_Nz,' leidas ',Nx,Ny,Nz,' eperadas'
      print*,'ERROR DETECTADO LEYENDO EN DIMENSIONES, CAMPOS'
      stop
 1004 print*,'ERROR DETECTADO LEYENDO LOS CAMPOS'
      stop
      END SUBROUTINE LEER_CAMPOS

      SUBROUTINE leer_datos()
      use dimensiones
      use flow
      use inputname
      use tiempo
      use sgdmodel
      use consfiltro
      use fieldrec
      use vieu
      IMPLICIT NONE
      character (1) char1


      open(10,file='data.in',form='formatted')

      read(10,'(a1)')char1
      read(10,101)Nx
      read(10,101)Ny
      read(10,101)Nz
      read(10,102)Reynolds
      read(10,102)Prandtl
      read(10,107)sgs_model
      read(10,102)Prandtl_t
      read(10,103)inputgrid
      read(10,103)inputfield
      read(10,102)Mach
      read(10,102)gamma
      read(10,102)Tmax
      read(10,102)CFL
      read(10,102)Cvisc
      read(10,102)cfilt
      read(10,109)n_serie
      read(10,102)dgrava1
      read(10,102)tvieu1

        write(6,703)Reynolds,Prandtl,sgs_model,Prandtl_t &
                   ,inputgrid,inputfield,Mach,gamma &
                   ,Tmax,CFL,Cvisc,cfilt,n_serie,dgrava1, &
                    tvieu1 
 703  format('Reynolds number-------------: ',e10.4/ &
             'Prandtl number--------------: ',e10.4/ &
             'sgs model-------------------: ',a5/ &
             'Turbulent Prandtl number----: ',e10.4/ &
             'grid file-------------------: ',a32/ &
             'restart file----------------: ',a20/ &
             'Mach number-----------------: ',e10.4/ &
             'Gamma ----------------------: ',e10.4/ &
             'end time--------------------: ',e10.4/ &
             'CFL-------------------------: ',e10.4/ &
             'Cvisc-----------------------: ',e10.4/ &
             'Cfilt-----------------------: ',e10.4/ &
             'Serie-----------------------: ',i3/    &
             'T Gravado-------------------: ',e10.4/ &
             'T VIS. MIN/MAX--------------: ',e10.4)

 101  format(20x,i4)
 102  format(20x,e12.4)
 103  format(20x,a32)
 104  format(20x,i2)
 105  format(20x,i6)
 106  format(20x,a3)
 107  format(20x,a5)
 108  format(20x,a4)
 109  format(20x,i3)



      RETURN
      END
