!************************************************************************
!   GENERADOR DE MALLAS USANDO LA TRANSFORMACION DE 'TANH'
!	(LE & MOIN)  "SEGUNDO"
!   nz(x): numero de puntos del dominio
!   lz(x): longitud del dominio
!   formulas de transformacion:
!	z1=e(1)*(1-tanh(g(1)*(e(1)-z')/tanh(g(1)*e(1))
!   z'=lz/(nz-1)*i i:iteracion
!    e1=eta,(debe ser la mitad de lz(x) <h>)
!    g(1)=gama,Indica la abertura del mallado.
!    Es una version modificada del generador para capa de mezcla. Ob
!    tiene una malla para capa limite en un canal de area transversal
!    rectangular.
!**********************************************************************
        program mallado
        implicit none
        real dx,dy,dz,ly,lx,lz,e1,e2,g1,g2
        integer nv,nx,ny,nz,i,j,k
        integer nx1,nx2,nx3
        real potx1,potx2,potx3
        real,allocatable,dimension(:,:,:)::xx,yy,zz
        real,allocatable,dimension(:)::x,y,z
        character(2)::regularx,regulary,regularz,escribir
    
        regularx='no'
        regulary='no'
        regularz='no'
        escribir='si'
        lx=11.0
        ly=5.5
        lz=5.5
        write(*,*)'lx,ly,lz'
        write(*,*)lx,ly,lz

        nx=150
        ny=109
        nz=109
        nv=5
        write(*,*)
        write(*,*)'nodos'
        write(*,*)nx,ny,nz

        allocate(x(nx))
        allocate(y(ny))
        allocate(z(nz))
        allocate(xx(nx,ny,nz))
        allocate(yy(nx,ny,nz))
        allocate(zz(nx,ny,nz))

        dx=lx/(FLOAT(nx-1))
        dy=ly/(FLOAT(ny-1))
        dz=lz/(FLOAT(nz-1))
!--------------------------------------------------------------
        if(regulary.eq.'si')then
           write(*,*)'malla regular en y'
           DO i=1,ny
              y(i)=dy*FLOAT(i-1)
           END DO
        else
           write(*,*)'malla irregular en y'
           e1=ly/2.
           g1=0.6
           write(*,*)'Contraccion y'
!          read(*,*)g1
           write(*,*)g1
           DO j=1,ny
!           y(j)=e1*(1.-TANH(g1*(e1-FLOAT(j-1)*dy))/TANH(g1*e1))
           y(j)=e1-1./g1*atanh((1.-float(j-1)*dy/e1)*tanh(g1*e1))
           END DO
        end if
!--------------------------------------------------------------
        if(regularz.eq.'si')then
           write(*,*)'malla regular en z'
           DO i=1,nz
              z(i)=dz*FLOAT(i-1)
           END DO
        else
           e2=lz/2.
           g2=0.6
           write(*,*)'Contraccion z'
!          read(*,*)g2
           write(*,*)g2
           DO k=1,nz
!           z(k)=e2*(1.-TANH(g2*(e2-FLOAT(k-1)*dz))/TANH(g2*e2))
              z(k)=e2-1./g2*atanh((1.-float(k-1)*dz/e2)*tanh(g2*e2))
           END DO
        end if
!--------------------------------------------------------------
        if(regularx.eq.'si')then
           write(*,*)'malla regular en x'
           DO i=1,nx
              x(i)=dx*FLOAT(i-1)
           ENDDO
        else
           write(*,*)'malla irregular en x'
           nx1=50
           nx2=35
           nx3=nx-nx1-nx2
           write(*,*)
           write(*,*)'nx1,nx2,nx3'
           write(*,*)nx1,nx2,nx3
           potx1=1.3
           potx2=1.1
           potx3=1.2
           write(*,*)
           write(*,*)'Potencia x1 x2 x3'
           write(*,*)potx1,potx2,potx3

           DO i=1,nx1
              x(i)=(1.5)*(FLOAT(i-1)/float(nx1-1))**potx1
           ENDDO
           DO i=1,nx2
              x(i+nx1)=(2.0)*(FLOAT(i)/float(nx2))**potx2+x(nx1)
           ENDDO
           DO i=1,nx3
              x(i+nx1+nx2)=(lx-x(nx1+nx2))*(FLOAT(i)/float(nx3))**potx3+x(nx2+nx1)
           ENDDO   

        end if


!--------------------------------------------------------------
        DO k=1,nz
           DO j=1,ny
              DO i=1,nx
              xx(i,j,k) = x(i)
              yy(i,j,k) = y(j)
              zz(i,j,k) = z(k)
              ENDDO
           ENDDO
        ENDDO
      if(escribir.eq.'si')then
      open(54,file='jet_malla.dat')
        write(54,*)'variables="x","y","z"'
        write(54,*)'zone i=',nx,'j=',ny,'k=',nz,'datapacking=point'
        DO k=1,nz
           DO j=1,ny
              DO i=1,nx
              write(54,*)x(i),y(j),z(k)
              ENDDO
           ENDDO
        ENDDO
      close(54)
      Write(*,*)'Malla creada--> ass_malla.dat'
      end if

      open(11,file='grid.bin',form='unformatted')
      write(11)nv
      write(11)nx,ny,nz
    
      WRITE(11)xx
      WRITE(11)yy
      WRITE(11)zz
      close(11)
      Write(*,*)'Malla creada--> grid.bin'
      deallocate(xx)
      deallocate(yy)
      deallocate(zz)
      deallocate(x)
      deallocate(y)
      deallocate(z)
      END
