        program fronterax
        implicit none
        integer nx,ny,nz
        integer i,j,k,nv
        real radio,r,u1,u2,y0,z0
        real,allocatable,DIMENSION(:,:,:)::xx,yy,zz,uu,uv,vv,ww,temp,press,cc
        real,allocatable,DIMENSION(:,:,:)::frontin

       open(69,file='grid.bin',form='unformatted')
        read(69)nv
        read(69)nx,ny,nz

        allocate(xx(nx,ny,nz))
        allocate(yy(nx,ny,nz))
        allocate(zz(nx,ny,nz))
        allocate(uu(nx,ny,nz))
        allocate(uv(nx,ny,nz))
        allocate(vv(nx,ny,nz))
        allocate(ww(nx,ny,nz))
        allocate(temp(nx,ny,nz))
        allocate(press(nx,ny,nz))
        allocate(cc(nx,ny,nz))
        allocate(frontin(ny,nz,6))
      !        real,DIMENSION(nx,ny,nz)::uv,vv,ww,temp
      !        real,DIMENSION(nx,ny,nz)::press,cc
      !        real,DIMENSION(ny,nz,6)::frontin
      !        real,DIMENSION(109,109,6)::frontin
           read(69)xx
           read(69)yy
           read(69)zz
       close(69)
radio=0.5
u1=0.07
u2=1.0
y0=yy(1,ny,1)/2.0
z0=zz(1,1,nz)/2.0
    DO k=1,nz
        DO j=1,ny
            DO i=1,nx

                r=sqrt((yy(i,j,k)-y0)**2.0+(zz(i,j,k)-z0)**2.0)
                uu(i,j,k)=0.5*((u1+u2)+(u1-u2)*tanh(10.0*(r/radio-radio/r)))
                uv(i,j,k)=uu(i,j,k)
                cc(i,j,k)=1.0
                IF(i > 10) then 
                    uu(i,j,k) = uu(0,1,1)
                endif
            END DO
         END DO
      END DO
      
    DO k=1,nz
        DO j=1,ny
            DO i=1,nx
                vv(i,j,k)=0.0
                ww(i,j,k)=0.0
                press(i,j,k)=1.0
                temp(i,j,k)=1.0
            ENDDO
        ENDDO
    ENDDO
    
close(61)
      DO k=1,nz
         DO j=1,ny
          frontin(j,k,1)=temp(1,j,k)
          frontin(j,k,2)=uv(1,j,k)
          frontin(j,k,3)=vv(1,j,k)
          frontin(j,k,4)=ww(1,j,k)
          frontin(j,k,5)=press(1,j,k)
          !frontin(j,k,5)=press(1,j,k)
          frontin(j,k,6)=cc(1,j,k)
         ENDDO
      ENDDO

      nv=5
      open(11,file='field_01.000',form='unformatted')
      write(11)nv
      write(11)nx,ny,nz
      write(11) uu
      write(11) vv
      write(11) ww
      write(11) temp
      write(11) press
      write(11) cc
      close(11)
      write(*,*)'Condicion de inicial creada->field_01.000'
      !write(*,*)frontin(55,55,5)
      open(11,file='frontx.in',form='unformatted')
      write(11)nv
      write(11)nx,ny,nz
      write(11)frontin
      close(11)
      write(*,*)'Condicion de frontera creada->frontx.in'


        deallocate(xx)
        deallocate(yy)
        deallocate(zz)
        deallocate(uu)
        deallocate(uv)
        deallocate(vv)
        deallocate(ww)
        deallocate(temp)
        deallocate(press)
        deallocate(cc)
        deallocate(frontin)
      END
