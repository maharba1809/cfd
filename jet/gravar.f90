!______________________________________________________________
       SUBROUTINE GRAVA()
!______________________________________________________________

       use dimensiones
       use velocidades
       use fieldrec
       IMPLICIT NONE
       character (12) outputfield
       real, allocatable, dimension (:,:,:) :: buffer
       integer i,j,k,l,m
       


       allocate(buffer(nx,ny,nz))
       WRITE(sample,'(i3.3)')i_sample
       WRITE(serie,'(i2.2)')n_serie
       outputfield='field_'//serie//'.'//sample
       write(6,*)'ESCRIBIENDO ARCHIVO','  ',outputfield
       OPEN(11,file=outputfield,form='unformatted')
       WRITE(11)nd
       WRITE(11)nx,ny,nz
       call scatter3d(u(:,:,:,1),buffer)
       WRITE(11)buffer
       call scatter3d(u(:,:,:,2),buffer)
       WRITE(11)buffer
       call scatter3d(u(:,:,:,3),buffer)
       WRITE(11)buffer
       WRITE(11)temp
       WRITE(11)pres
       i_sample=i_sample+1
       igrava=igrava+1
       dgrava=float(igrava)*dgrava1
       CLOSE(11)
       RETURN
       END SUBROUTINE GRAVA
