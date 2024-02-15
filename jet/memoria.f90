!-----------------------------------------------------------------------
      subroutine memoria()
!-----------------------------------------------------------------------

      use dimensiones
      use consderper
      use consdernper
      use consderper_f
      use consdernper_f
      use funciones
      IMPLICIT NONE
    


      allocate(axp(nx))
      allocate(bxp(nx))
      allocate(cxp(nx))
      allocate(ayp(ny))
      allocate(byp(ny))
      allocate(cyp(ny))
      allocate(azp(nz))
      allocate(bzp(nz))
      allocate(czp(nz))
      allocate(axnp(nx))
      allocate(bxnp(nx))
      allocate(cxnp(nx))
      allocate(aynp(ny))
      allocate(bynp(ny))
      allocate(cynp(ny))
      allocate(aznp(nz))
      allocate(bznp(nz))
      allocate(cznp(nz))
      allocate(axpf(nx))
      allocate(bxpf(nx))
      allocate(cxpf(nx))
      allocate(aypf(ny))
      allocate(bypf(ny))
      allocate(cypf(ny))
      allocate(azpf(nz))
      allocate(bzpf(nz))
      allocate(czpf(nz))
      allocate(axnpf(nx))
      allocate(bxnpf(nx))
      allocate(cxnpf(nx))
      allocate(aynpf(ny))
      allocate(bynpf(ny))
      allocate(cynpf(ny))
      allocate(aznpf(nz))
      allocate(bznpf(nz))
      allocate(cznpf(nz))
      allocate(u(nx))
      allocate(u2(nx))

      
      return
      end subroutine memoria
!-----------------------------------------------------------------------
      subroutine memoria_dealloc()
!-----------------------------------------------------------------------
      use dimensiones
      use consderper
      use consdernper
      use consderper_f
      use consdernper_f
      use funciones
      IMPLICIT NONE

      deallocate(axp)
      deallocate(bxp)
      deallocate(cxp)
      deallocate(ayp)
      deallocate(byp)
      deallocate(cyp)
      deallocate(azp)
      deallocate(bzp)
      deallocate(czp)
      deallocate(axnp)
      deallocate(bxnp)
      deallocate(cxnp)
      deallocate(aynp)
      deallocate(bynp)
      deallocate(cynp)
      deallocate(aznp)
      deallocate(bznp)
      deallocate(cznp)
      deallocate(axpf)
      deallocate(bxpf)
      deallocate(cxpf)
      deallocate(aypf)
      deallocate(bypf)
      deallocate(cypf)
      deallocate(azpf)
      deallocate(bzpf)
      deallocate(czpf)
      deallocate(axnpf)
      deallocate(bxnpf)
      deallocate(cxnpf)
      deallocate(aynpf)
      deallocate(bynpf)
      deallocate(cynpf)
      deallocate(aznpf)
      deallocate(bznpf)
      deallocate(cznpf)
      deallocate(u)
      deallocate(u2)
      return
      end subroutine memoria_dealloc
