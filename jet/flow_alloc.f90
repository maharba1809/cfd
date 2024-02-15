!-----------------------------------------------------------------------
      subroutine flow_alloc()
!-----------------------------------------------------------------------

      use dimensiones
      use consderper
      use consdernper
      use derivtools
      use jacobtools
      use mallagrid
      use derivvel
      use mflujos
      use dmflujos
      use viscosidades
      use velocidades
      use variables
      use right
      use tiempo
      use consderper_f
      use consdernper_f
      use sgdmodel
      use vorticidad
      use stat
      use fronterapl
      use bob
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
      allocate(du(nx))
      allocate(dv(ny))
      allocate(dw(nz))
      allocate(jbn(nx,ny,nz,11))
      allocate(x(nx,ny,nz,nd))
      allocate(dvel(nx,ny,nz,9))
      allocate(dconc(nx,ny,nz,nd))
      allocate(dtemp(nx,ny,nz,nd))
      allocate(e(nx,ny,nz))
      allocate(f(nx,ny,nz))
      allocate(g(nx,ny,nz))
      allocate(ev(nx,ny,nz))
      allocate(fv(nx,ny,nz))
      allocate(gv(nx,ny,nz))
      allocate(dere(nx,ny,nz))
      allocate(derf(nx,ny,nz))
      allocate(derg(nx,ny,nz))
      allocate(derev(nx,ny,nz))
      allocate(derfv(nx,ny,nz))
      allocate(dergv(nx,ny,nz))
      allocate(vis(nx,ny,nz))
      allocate(vis2(nx,ny,nz))
      allocate(vis6(nx,ny,nz))
      allocate(u(nx,ny,nz,nd))
      allocate(conc(nx,ny,nz))
      allocate(temp(nx,ny,nz))
      allocate(pres(nx,ny,nz))
      allocate(um(nx,ny,nz,nd))
      allocate(um0(nx,ny,nz,nd))
      allocate(um1(nx,ny,nz,nd))
      allocate(rs(nx,ny,nz))
      allocate(rsv(nx,ny,nz))
      allocate(dx(nx,ny,nz))
      allocate(dy(nx,ny,nz))
      allocate(dz(nx,ny,nz))
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
      allocate(amut(nx,ny,nz))
      allocate(dxmsgd(nx,ny,nz))
      allocate(dxpsgd(nx,ny,nz))
      allocate(dymsgd(nx,ny,nz))
      allocate(dypsgd(nx,ny,nz))
      allocate(dzmsgd(nx,ny,nz))
      allocate(dzpsgd(nx,ny,nz))
      allocate(dxyzsgd(nx,ny,nz))
      allocate(wx(nx,ny,nz))
      allocate(wy(nx,ny,nz))
      allocate(wz(nx,ny,nz))
      allocate(st(nx,ny,nz,15))
      allocate(frontx(2,nd,ny,nz))
      allocate(frontz(2,nd,nx,ny))
      allocate(fronty(2,nd,nx,nz))
      allocate(frontin(ny,nz,nd))
      allocate(imask(2,nx,ny))
      allocate(esplayer(ny,nz,nd))
      allocate(esource(nx))

      
      return
      end subroutine flow_alloc
!-----------------------------------------------------------------------
      subroutine flow_dealloc()
!-----------------------------------------------------------------------
      use dimensiones
      use consderper
      use consdernper
      use derivtools
      use jacobtools
      use mallagrid
      use derivvel
      use mflujos
      use dmflujos
      use viscosidades
      use velocidades
      use variables
      use right
      use tiempo
      use consderper_f
      use consdernper_f
      use sgdmodel
      use vorticidad
      use stat
      use fronterapl
      use bob
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
      deallocate(du)
      deallocate(dv)
      deallocate(dw)
      deallocate(jbn)
      deallocate(x)
      deallocate(dvel)
      deallocate(e)
      deallocate(f)
      deallocate(g)
      deallocate(dere)
      deallocate(derf)
      deallocate(derg)
      deallocate(derev)
      deallocate(derfv)
      deallocate(dergv)
      deallocate(ev)
      deallocate(fv)
      deallocate(gv)
      deallocate(vis)
      deallocate(vis2)
      deallocate(vis6)
      deallocate(u)
      deallocate(conc)
      deallocate(temp)
      deallocate(pres)
      deallocate(um)
      deallocate(um0)
      deallocate(um1)
      deallocate(rs)
      deallocate(rsv)
      deallocate(dx)
      deallocate(dy)
      deallocate(dz)
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
      deallocate(amut)
      deallocate(dxmsgd)
      deallocate(dxpsgd)
      deallocate(dymsgd)
      deallocate(dypsgd)
      deallocate(dzmsgd)
      deallocate(dzpsgd)
      deallocate(dxyzsgd)
      deallocate(wx)
      deallocate(wy)
      deallocate(wz)
      deallocate(st)
      deallocate(frontx)
      deallocate(frontz)
      deallocate(fronty)
      deallocate(frontin)
      deallocate(imask)
      deallocate(esplayer)
      deallocate(esource)
      return
      end subroutine flow_dealloc
