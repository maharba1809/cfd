      SUBROUTINE INI_CONSTANTES ()
      use dimensiones
      use consadim
      use flow
      IMPLICIT NONE

      c1=gamma*mach*mach/2.e0
      c2=gamma-1.e0
      c3=1./(gamma*mach*mach)
      c4=1./reynolds
      c5=gamma/(reynolds*prandtl)
      c6=gamma*mach*mach


      return
      end

