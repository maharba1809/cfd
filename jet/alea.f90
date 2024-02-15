!_________________________________________________
      FUNCTION ALEA (R)
!_________________________________________________

!      use random
!      implicit none
!      real :: R ,rand,IY0,IY1,alea
      DATA IA1, IA0, IA1MA0 /1536, 1029, 507/
      DATA IC /1731/
      DATA IX1, IX0 /0, 0/
!***FIRST EXECUTABLE STATEMENT  RAND
      IF (R.LT.0.) GO TO 10
      IF (R.GT.0.) GO TO 20

      IY0 = IA0*IX0
      IY1 = IA1*IX1 + IA1MA0*(IX0-IX1) + IY0
      IY0 = IY0 + IC
      IX0 = MOD (IY0, 2048)
      IY1 = IY1 + (IY0-IX0)/2048
      IX1 = MOD (IY1, 2048)
!
 10   RAND = IX1*2048 + IX0
      alea = RAND / 4194304.
      RETURN 
!
 20   IX1 = AMOD(R,1.)*4194304. + 0.5
      IX0 = MOD (IX1, 2048)
      IX1 = (IX1-IX0)/2048
      GO TO 10

      END FUNCTION ALEA

!_________________________________________________
      SUBROUTINE INI_ALEA ()
!_________________________________________________

      use random
      implicit none
!
      IA1=1536.
      IA0=1029.
      IA1MA0=507.
      IC=1731.
      IX1=0.
      IX0=0.
      s1=0.
      s2=0.
      s3=0.

      RETURN
      END SUBROUTINE INI_ALEA
