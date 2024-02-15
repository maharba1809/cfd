!-----------------------------------------------------------------------
      program ehecatl
!-----------------------------------------------------------------------

      use dimensiones
      use tiempo
      use deltas
      use ranaleo
      IMPLICIT NONE
      idum=565878
      nd=5
      ne=5
      itr=20
      CALL LEER_DATOS ()
      CALL FLOW_ALLOC ()
      CALL INI_CONSTANTES ()
      CALL LEER_MALLA ()
      CALL LEER_CAMPOS ()
      CALL INI_SOURCE()
      CALL INIESTADISTICAS()
      CALL INIFRONTERA_X_VAL ()
      CALL INI_ESPONJA()

      deltax=float(nx-1)
      deltay=float(ny-1)
      deltaz=float(nz-1)

      CALL INIDERIV ()
      CALL INIFILTRO ()
      CALL INITSTEP ()
      CALL INISGDM()
      CALL JACOBEANO ()
      CALL RUNGEK ()
      CALL FINESTADISTICAS()
      CALL tp_maker ()
      CALL FLOW_DEALLOC ()

      PRINT *,'FIN DEL CALCULO EN EL TIEMPO',' ',dto,' ','ITERACION',' ',it
      END

