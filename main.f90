PROGRAM SYS15F
   USE, INTRINSIC :: ISO_C_BINDING
   USE FUN
   USE IFPORT

   IMPLICIT NONE

   INTEGER(C_INT) I, J, HOURS, MINUTES, SECONDS, RES
   REAL(C_DOUBLE) START_TIME, STOP_TIME, CALC_TIME
   COMPLEX(C_DOUBLE_COMPLEX) PC

   CALL INIT()

   WRITE (*, '(/)')

   START_TIME = DCLOCK()
   DO I = 1, NDTR
      !DCIR(1) = DCIR1
      !DCIR(2) = DCIR0 + (I - 1)*DCIRH

      !DTR(1) = DTR0 + (I - 1)*DTRH
      DTR(2) = DTR0 + (I - 1)*DTRH

      WRITE (PATH, '(F6.4,A)') DTR(2), '/'
      PRINT *, 'DTR1 = ', DTR(1),  'DTR2 = ', DTR(2)
      PRINT *, 'DCIR1 = ', DCIR(1), 'DCIR2 = ', DCIR(2)
      !WRITE (PATH, '(A)') './'          

      !WRITE (PATH, '(F7.5,A)') DCIR(2), '/'
      !PRINT *, 'DCIR = ', DCIR(1)

      RES = MAKEDIRQQ(PATH)

      CALL ODE4F()
      WRITE (*, '(/)')
      PRINT *, 'WRITING...'

      IF (INHER .EQ. .TRUE.) THEN
         F10 = F(1, NT)
         F20 = F(3, NT)
         F30 = F(5, NT)
         P10 = MOD(F(2, NT), 2*PI)
         P20 = MOD(F(4, NT), 2*PI)
         P30 = MOD(F(6, NT), 2*PI)
      END IF

      CALL WRITE_RESULTS()
   END DO
   STOP_TIME = DCLOCK()

   CALC_TIME = STOP_TIME - START_TIME

   HOURS = CALC_TIME/3600
   MINUTES = (CALC_TIME - HOURS*3600)/60
   SECONDS = CALC_TIME - HOURS*3600 - MINUTES*60

   PRINT *, 'EXECUTION TIME:', HOURS, 'H :', MINUTES, 'M :', SECONDS, 'S'
   WRITE (*, '(/)')

   PAUSE
   STOP
END PROGRAM
