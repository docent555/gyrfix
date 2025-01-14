MODULE D02PVF_D02PCF_D02PDF
   USE, INTRINSIC :: ISO_C_BINDING

   INTEGER(C_INT) NEQP, LENWRK, METHOD, IFAILP
   REAL(C_DOUBLE) ZSTART, HSTART
   LOGICAL(C_BOOL) ERRASS

   REAL(C_DOUBLE), ALLOCATABLE, TARGET :: THRES(:), WORKP(:), PGOT(:), PPGOT(:), PMAX(:)

   PRIVATE ALLOCATE_ARRAYS, DEALLOCATE_ARRAYS, NEQP

CONTAINS

   SUBROUTINE D02PVF_D02PCF_D02PDF_INIT(NEE)
      IMPLICIT NONE

      INTEGER(C_INT) NEE, L

      NEQP = 4*NEE
      LENWRK = 32*NEQP
      ZSTART = 0.0D0
      ERRASS = .FALSE.
      HSTART = 0.0D0
      METHOD = 2
      IFAILP = 1

      CALL ALLOCATE_ARRAYS()

      DO L = 1, NEQP
         THRES(L) = 1.0D-8
      END DO
   END SUBROUTINE D02PVF_D02PCF_D02PDF_INIT

   SUBROUTINE ALLOCATE_ARRAYS()
      IMPLICIT NONE

      INTEGER(C_INT) ERR_ALLOC

      ALLOCATE (THRES(NEQP), WORKP(LENWRK), PGOT(NEQP), PPGOT(NEQP), PMAX(NEQP), STAT=ERR_ALLOC)

      IF (ERR_ALLOC /= 0) THEN
         PRINT *, "ALLOCATION ERROR"
         PAUSE
         STOP
      END IF
   END SUBROUTINE ALLOCATE_ARRAYS

   SUBROUTINE DEALLOCATE_ARRAYS()
      IMPLICIT NONE

      INTEGER(C_INT) ERR_DEALLOC

      DEALLOCATE (THRES, WORKP, PGOT, PPGOT, PMAX, STAT=ERR_DEALLOC)

      IF (ERR_DEALLOC /= 0) THEN
         PRINT *, "DEALLOCATION ERROR"
         PAUSE
         STOP
      END IF
   END SUBROUTINE DEALLOCATE_ARRAYS

END MODULE D02PVF_D02PCF_D02PDF
