MODULE RADAU
   USE, INTRINSIC :: ISO_C_BINDING

   INTEGER(C_INT) NEQF, NRDF, LWORKF, LIWORKF, AIPARF(1), IPARF, IOUTF, IDIDF, ITOLF, &
      IJAC, MLJAC, MUJAC, IMAS, MLMAS, MUMAS, IOUT, LJAC, LMAS, LE, IDID
   REAL(C_DOUBLE) RTOLF, ATOLF, RPARF, RTOL, ATOL, ARTOLF(1), AATOLF(1), ARPARF(1), H

   REAL(C_DOUBLE), ALLOCATABLE, TARGET :: WORKF(:), YF(:)
   INTEGER(C_INT), ALLOCATABLE, TARGET :: IWORKF(:)

   PRIVATE ALLOCATE_ARRAYS, DEALLOCATE_ARRAYS
CONTAINS

   SUBROUTINE RADAU_INIT(NF, FTOLL)
      IMPLICIT NONE

      INTEGER(C_INT), INTENT(IN) :: NF
      REAL(C_DOUBLE), INTENT(IN) :: FTOLL

      INTEGER I

      NEQF = NF
      NRDF = NF

      IJAC = 1
      MLJAC = NEQF
      IMAS = 0
      !MLMAS = 0
      !MUMAS = 0
      IOUT = 1
      LJAC = NEQF
      LMAS = 0
      LE = NEQF
      RTOLF = FTOLL
      ATOLF = 1.0D0*RTOLF
      ITOLf = 0

      H = 1.0D-6

      LWORKF = NEQF*(LJAC + LMAS + 3*LE + 12) + 20
      LIWORKF = 3*NEQF + 20

      CALL ALLOCATE_ARRAYS()

      !SET DEFAULT VALUES
      DO I = 1, 20
         IWORKF(I) = 0
         WORKF(I) = 0.D0
      END DO
   END SUBROUTINE RADAU_INIT

   SUBROUTINE ALLOCATE_ARRAYS()
      IMPLICIT NONE

      INTEGER(C_INT) ERR_ALLOC

      ALLOCATE (WORKF(LWORKF), IWORKF(LIWORKF), YF(NEQF), STAT=ERR_ALLOC)

      IF (ERR_ALLOC /= 0) THEN
         PRINT *, "ALLOCATION ERROR"
         PAUSE
         STOP
      END IF
   END SUBROUTINE ALLOCATE_ARRAYS

   SUBROUTINE DEALLOCATE_ARRAYS()
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE

      INTEGER(C_INT) ERR_DEALLOC

      DEALLOCATE (WORKF, IWORKF, YF, STAT=ERR_DEALLOC)

      IF (ERR_DEALLOC /= 0) THEN
         PRINT *, "DEALLOCATION ERROR"
         PAUSE
         STOP
      END IF
   END SUBROUTINE DEALLOCATE_ARRAYS

END MODULE RADAU