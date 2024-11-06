MODULE FUN
   USE, INTRINSIC :: ISO_C_BINDING
   USE IFCORE
   USE IFPORT
   USE CSPLINE

   INTEGER(C_INT) NE, NEQP, NEQF, NT, NZ, FREQ_OUT, L, INHARM, NDCIR, METH, NTEND, NDTR
   REAL(C_DOUBLE) ZEX, DZ, TEND, DTR(2), Q(3), ICU(2), TH(2), A(2), DCIR(2), DCIR1, R(2), &
      F0(6), DT, F10, F20, F30, P10, P20, P30, FTOL, PTOL, NHARM, DTRH, DTRB(2), DTR0
   COMPLEX(C_DOUBLE_COMPLEX) FP(2)
   LOGICAL(C_BOOL) WC, FOK, LENSM, BTOD, IATOI, INHER, SQR
   CHARACTER(LEN=9) PATH

   INTEGER(C_INT) BREAKNUM(3)
   REAL(C_DOUBLE) PHITMP0(3), PHITMP1(3)
   COMPLEX(C_DOUBLE_COMPLEX) FC, FCOMP(3)

   INTEGER(C_INT), ALLOCATABLE, TARGET :: IDXRE(:, :), IDXIM(:, :), IDXP(:, :)
   COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, TARGET :: U(:), MEAN(:), XI1(:), XI2(:), XI1_KR(:), XI2_KR(:)
   REAL(C_DOUBLE), ALLOCATABLE, TARGET :: TAX(:), ZAX(:), ETA(:, :), W(:, :), F(:, :), P(:, :), &
                                          PHI(:, :), PHIOS(:, :), WOS(:, :), &
                                          CL1(:), LHS1(:), RHS1(:), CL2(:), LHS2(:), RHS2(:), &
                                          REB(:), REC(:), RED(:), IMB(:), IMC(:), IMD(:)

   INCLUDE 'Z.INC'
   INCLUDE 'RE.INC'
   INCLUDE 'IM.INC'

   INCLUDE 'ZA_22-09-23.INC'
   INCLUDE 'REA_22-09-23.INC'
   INCLUDE 'IMA_22-09-23.INC'

   INCLUDE 'Z_ANG.INC'
   INCLUDE 'RE_ANG.INC'
   INCLUDE 'IM_ANG.INC'

   COMPLEX(C_DOUBLE_COMPLEX), PARAMETER :: IC = (0.0D0, 1.0D0)
   REAL(C_DOUBLE), PARAMETER :: PI = 2.0D0*DACOS(0.0D0)

   PRIVATE FREQ_OUT, TEND, Q, I, TH, A, R, F0, PITCH

CONTAINS

   SUBROUTINE READ_PARAM() BIND(C, NAME='READ_PARAM')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPORT
      IMPLICIT NONE

      NAMELIST /PARAM/ NE, TEND, ZEX, Q1, Q2, Q3, I1, I2, TH1, TH2, A1, A2, &
         R1, R2, F10, F20, F30, P10, P20, P30, DT, DZ, FTOL, PTOL, WC, FOK, INHARM, &
         DCIR1, DCIR2, INHER, METH, SQR, NZ, DTRH, DTRB, DTR1

      REAL(C_DOUBLE) Q1, Q2, Q3, I1, I2, TH1, TH2, A1, A2, DTR1, DTR2, R1, R2, DCIR1, DCIR2

      OPEN (UNIT=1, FILE='INPUT_FORTRAN.IN', STATUS='OLD', ERR=101)
      READ (UNIT=1, NML=PARAM, ERR=102)
      CLOSE (UNIT=1)

      Q(1) = Q1
      Q(2) = Q2
      Q(3) = Q3
      ICU(1) = I1
      ICU(2) = I2
      TH(1) = TH1
      TH(2) = TH2
      A(1) = A1
      A(2) = A2
      DTR(1) = DTR1
      !DTR(2) = DTR2
      R(1) = R1
      R(2) = R2
      DCIR(1) = DCIR1
      DCIR(2) = DCIR2

      WRITE (*, NML=PARAM)

      RETURN
101   PRINT *, "ERROR OF FILE OPEN"; PAUSE; STOP
102   PRINT *, 'ERROR OF READING FILE "INPUT_FORTRAN.IN"'; PAUSE; STOP
   END SUBROUTINE READ_PARAM

   SUBROUTINE INIT()
      USE D02PVF_D02PCF_D02PDF, ONLY: D02PVF_D02PCF_D02PDF_INIT
      USE D02NVF_D02M_N, ONLY: D02NVF_D02M_N_INIT
      USE DOPRF, ONLY: DOPRF_INIT
      USE DOPRP, ONLY: DOPRP_INIT
      USE RADAU

      IMPLICIT NONE

      CALL READ_PARAM()

      NEQP = 4*NE
      NEQF = 6

      !ZEX = 47.1343

      PRINT *, 'ICU1 = ', ICU(1), 'ICU2 = ', ICU(2)

      NHARM = DBLE(INHARM)

      NT = TEND/DT + 1
      IF (NZ > 0) THEN
         DZ = ZEX/(NZ - 1)
         PRINT *, 'NZ > 0'
         PRINT *, 'NZ = ', NZ
         PRINT *, 'DZ = ', DZ
      ELSE
         NZ = ZEX/DZ + 1
         PRINT *, 'NZ == 0'
         PRINT *, 'NZ = ', NZ
         PRINT *, 'DZ = ', DZ
      END IF

      IF (METH .EQ. 1) THEN

         CALL DOPRP_INIT(NE, PTOL)
         CALL DOPRF_INIT(NE, FTOL)

      ELSE IF (METH .EQ. 2) THEN

         CALL D02PVF_D02PCF_D02PDF_INIT(NE)

         F0(1) = F10
         F0(2) = P10
         F0(3) = F20
         F0(4) = P20
         F0(5) = F30
         F0(6) = P30

         CALL D02NVF_D02M_N_INIT(6, NT, DT, TEND, FTOL, F0)

      ELSE IF (METH .EQ. 3) THEN

         !CALL D02PVF_D02PCF_D02PDF_INIT(NE)
         CALL DOPRP_INIT(NE, PTOL)

      ELSE IF (METH .EQ. 4) THEN

         !CALL D02PVF_D02PCF_D02PDF_INIT(NE)
         CALL DOPRP_INIT(NE, PTOL)
         CALL RADAU_INIT(NEQF, FTOL)

      END IF

      CALL ALLOCATE_ARRAYS()

      F(1, 1) = F10
      F(2, 1) = P10
      F(3, 1) = F20
      F(4, 1) = P20
      F(5, 1) = F30
      F(6, 1) = P30

      DO I = 1, NT
         TAX(I) = (I - 1)*DT
      END DO

      DO I = 1, NZ
         ZAX(I) = (I - 1)*DZ
      END DO

      CALL CMPLXSPLINE(918, ZEX, Z_ANG, RE_ANG, IM_ANG, REB, REC, RED, IMB, IMC, IMD)

      CALL CALC_U(U, ZEX, NZ, ZAX)

      DO I = 1, 2
         IDXRE(I, :) = (/2*(I - 1)*NE + 1:(2*I - 1)*NE/)
         IDXIM(I, :) = (/(2*I - 1)*NE + 1:2*I*NE/)
      END DO

      IDXP(1, :) = (/1:NE/)
      IDXP(2, :) = (/NE + 1:2*NE/)

      DO I = 1, NE
         P(I, 1) = DREAL(CDEXP(IC*(I - 1)/DBLE(NE)*2*PI))
         P(NE + I, 1) = DIMAG(CDEXP(IC*(I - 1)/DBLE(NE)*2*PI))
         P(2*NE + I, 1) = DREAL(CDEXP(IC*(I - 1)/DBLE(NE)*2*PI))
         P(3*NE + I, 1) = DIMAG(CDEXP(IC*(I - 1)/DBLE(NE)*2*PI))
      END DO

      !IF (DCIRH .NE. 0) THEN
      !   NDCIR = (DCIRB(2) - DCIRB(1))/DCIRH + 1
      !ELSE
      !   NDCIR = 1
      !END IF

      !DCIR0 = DCIRB(1)

      IF (DTRH .NE. 0) THEN
         NDTR = (DTRB(2) - DTRB(1))/DTRH + 1
      ELSE
         NDTR = 1
      END IF

      DTR0 = DTRB(1)

      W(:, :) = 0

   END SUBROUTINE INIT

   FUNCTION SQUVAL(ZZ)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE, C_DOUBLE_COMPLEX, C_INT
      IMPORT, ONLY:ZEX, ZA, REA, IMA
      IMPLICIT NONE

      REAL(C_DOUBLE), INTENT(IN) :: ZZ

      COMPLEX(C_DOUBLE_COMPLEX) SQUVAL
      REAL(C_DOUBLE) Z, RE, IM, DZ, Z1
      INTEGER(C_INT) L

      DZ = 0.280211
      Z = ZZ/ZEX*185.5D0
      L = Z/DZ

      IF (L .EQ. 0) THEN
         L = 2
      ELSEIF (L .GE. 662) THEN
         L = 662
      ELSE
         IF ((Z - ZA(L)) .GT. 0.5*DZ) L = L + 1
      END IF

      Z1 = ZA(L)
      Z = Z - 8.5D0

      RE = REA(L - 1) + ((-REA(L - 1) + REA(L))*(DZ + Z - Z1))/DZ + ((REA(L - 1)/2.0D0 - REA(L) + &
                                                                      REA(L + 1)/2.0D0)*(Z - Z1)*(DZ + Z - Z1))/DZ/DZ
      IM = IMA(L - 1) + ((-IMA(L - 1) + IMA(L))*(DZ + Z - Z1))/DZ + ((IMA(L - 1)/2.0D0 - IMA(L) + &
                                                                      IMA(L + 1)/2.0D0)*(Z - Z1)*(DZ + Z - Z1))/DZ/DZ
      !!!NE RABOTAET
      !RE = ((REA(L - 1) - 2*REA(L) + REA(L + 1))*Z**2)/(2.*DZ**2) &
      !     + (Z*(DZ*(-REA(L - 1) + REA(L + 1)) - 2*(REA(L - 1) - 2*REA(L) + REA(L + 1))*Z1))/(2.*DZ**2) + &
      !     -(2*DZ**2*REA(L) + DZ*(REA(L - 1) - REA(L + 1))*Z1 + (REA(L - 1) - 2*REA(L) + REA(L + 1))*Z1**2)/(2.*DZ**2)
      !IM = ((IMA(L - 1) - 2*IMA(L) + IMA(L + 1))*Z**2)/(2.*DZ**2) + &
      !     (Z*(DZ*(-IMA(L - 1) + IMA(L + 1)) - 2*(IMA(L - 1) - 2*IMA(L) + IMA(L + 1))*Z1))/(2.*DZ**2) + &
      !     -(2*DZ**2*IMA(L) + DZ*(IMA(L - 1) - IMA(L + 1))*Z1 + (IMA(L - 1) - 2*IMA(L) + IMA(L + 1))*Z1**2)/(2.*DZ**2)

      SQUVAL = DCMPLX(RE, IM)

   END FUNCTION SQUVAL

   FUNCTION UVAL(ZZ)

      IMPLICIT NONE

      REAL(C_DOUBLE), INTENT(IN) :: ZZ

      COMPLEX(C_DOUBLE_COMPLEX) UVAL
      REAL(C_DOUBLE) Z, RE, IM, D
      INTEGER(C_INT) L

      Z = ZZ/ZEX*185.5 - 8.5
      L = (Z + 8.5)/0.28021 + 1
      D = Z - ZA(L)

      !PRINT *, Z, L, D

      IF (D .GT. 0.0 .AND. L /= 663) THEN
         RE = (REA(L)*ZA(L + 1) - REA(L + 1)*ZA(L))/(ZA(L + 1) - ZA(L)) + &
              (REA(L + 1) - REA(L))/(ZA(L + 1) - ZA(L))*Z
         IM = (IMA(L)*ZA(L + 1) - IMA(L + 1)*ZA(L))/(ZA(L + 1) - ZA(L)) + &
              (IMA(L + 1) - IMA(L))/(ZA(L + 1) - ZA(L))*Z
      ELSE IF (D .LT. 0.0 .AND. L /= 1) THEN
         RE = (REA(L - 1)*ZA(L) - REA(L)*ZA(L - 1))/(ZA(L) - ZA(L - 1)) + &
              (REA(L) - REA(L - 1))/(ZA(L) - ZA(L - 1))*Z
         IM = (IMA(L - 1)*ZA(L) - IMA(L)*ZA(L - 1))/(ZA(L) - ZA(L - 1)) + &
              (IMA(L) - IMA(L - 1))/(ZA(L) - ZA(L - 1))*Z
      ELSE
         RE = REA(L)
         IM = IMA(L)
      END IF

      UVAL = DCMPLX(RE, IM)

   END FUNCTION UVAL

   FUNCTION SQUVAL22(ZZ)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE, C_DOUBLE_COMPLEX, C_INT
      IMPORT, ONLY:ZEX, ZA22, REA22, IMA22
      IMPLICIT NONE

      REAL(C_DOUBLE), INTENT(IN) :: ZZ

      COMPLEX(C_DOUBLE_COMPLEX) SQUVAL22
      REAL(C_DOUBLE) Z, RE, IM, DZ, ZL
      INTEGER(C_INT) L, NS

      Z = ZZ/ZEX*52.43D0

      L = 1
      DO WHILE (ZA22(L) < Z)
         L = L + 1
      END DO

      ZL = ZA22(L)
      IF (L .EQ. 1) THEN
         L = 2
         ZL = ZA22(L)
         DZ = (ZA22(L + 1) - ZA22(L - 1))/2.0D0
      ELSEIF (L .GE. 1122) THEN
         L = 1122
         ZL = ZA22(L)
         DZ = (ZA22(L + 1) - ZA22(L - 1))/2.0D0
      ELSE
         IF ((Z - ZL)/(ZA22(L + 1) - ZL) .GT. 0.5) THEN
            L = L + 1
            ZL = ZA22(L)
         END IF
         DZ = (ZA22(L + 1) - ZA22(L - 1))/2.0D0
      END IF

      !Z = Z - 0.0D0

      RE = REA22(L - 1) + ((-REA22(L - 1) + REA22(L))*(DZ + Z - ZL))/DZ + ((REA22(L - 1)/2.0D0 - REA22(L) + &
                                                                            REA22(L + 1)/2.0D0)*(Z - ZL)*(DZ + Z - ZL))/DZ/DZ
      IM = IMA22(L - 1) + ((-IMA22(L - 1) + IMA22(L))*(DZ + Z - ZL))/DZ + ((IMA22(L - 1)/2.0D0 - IMA22(L) + &
                                                                            IMA22(L + 1)/2.0D0)*(Z - ZL)*(DZ + Z - ZL))/DZ/DZ
      !!!NE RABOTAET
      !RE = ((REA22(L - 1) - 2*REA22(L) + REA22(L + 1))*Z**2)/(2.*DZ**2) &
      !     + (Z*(DZ*(-REA22(L - 1) + REA22(L + 1)) - 2*(REA22(L - 1) - 2*REA22(L) + REA22(L + 1))*ZL))/(2.*DZ**2) + &
      !     -(2*DZ**2*REA22(L) + DZ*(REA22(L - 1) - REA22(L + 1))*ZL + (REA22(L - 1) - 2*REA22(L) + REA22(L + 1))*ZL**2)/(2.*DZ**2)
      !IM = ((IMA22(L - 1) - 2*IMA22(L) + IMA22(L + 1))*Z**2)/(2.*DZ**2) + &
      !     (Z*(DZ*(-IMA22(L - 1) + IMA22(L + 1)) - 2*(IMA22(L - 1) - 2*IMA22(L) + IMA22(L + 1))*ZL))/(2.*DZ**2) + &
      !     -(2*DZ**2*IMA22(L) + DZ*(IMA22(L - 1) - IMA22(L + 1))*ZL + (IMA22(L - 1) - 2*IMA22(L) + IMA22(L + 1))*ZL**2)/(2.*DZ**2)

      SQUVAL22 = DCMPLX(RE, IM)
   END FUNCTION SQUVAL22

   FUNCTION UVAL22(ZZ)

      IMPLICIT NONE

      REAL(C_DOUBLE), INTENT(IN) :: ZZ

      COMPLEX(C_DOUBLE_COMPLEX) UVAL22
      REAL(C_DOUBLE) Z, RE, IM, D
      INTEGER(C_INT) L

      Z = ZZ/ZEX*52.43

      L = 1
      DO WHILE (ZA22(L) < Z)
         L = L + 1
      END DO
      D = Z - ZA22(L)

      !PRINT *, Z, L, D

      IF (D .GT. 0.0 .AND. L /= 1123) THEN
         RE = (REA22(L)*ZA22(L + 1) - REA22(L + 1)*ZA22(L))/(ZA22(L + 1) - ZA22(L)) + &
              (REA22(L + 1) - REA22(L))/(ZA22(L + 1) - ZA22(L))*Z
         IM = (IMA22(L)*ZA22(L + 1) - IMA22(L + 1)*ZA22(L))/(ZA22(L + 1) - ZA22(L)) + &
              (IMA22(L + 1) - IMA22(L))/(ZA22(L + 1) - ZA22(L))*Z
      ELSE IF (D .LT. 0.0 .AND. L /= 1) THEN
         RE = (REA22(L - 1)*ZA22(L) - REA22(L)*ZA22(L - 1))/(ZA22(L) - ZA22(L - 1)) + &
              (REA22(L) - REA22(L - 1))/(ZA22(L) - ZA22(L - 1))*Z
         IM = (IMA22(L - 1)*ZA22(L) - IMA22(L)*ZA22(L - 1))/(ZA22(L) - ZA22(L - 1)) + &
              (IMA22(L) - IMA22(L - 1))/(ZA22(L) - ZA22(L - 1))*Z
      ELSE
         RE = REA22(L)
         IM = IMA22(L)
      END IF

      UVAL22 = DCMPLX(RE, IM)

   END FUNCTION UVAL22

   SUBROUTINE ALLOCATE_ARRAYS()
      !USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE

      INTEGER(C_INT) ERR_ALLOC

      ALLOCATE (F(6, NT), P(NEQP, NZ), U(NZ), TAX(NT), ZAX(NZ), MEAN(NZ), ETA(2, NT), W(3, NT), &
                IDXRE(2, NE), IDXIM(2, NE), WOS(3, NT), PHI(3, NT), PHIOS(3, NT), IDXP(2, NE), &
                CL1(NT), LHS1(NT), RHS1(NT), CL2(NT), LHS2(NT), RHS2(NT), XI1(NT), XI2(NT), XI1_KR(NT), XI2_KR(NT), &
                REB(NZ), REC(NZ), RED(NZ), IMB(NZ), IMC(NZ), IMD(NZ), STAT=ERR_ALLOC)

      IF (ERR_ALLOC /= 0) THEN
         PRINT *, "ALLOCATION ERROR"
         PAUSE
         STOP
      END IF
   END SUBROUTINE ALLOCATE_ARRAYS

   SUBROUTINE DEALLOCATE_ARRAYS()
      !USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE

      INTEGER(C_INT) ERR_DEALLOC

      DEALLOCATE (F, P, U, TAX, ZAX, MEAN, ETA, W, &
                  IDXRE, IDXIM, WOS, PHI, PHIOS, IDXP, CL1, LHS1, RHS1, CL2, LHS2, RHS2, STAT=ERR_DEALLOC)

      IF (ERR_DEALLOC /= 0) THEN
         PRINT *, "DEALLOCATION ERROR"
         PAUSE
         STOP
      END IF
   END SUBROUTINE DEALLOCATE_ARRAYS

   SUBROUTINE WRITE_RESULTS()
      IMPLICIT NONE

      INTEGER I, J
      IF (METH .NE. 2) THEN
         IF (WC .EQ. .TRUE.) THEN
            W(:, 1) = 0
            DO I = 2, NTEND
               DO J = 1, 3
                  !W(J, I - 1) = DIMAG(LOG(F(2*J - 1, I)*CDEXP(IC*F(2*J, I))/(F(2*J - 1, I - 1)*CDEXP(IC*F(2*J, I - 1)))))/DT
                  W(J, I) = (F(2*J, I) - F(2*J, I - 1))/DT
               END DO
            END DO
            PRINT *, 'FREQUENCY CALCULATED FROM PHASE. ( WC = ', WC, ')'
         ELSEIF (WC .EQ. .FALSE.) THEN
            CALL FREQ()
            PRINT *, 'FREQUENCY CALCULATED FROM RHS. ( WC = ', WC, ')'
         END IF
      END IF

      PHI(:, 1) = 0; 
      DO I = 2, NTEND
         DO J = 1, 3
            PHI(J, I) = PHI(J, I - 1) + DIMAG(LOG(F(2*J - 1, I)*CDEXP(IC*F(2*J, I))/(F(2*J - 1, I - 1)*CDEXP(IC*F(2*J, I - 1)))))
         END DO
      END DO

      BREAKNUM(:) = 0
      FCOMP(1) = F(2*1 - 1, 1)*CDEXP(IC*F(2*1, 1))
      FCOMP(2) = F(2*2 - 1, 1)*CDEXP(IC*F(2*2, 1))
      FCOMP(3) = F(2*3 - 1, 1)*CDEXP(IC*F(2*3, 1))
      PHITMP0(:) = DATAN2(DIMAG(FCOMP(:)), DREAL(FCOMP(:)))
      !PHITMP0(:) = DATAN2(DIMAG(F(:, 1)), DREAL(F(:, 1)))
      PHIOS(:, 1) = PHITMP0(:)
      DO I = 2, NTEND
         DO J = 1, 3
            FC = F(2*J - 1, I)*CDEXP(IC*F(2*J, I))
            PHITMP1(J) = DATAN2(DIMAG(FC), DREAL(FC))
            IF ((PHITMP1(J) - PHITMP0(J)) .GT. PI) BREAKNUM(J) = BREAKNUM(J) - 1
            IF ((PHITMP1(J) - PHITMP0(J)) .LT. -PI) BREAKNUM(J) = BREAKNUM(J) + 1
            PHIOS(J, I) = PHITMP1(J) + 2.*PI*BREAKNUM(J)
            !PHIOS(J, I) = PHITMP1(J)
            PHITMP0(J) = PHITMP1(J)
         END DO
      END DO

      DO I = 1, NTEND - 1
         DO J = 1, 3
            WOS(J, I) = (PHIOS(J, I + 1) - PHIOS(J, I))/DT
         END DO
      END DO

      WRITE (*, '(/)')

      !PAUSE

      OPEN (3, FILE=TRIM(ADJUSTL(PATH))//'CL1.DAT')
      DO I = 1, NTEND
         WRITE (3, '(5F14.6,A)') TAX(I), CL1(I), LHS1(I), RHS1(I), ABS(CL1(I)/LHS1(I))*100, ' %'
      END DO
      CLOSE (3)

      OPEN (3, FILE=TRIM(ADJUSTL(PATH))//'CL2.DAT')
      DO I = 1, NTEND
         WRITE (3, '(5F14.6,A)') TAX(I), CL2(I), LHS2(I), RHS2(I), ABS(CL2(I)/LHS2(I))*100, ' %'
      END DO
      CLOSE (3)

      OPEN (1, FILE=TRIM(ADJUSTL(PATH))//'F.DAT')
      DO I = 1, NTEND
         !WRITE (1, '(4E17.8)') TAX(I), DABS(F(1, I)), DABS(F(3, I)), DABS(F(5, I))
         WRITE (1, '(4F14.6)') TAX(I), F(1, I), F(3, I), F(5, I)
      END DO
      CLOSE (1)

      OPEN (13, FILE=TRIM(ADJUSTL(PATH))//'FCMPLX.DAT')
      DO I = 1, NTEND
         FCOMP(1) = F(2*1 - 1, I)*CDEXP(IC*F(2*1, I))
         FCOMP(2) = F(2*2 - 1, I)*CDEXP(IC*F(2*2, I))
         FCOMP(3) = F(2*3 - 1, I)*CDEXP(IC*F(2*3, I))
         WRITE (13, '(7F14.6)') TAX(I), DREAL(FCOMP(1)), DIMAG(FCOMP(1)), DREAL(FCOMP(2)), DIMAG(FCOMP(2)), &
            DREAL(FCOMP(3)), DIMAG(FCOMP(3))
      END DO
      CLOSE (13)

      OPEN (2, FILE=TRIM(ADJUSTL(PATH))//'E.DAT')
      DO I = 1, NTEND
         WRITE (2, '(3F14.6)') TAX(I), ETA(1, I), ETA(2, I)
      END DO
      CLOSE (2)
    
      OPEN (3, FILE=TRIM(ADJUSTL(PATH))//'W.DAT')
      DO I = 1, NTEND
         WRITE (3, '(4F14.6)') TAX(I), W(1, I), W(2, I), W(3, I)
      END DO
      CLOSE (3)

      OPEN (1, FILE=TRIM(ADJUSTL(PATH))//'P.DAT')
      DO I = 1, NTEND
         !WRITE (1, '(4E17.8)') TAX(I), PHI(1, I), PHI(2, I), PHI(3, I)
         WRITE (1, '(4F14.6)') TAX(I), F(2, I), F(4, I), F(6, I)
      END DO
      CLOSE (1)

      OPEN (1, FILE=TRIM(ADJUSTL(PATH))//'POS.DAT')
      DO I = 1, NTEND
         WRITE (1, '(4F14.6)') TAX(I), PHIOS(1, I), PHIOS(2, I), PHIOS(3, I)
      END DO
      CLOSE (1)

      OPEN (3, FILE=TRIM(ADJUSTL(PATH))//'WOS.DAT')
      DO I = 1, NTEND - 1
         WRITE (3, '(4F14.6)') TAX(I + 1), WOS(1, I), WOS(2, I), WOS(3, I)
      END DO
      CLOSE (3)

      OPEN (1, FILE=TRIM(ADJUSTL(PATH))//'XI1.DAT')
      DO I = 1, NTEND - 1
         WRITE (1, '(3F14.6)') TAX(I + 1), XI1(I)
      END DO
      CLOSE (1)

      OPEN (1, FILE=TRIM(ADJUSTL(PATH))//'XI2.DAT')
      DO I = 1, NTEND - 1
         WRITE (1, '(3F14.6)') TAX(I + 1), XI2(I)
      END DO
      CLOSE (1)

      OPEN (1, FILE=TRIM(ADJUSTL(PATH))//'XI1_KR.DAT')
      DO I = 1, NTEND - 1
         WRITE (1, '(3F14.6)') TAX(I + 1), XI1_KR(I)
      END DO
      CLOSE (1)

      OPEN (1, FILE=TRIM(ADJUSTL(PATH))//'XI2_KR.DAT')
      DO I = 1, NTEND - 1
         WRITE (1, '(3F14.6)') TAX(I + 1), XI2_KR(I)
      END DO
      CLOSE (1)

      !CALL WRITE_PARAM(PATH)

      RETURN

101   PRINT *, 'ERROR OF FILE OPEN.'
      PAUSE
      STOP
102   PRINT *, 'ERROR OF FILE READING.'
      PAUSE
      STOP
103   PRINT *, 'ERROR OF FILE WRITING.'
      PAUSE
      STOP
   END SUBROUTINE WRITE_RESULTS

   SUBROUTINE INIT_DOPR_P()
      USE DOPRP, ONLY: RPARP, IPARP, ITOLP, RTOLP, ATOLP, IWORKP, WORKP, PTOL, &
                       ARTOLP, AATOLP, ARPARP, AIPARP
      IMPORT, ONLY:NZ, DZ

      IMPLICIT NONE

      RPARP = 0.0
      IPARP = 0
      ITOLP = 0
      RTOLP = PTOL
      ATOLP = RTOLP
      IWORKP(:) = 0
      WORKP(:) = 0.0D0
      WORKP(6) = DZ

      ARTOLP(1) = RTOLP
      AATOLP(1) = ATOLP
      ARPARP(1) = RPARP
      AIPARP(1) = IPARP
   END SUBROUTINE INIT_DOPR_P

   SUBROUTINE SOLVEF_RADAU()
      USE, INTRINSIC :: ISO_C_BINDING
      USE RADAU
      IMPLICIT NONE

      REAL(C_DOUBLE) T, XOUTF
      INTEGER(C_INT) ITF, J

      REAL(C_DOUBLE) PEX(NEQP)

      COMMON/INTERNF/XOUTF, ITF
      !COMMON/INTERNF/ ITF

      RPARF = 0.0
      IPARF = 0
      T = TAX(1)
      XOUTF = T
      ITF = 0
      YF = F(:, 1)
      !ITOLF = 0
      !RTOLF = FTOL
      !ATOLF = RTOLF
      !IOUTF = 6

      ARTOLF(1) = RTOLF
      AATOLF(1) = ATOLF
      ARPARF(1) = RPARF
      AIPARF(1) = IPARF

      IWORKF(:) = 0
      WORKF(:) = 0.0D0

      IWORKF(5) = 6

      CALL RADAU5(NEQF, DFDT_RADAU, T, YF, TEND, H, &
                  RTOLF, ATOLF, ITOLF, &
                  JAC_RADAU, IJAC, MLJAC, MUJAC, &
                  MAS, IMAS, MLMAS, MUMAS, &
                  SOLOUT_RADAU, IOUT, &
                  WORKF, LWORKF, IWORKF, LIWORKF, RPARF, IPARF, IDID)

      DO J = 1, NEQF
         F(J, NT) = YF(J)
      END DO

      CALL CALCPEX(F(:, NT), P, CL1(NT), LHS1(NT), RHS1(NT), CL2(NT), LHS2(NT), RHS2(NT))
      ETA(:, NT) = EFF(P(:, NZ))
      !ETA(:, NT) = EFF(P(:, NZ))
      !ETAG(:, NT) = PITCH**2/(PITCH**2 + 1)*ETA(:, NT)
      XI1(NT) = XI(P, 1)
      XI2(NT) = XI(P, 2)
      XI1_KR(NT) = XI1(NT)/(F(1, NT)*CDEXP(IC*F(2, NT)))
      XI2_KR(NT) = XI2(NT)/(F(3, NT)*CDEXP(IC*F(4, NT)))
   END SUBROUTINE SOLVEF_RADAU

   SUBROUTINE SOLVEF_DOPR()
      USE, INTRINSIC :: ISO_C_BINDING
      USE DOPRF, ONLY: YF, IOUTF, AATOLF, ARPARF, ITOLF, WORKF, LWORKF, LIWORKF, AIPARF, IDIDF, ARTOLF, IWORKF
      IMPLICIT NONE

      REAL(C_DOUBLE) RPARF, T, XOUTF, RTOLF, ATOLF
      INTEGER(C_INT) IPARF, ITF, J

      REAL(C_DOUBLE) PEX(NEQP)

      COMMON/INTERNF/XOUTF, ITF
      !COMMON/INTERNF/ ITF

      RPARF = 0.0
      IPARF = 0
      T = TAX(1)
      XOUTF = T
      ITF = 0
      YF = F(:, 1)
      ITOLF = 0
      RTOLF = FTOL
      ATOLF = RTOLF
      IOUTF = 6

      ARTOLF(1) = RTOLF
      AATOLF(1) = ATOLF
      ARPARF(1) = RPARF
      AIPARF(1) = IPARF

      IWORKF(:) = 0
      WORKF(:) = 0.0D0

      IWORKF(5) = 6

      CALL DOPRI5_F(6, DFDT_DOPR, T, YF, TEND, ARTOLF, AATOLF, ITOLF, SOLOUTF, IOUTF, &
                    WORKF, LWORKF, IWORKF, LIWORKF, ARPARF, AIPARF, IDIDF)

      DO J = 1, NEQF
         F(J, NT) = YF(J)
      END DO

      CALL CALCPEX(F(:, NT), P, CL1(NT), LHS1(NT), RHS1(NT), CL2(NT), LHS2(NT), RHS2(NT))
      ETA(:, NT) = EFF(P(:, NZ))
      !ETA(:, NT) = EFF(P(:, NZ))
      !ETAG(:, NT) = PITCH**2/(PITCH**2 + 1)*ETA(:, NT)
      XI1(NT) = XI(P, 1)
      XI2(NT) = XI(P, 2)
      XI1_KR(NT) = XI1(NT)/(F(1, NT)*CDEXP(IC*F(2, NT)))
      XI2_KR(NT) = XI2(NT)/(F(3, NT)*CDEXP(IC*F(4, NT)))
   END SUBROUTINE SOLVEF_DOPR

   SUBROUTINE SOLVEP_DOPR(PIN, PEX, C)
      USE DOPRP, ONLY: YP, IWORKP, IOUTP, ARTOLP, AATOLP, ARPARP, ITOLP, WORKP, LWORKP, LIWORKP, AIPARP, IDIDP

      IMPLICIT NONE

      REAL(C_DOUBLE), INTENT(IN) :: PIN(:)
      REAL(C_DOUBLE), INTENT(INOUT) :: PEX(:, :)
      CHARACTER(C_CHAR), INTENT(IN) :: C

      INTEGER(C_INT) ITP
      REAL(C_DOUBLE) Z, XOUTP
      COMMON/INTERNP/XOUTP, ITP

      CALL INIT_DOPR_P()

      IF (C .EQ. 'P') THEN
         IOUTP = 2
         Z = ZAX(1)
         XOUTP = Z
         ITP = 0
         YP(:) = PIN(:)
         IWORKP(5) = NEQP

         CALL DOPRI5_P(NEQP, DPDZ_DOPR, Z, YP, ZEX, ARTOLP, AATOLP, ITOLP, SOLOUTP, IOUTP, &
                       WORKP, LWORKP, IWORKP, LIWORKP, ARPARP, AIPARP, IDIDP)

         PEX(:, NZ) = YP(:)

      ELSE
         IOUTP = 0
         Z = ZAX(1)
         XOUTP = Z
         ITP = 0
         YP(:) = PIN(:)
         IWORKP(5) = 0

         CALL DOPRI5_P(NEQP, DPDZ_DOPR, Z, YP, ZEX, ARTOLP, AATOLP, ITOLP, SOLOUT_FICTION, 0, &
                       WORKP, LWORKP, IWORKP, LIWORKP, ARPARP, AIPARP, IDIDP)

         PEX(:, 1) = YP(:)
      END IF

   END SUBROUTINE SOLVEP_DOPR

   SUBROUTINE SOLVEP_NAG(PIN, PEX, C)
      USE D02PVF_D02PCF_D02PDF, ONLY: ZSTART, THRES, METHOD, ERRASS, HSTART, WORKP, LENWRK, IFAILP, PGOT, PPGOT, PMAX

      IMPLICIT NONE

      REAL(C_DOUBLE), INTENT(IN) :: PIN(:)
      REAL(C_DOUBLE), INTENT(INOUT) :: PEX(:, :)
      CHARACTER(C_CHAR), INTENT(IN) :: C

      INTEGER(C_INT) I
      REAL(C_DOUBLE) :: ZWANT, ZGOT

      CALL D02PVF(NEQP, ZSTART, PIN(:), ZEX, PTOL, THRES, METHOD, 'USUAL TASK', ERRASS, HSTART, WORKP, LENWRK, IFAILP)

      IF (C .EQ. 'P') THEN
         DO I = 1, NZ - 1
            !ZWANT = I*DZ
            ZWANT = ZAX(I + 1)
            CALL D02PCF(DPDZ_NAG, ZWANT, ZGOT, PGOT, PPGOT, PMAX, WORKP, IFAILP)

            IF (IFAILP .NE. 0) THEN
               WRITE (*, *)
               WRITE (*, 99998) 'EXIT D02PCF WITH IFAIL = ', IFAILP, '  AND Z = ', ZWANT
               PAUSE
               STOP
            END IF

            P(:, I + 1) = PGOT
         END DO
      ELSE
         CALL D02PCF(DPDZ_NAG, ZEX, ZGOT, PEX, PPGOT, PMAX, WORKP, IFAILP)
      END IF
99998 FORMAT(1X, A, I2, A, D12.5)
   END SUBROUTINE SOLVEP_NAG

   SUBROUTINE SOLVEF_NAG()
      USE D02NVF_D02M_N, ONLY:RTOL, NEQMAX, NY2DIM, MAXORD, PETZLD, CONST, TCRIT, HMIN, HMAX, H0, MAXSTP, MXHNIL, NWKJAC, RWORK, IFAILF, NEQ, &
                               T, TOUT, ITASK, ITOL, ATOL, Y, ITRACE, XOUT, IT, YDOT, INFORM, YSAVE, WKJAC, NWKJAC, W_MONITR

      IMPLICIT NONE

      CALL D02NVF(NEQMAX, NY2DIM, MAXORD, 'DEFAULT', PETZLD, CONST, TCRIT, HMIN, HMAX, &
                  H0, MAXSTP, MXHNIL, 'DEFAULT', RWORK, IFAILF)
      CALL D02NSF(NEQ, NEQMAX, 'A', NWKJAC, RWORK, IFAILF)

      T = 0.0D0
      TOUT = TEND
      ITASK = 1
      ITOL = 1
      RTOL(1) = FTOL
      ATOL(1) = 0.001*FTOL
      DO I = 1, 6
         CONST(I) = 0.0D0
      END DO
      TCRIT = TOUT
      Y(:) = F(:, 1)
      ITRACE = -1
      XOUT = DT
      IT = 2
      IFAILF = 1

      !PRINT *, T, TOUT, ITASK, RTOL(1), ATOL(1), CONST, TCRIT, Y, ITRAICE, XOUT, IT, IFAIL
      !STOP

      CALL D02NBF(NEQ, NEQMAX, T, TOUT, Y, YDOT, RWORK, RTOL, ATOL, ITOL, INFORM, DFDT_NAG, YSAVE, &
                  NY2DIM, JAC_NAG, WKJAC, NWKJAC, MONITR, ITASK, ITRACE, IFAILF)

      IF (IFAILF .NE. 0) THEN
         WRITE (*, *)
         WRITE (*, 99998) 'EXIT D02NBF WITH IFAIL = ', IFAILF, '  AND T = ', T
         !PAUSE
         !STOP
      END IF
99998 FORMAT(1X, A, I2, A, D12.5)
   END SUBROUTINE SOLVEF_NAG

   SUBROUTINE ODE4F()
      !IMPORT, ONLY:FP, F, IC, P, ETA, ETAG, NZ, PITCH
      IMPLICIT NONE

      FP(1) = F(1, 1)*CDEXP(IC*F(2, 1))
      FP(2) = F(3, 1)*CDEXP(IC*F(4, 1))

      !SOLVE EQ. AT T=0
      IF (METH .EQ. 1) THEN
         CALL SOLVEP_DOPR(P(:, 1), P, 'P')
      ELSEIF (METH .EQ. 2) THEN
         CALL SOLVEP_NAG(P(:, 1), P, 'P')
      ELSEIF (METH .EQ. 3) THEN
         CALL SOLVEP_RK(P(:, 1), P, 'P')
      ELSEIF (METH .EQ. 4) THEN
         CALL SOLVEP_DOPR(P(:, 1), P, 'P')
      END IF

      !OPEN (1, FILE='TEST.DAT')
      !DO I = 1, NZ
      !   WRITE (1, '(F14.6,A,\)') ZAX(I), ' '
      !   DO J = 1, NE
      !      WRITE (1, '(F14.6,A,\)') CDABS(DCMPLX(P(J, I), P(NE + J, I))), ' '
      !   END DO
      !   WRITE (1, '(/\)')
      !END DO
      !CLOSE (1)
      !STOP

      ETA(:, 1) = EFF(P(:, NZ))
      !ETAG(:, 1) = PITCH**2/(PITCH**2 + 1.0D0)*ETA(:, 1)
      XI1(1) = XI(P, 1)
      XI2(1) = XI(P, 2)
      XI1_KR(1) = XI1(1)/FP(1)
      XI2_KR(1) = XI2(1)/FP(2)

      !TIME SOLVER
      IF (METH .EQ. 1) THEN
         CALL SOLVEF_DOPR()
      ELSEIF (METH .EQ. 2) THEN
         CALL SOLVEF_NAG()
      ELSEIF (METH .EQ. 3) THEN
         CALL SOLVEF_RK()
      ELSEIF (METH .EQ. 4) THEN
         CALL SOLVEF_RADAU()
      END IF
   END SUBROUTINE ODE4F

   FUNCTION EFF(PEX) RESULT(ETA)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE, C_INT
      IMPORT, ONLY:NE, IDXRE, IDXIM

      IMPLICIT NONE

      INTEGER(C_INT) I
      REAL(C_DOUBLE) ETA(2)
      REAL(C_DOUBLE), INTENT(IN) :: PEX(:)

      DO I = 1, 2
         ETA(I) = 1 - SUM(CDABS(DCMPLX(PEX(IDXRE(I, :)), PEX(IDXIM(I, :))))**2)/NE
      END DO
   END FUNCTION EFF

   SUBROUTINE DPDZ_COMMON(Z, P, PRHS)
      IMPORT :: NE, ZEX, F, IC, DTR
      IMPLICIT NONE

      REAL(C_DOUBLE) Z, P(*), PRHS(*)

      INTEGER(C_INT) I, REIDX(NE), IMIDX(NE)
      COMPLEX(C_DOUBLE_COMPLEX) :: U
      COMPLEX(C_DOUBLE_COMPLEX) S(NE), PTMP(NE)

      IF (FOK .EQ. .FALSE.) THEN
         U = DEXP(-3.0D0*((Z - ZEX/2)/(ZEX/2))**2)
      ELSEIF (SQR .EQ. .TRUE.) THEN
         IF (INHARM .EQ. 1) THEN
            U = SQUVAL22(Z)
         ELSE IF (INHARM .EQ. 2) THEN
            U = SQUVAL(Z)
         END IF
      ELSE
         IF (INHARM .EQ. 1) THEN
            !U = UVAL22(Z)
            U = CMPLXEVAL(918, Z, Z_ANG, RE_ANG, IM_ANG, REB, REC, RED, IMB, IMC, IMD)            
         ELSE IF (INHARM .EQ. 2) THEN
            U = UVAL(Z)
         END IF
      END IF

      DO I = 1, 2
         PTMP = DCMPLX(P(IDXRE(I, :)), P(IDXIM(I, :)))

         S = IC*(FP(I)*U*DCONJG(PTMP)**(INHARM - 1) - (DTR(I) + CDABS(PTMP)**2 - 1)*PTMP)

         PRHS(IDXRE(I, :)) = DREAL(S)
         PRHS(IDXIM(I, :)) = DIMAG(S)
      END DO
   END SUBROUTINE DPDZ_COMMON

   SUBROUTINE DPDZ_DOPR(NEQP, Z, P, PRHS, RPARP, IPARP)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_DOUBLE_COMPLEX

      IMPLICIT NONE

      REAL(C_DOUBLE) Z, P(*), PRHS(*)
      INTEGER(C_INT) NEQP, IPARP
      REAL(C_DOUBLE) RPARP

      CALL DPDZ_COMMON(Z, P, PRHS)
   END SUBROUTINE DPDZ_DOPR

   SUBROUTINE DPDZ_NAG(Z, P, PRHS)
      IMPLICIT NONE

      REAL(C_DOUBLE) Z, P(*), PRHS(*)

      CALL DPDZ_COMMON(Z, P, PRHS)
   END SUBROUTINE DPDZ_NAG

   COMPLEX(C_DOUBLE_COMPLEX) FUNCTION XI(P, NUM)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE, C_DOUBLE_COMPLEX
      IMPORT, ONLY:NE, NZ, MEAN, U, DZ, IDXRE, IDXIM, INHARM

      IMPLICIT NONE

      INTEGER(C_INT) I, NUM, J
      REAL(C_DOUBLE) P(:, :)

      DO I = 1, NZ
         MEAN(I) = SUM(DCMPLX(P(IDXRE(NUM, :), I), P(IDXIM(NUM, :), I))**INHARM, 1)/NE
      END DO

      !OPEN (10, FILE='TEST.DAT')
      !DO I = 1, NZ
      !   WRITE(10,'(F12.6,A,\)') (I - 1)*DZ, '   '
      !   DO J = 1, NE
      !      WRITE (10, '(F12.6,A,\)') CDABS(DCMPLX(P(J, I), P(NE + J, I))), '    '
      !   END DO
      !   WRITE(10,'(/\)')
      !END DO
      !CLOSE (10)
      !STOP

      MEAN = DCONJG(U)*MEAN

      XI = (0.5D0*(MEAN(1) + MEAN(NZ)) + SUM(MEAN(2:(NZ - 1))))*DZ

   END FUNCTION

   SUBROUTINE DFDT_NAG(NEQF, T, F, S, IRES)
      IMPLICIT NONE

      INTEGER(C_INT) NEQF, IRES
      REAL(C_DOUBLE) T, F(NEQF), S(NEQF)

      FP(1) = F(1)*CDEXP(IC*F(2))
      FP(2) = F(3)*CDEXP(IC*F(4))

      !IF (T == 0.5) THEN
      !   OPEN (1, FILE='TEST.DAT')
      !   DO I = 1, NZ
      !      WRITE (1, '(F14.6,A,\)') ZAX(I), ' '
      !      DO J = 1, NE
      !         WRITE (1, '(F14.6,A,\)') CDABS(DCMPLX(P(J, I), P(NE + J, I))), ' '
      !      END DO
      !      WRITE (1, '(/\)')
      !   END DO
      !   CLOSE (1)
      !   STOP
      !END IF

      CALL SOLVEP_NAG(P(:, 1), P, 'P')

      CALL DFDT_COMMON(NEQF, T, F, S)

   END SUBROUTINE DFDT_NAG

   SUBROUTINE DFDT_RADAU(NEQF, T, F, S, RPARF, IPARF)
      IMPLICIT NONE
      INTEGER(C_INT) :: NEQF, IPARF, ITP
      REAL(C_DOUBLE) T, F(NEQF), S(NEQF), RPARF

      FP(1) = F(1)*CDEXP(IC*F(2))
      FP(2) = F(3)*CDEXP(IC*F(4))

      CALL SOLVEP_DOPR(P(:, 1), P, 'P')

      CALL DFDT_COMMON(NEQF, T, F, S)

   END SUBROUTINE DFDT_RADAU

   SUBROUTINE DFDT_DOPR(NEQF, T, F, S, RPARF, IPARF)
      IMPLICIT NONE
      INTEGER(C_INT) :: NEQF, IPARF, ITP
      REAL(C_DOUBLE) T, F(NEQF), S(NEQF), RPARF

      FP(1) = F(1)*CDEXP(IC*F(2))
      FP(2) = F(3)*CDEXP(IC*F(4))

      CALL SOLVEP_DOPR(P(:, 1), P, 'P')

      !IF (T == 0.5) THEN
      !   OPEN (1, FILE='TEST.DAT')
      !   DO I = 1, NZ
      !      WRITE (1, '(F14.6,A,\)') ZAX(I), ' '
      !      DO J = 1, NE
      !         WRITE (1, '(F14.6,A,\)') CDABS(DCMPLX(P(J, I), P(NE + J, I))), ' '
      !      END DO
      !      WRITE (1, '(/\)')
      !   END DO
      !   CLOSE (1)
      !   STOP
      !END IF

      CALL DFDT_COMMON(NEQF, T, F, S)

   END SUBROUTINE DFDT_DOPR

   SUBROUTINE DFDT_COMMON(NEQF, T, F, S)
      IMPLICIT NONE

      INTEGER(C_INT) :: II, JJ, ITER_NUM = 1, TIME_NUM = 1, NEQF
      REAL(C_DOUBLE) T, F(NEQF), S(NEQF), &
         X1R, X1I, Q31, I1, R1, TH1, DCIR1, DCOS1, DSIN1, &
         X2R, X2I, Q32, I2, R2, TH2, DCIR2, DCOS2, DSIN2, Q3, &
         F1, F2, F3, PHI1, PHI2, PHI3, A1, A2, ZWANT, ZGOT, E1, E2
      COMPLEX(C_DOUBLE_COMPLEX) X1, X2
      LOGICAL BP

      X1 = XI(P, 1)
      X2 = XI(P, 2)

      X1R = DREAL(X1)
      X1I = DIMAG(X1)
      X2R = DREAL(X2)
      X2I = DIMAG(X2)

      F1 = F(1)
      PHI1 = F(2)
      F2 = F(3)
      PHI2 = F(4)
      F3 = F(5)
      PHI3 = F(6)

      Q31 = Q(3)/Q(1)
      I1 = ICU(1)
      R1 = R(1)
      TH1 = TH(1)
      DCIR1 = DCIR(1)
      DCOS1 = DCOS(PHI1)
      DSIN1 = DSIN(PHI1)

      Q32 = Q(3)/Q(2)
      I2 = ICU(2)
      R2 = R(2)
      TH2 = TH(2)
      DCIR2 = DCIR(2)
      DCOS2 = DCOS(PHI2)
      DSIN2 = DSIN(PHI2)

      Q3 = Q(3)
      A1 = A(1)
      A2 = A(2)

      S(1) = (-NHARM*F1 + I1*(-X1I*DCOS1 + X1R*DSIN1) + 2*R1*NHARM*F3*DCOS(PHI3 - PHI1 - TH1))*Q31
      S(2) = -2*DCIR1*Q3 + (I1/F1*(X1R*DCOS1 + X1I*DSIN1) + 2*R1*NHARM*(F3/F1)*DSIN(PHI3 - PHI1 - TH1))*Q31

      S(3) = (-NHARM*F2 + I2*(-X2I*DCOS2 + X2R*DSIN2) + 2*R2*NHARM*F3*DCOS(PHI3 - PHI2 - TH2))*Q32
      S(4) = -2*DCIR2*Q3 + (I2/F2*(X2R*DCOS2 + X2I*DSIN2) + 2*R2*NHARM*(F3/F2)*DSIN(PHI3 - PHI2 - TH2))*Q32

      S(5) = -F3 + A1*F1*DCOS(PHI1 - PHI3) + A2*F2*DCOS(PHI2 - PHI3)
      S(6) = A1*F1/F3*DSIN(PHI1 - PHI3) + A2*F2/F3*DSIN(PHI2 - PHI3)

   END SUBROUTINE DFDT_COMMON

   SUBROUTINE CALC_U(U, ZEX_W, NZ, ZAX)
      IMPORT
      IMPLICIT NONE

      INTEGER(C_INT), INTENT(IN) :: NZ
      REAL(C_DOUBLE), INTENT(IN) :: ZEX_W, ZAX(NZ)
      COMPLEX(C_DOUBLE_COMPLEX), INTENT(OUT) :: U(:)

      INTEGER(C_INT) I

      IF (FOK .EQ. .FALSE.) THEN
         DO I = 1, NZ
            U(I) = DEXP(-3.0D0*((ZAX(I) - ZEX_W/2)/(ZEX_W/2))**2)
         END DO
      ELSEIF (SQR .EQ. .TRUE.) THEN
         IF (INHARM .EQ. 1) THEN
            DO I = 1, NZ
               U(I) = SQUVAL22(ZAX(I))
            END DO
         ELSE IF (INHARM .EQ. 2) THEN
            DO I = 1, NZ
               U(I) = SQUVAL(ZAX(I))
            END DO
         END IF
      ELSE
         IF (INHARM .EQ. 1) THEN
            DO I = 1, NZ
               !U(I) = UVAL22(ZAX(I))
               U(I) = CMPLXEVAL(918, ZAX(I), Z_ANG, RE_ANG, IM_ANG, REB, REC, RED, IMB, IMC, IMD)             
            END DO
         ELSE IF (INHARM .EQ. 2) THEN
            DO I = 1, NZ
               U(I) = UVAL(ZAX(I))
            END DO
         END IF
      END IF

      OPEN (1, FILE='STRUC.DAT')
      DO I = 1, NZ
         IF (INHARM == 2) THEN
            WRITE (1, '(3F14.6)') ZAX(I)/ZAX(NZ)*185.5 - 8.5, DREAL(U(I)), DIMAG(U(I))
         ELSE
            !WRITE (1, '(3F14.6)') ZAX(I)/ZAX(NZ)*52.43, DREAL(U(I)), DIMAG(U(I))
            WRITE (1, '(3F14.6)') ZAX(I), DREAL(U(I)), DIMAG(U(I))
         END IF
      END DO
      CLOSE (1)
      !STOP

   END SUBROUTINE

   SUBROUTINE JAC_RADAU(N, X, Y, DFY, LDFY, RPAR, IPAR)
      DOUBLE PRECISION X, Y(N), DFY(LDFY, N)

      REAL(C_DOUBLE) X1R, X1I, Q31, I1, R1, TH1, DCIR1, DCOS1, DSIN1, &
         X2R, X2I, Q32, I2, R2, TH2, DCIR2, DCOS2, DSIN2, Q333, &
         F1, F2, F3, PHI1, PHI2, PHI3, A1, A2, PH311, PH322, PHI13, PHI23, &
         DSIN311, DSIN322, DCOS311, DCOS322, DCOS13, DSIN13, DCOS23, DSIN23, PHI131, PHI232
      COMPLEX(C_DOUBLE_COMPLEX) X1, X2

      FP(1) = Y(1)*CDEXP(IC*Y(2))
      FP(2) = Y(3)*CDEXP(IC*Y(4))

      CALL SOLVEP_DOPR(P(:, 1), P, 'P')

      X1 = XI(P, 1)
      X2 = XI(P, 2)

      X1R = DREAL(X1)
      X1I = DIMAG(X1)
      X2R = DREAL(X2)
      X2I = DIMAG(X2)

      F1 = Y(1)
      PHI1 = Y(2)
      F2 = Y(3)
      PHI2 = Y(4)
      F3 = Y(5)
      PHI3 = Y(6)

      Q31 = Q(3)/Q(1)
      I1 = ICU(1)
      R1 = R(1)
      TH1 = TH(1)
      DCIR1 = DCIR(1)
      DCOS1 = DCOS(PHI1)
      DSIN1 = DSIN(PHI1)

      Q32 = Q(3)/Q(2)
      I2 = ICU(2)
      R2 = R(2)
      TH2 = TH(2)
      DCIR2 = DCIR(2)
      DCOS2 = DCOS(PHI2)
      DSIN2 = DSIN(PHI2)

      PH311 = PHI3 - PHI1 - TH1
      PH322 = PHI3 - PHI2 - TH2
      PHI13 = PHI1 - PHI3
      PHI23 = PHI2 - PHI3
      PHI131 = PHI1 - PHI3 + TH1
      PHI232 = PHI2 - PHI3 + TH2

      DCOS311 = DCOS(PH311)
      DSIN311 = DSIN(PH311)
      DCOS322 = DCOS(PH322)
      DSIN322 = DSIN(PH322)
      DCOS13 = DCOS(PHI13)
      DSIN13 = DSIN(PHI13)
      DCOS23 = DCOS(PHI23)
      DSIN23 = DSIN(PHI23)

      DFY(1, 1) = -NHARM*Q31
      DFY(1, 2) = (I1*Q31*X1R*DCOS(PHI1) + I1*Q31*X1I*DSIN(PHI1) - 2*F3*NHARM*Q31*R1*DSIN(PHI1 - PHI3 + TH1))
      DFY(1, 3) = 0
      DFY(1, 4) = 0
      DFY(1, 5) = 2*NHARM*Q31*R1*DCOS(PHI1 - PHI3 + TH1)
      DFY(1, 6) = 2*F3*NHARM*Q31*R1*DSIN(PHI1 - PHI3 + TH1)

    DFY(2, 1) = -((I1*Q31*X1R*DCOS(PHI1))/F1**2) - (I1*Q31*X1I*DSIN(PHI1))/F1**2 + (2*F3*NHARM*Q31*R1*DSIN(PHI1 - PHI3 + TH1))/F1**2
      DFY(2, 2) = (I1*Q31*X1I*DCOS(PHI1))/F1 - (2*F3*NHARM*Q31*R1*DCOS(PHI1 - PHI3 + TH1))/F1 - (I1*Q31*X1R*DSIN(PHI1))/F1
      DFY(2, 3) = 0
      DFY(2, 4) = 0
      DFY(2, 5) = -2*NHARM*Q31*R1*DSIN(PHI1 - PHI3 + TH1)/F1
      DFY(2, 6) = 2*F3*NHARM*Q31*R1*DCOS(PHI1 - PHI3 + TH1)/F1

      DFY(3, 1) = 0
      DFY(3, 2) = 0
      DFY(3, 3) = -NHARM*Q32
      DFY(3, 4) = (I2*Q32*X2R*DCOS(PHI2) + I2*Q32*X2I*DSIN(PHI2) - 2*F3*NHARM*Q32*R2*DSIN(PHI2 - PHI3 + TH2))
      DFY(3, 5) = 2*NHARM*Q32*R2*DCOS(PHI2 - PHI3 + TH2)
      DFY(3, 6) = 2*F3*NHARM*Q32*R2*DSIN(PHI2 - PHI3 + TH2)

      DFY(4, 1) = 0
      DFY(4, 2) = 0
    DFY(4, 3) = -((I2*Q32*X2R*DCOS(PHI2))/F2**2) - (I2*Q32*X2I*DSIN(PHI2))/F2**2 + (2*F3*NHARM*Q32*R2*DSIN(PHI2 - PHI3 + TH2))/F2**2
      DFY(4, 4) = (I2*Q32*X2I*DCOS(PHI2))/F2 - (2*F3*NHARM*Q32*R2*DCOS(PHI2 - PHI3 + TH2))/F2 - (I2*Q32*X2R*DSIN(PHI2))/F2
      DFY(4, 5) = -2*NHARM*Q32*R2*DSIN(PHI2 - PHI3 + TH2)/F2
      DFY(4, 6) = 2*F3*NHARM*Q32*R2*DCOS(PHI2 - PHI3 + TH2)/F2

      DFY(5, 1) = (A1*DCOS(PHI1 - PHI3))
      DFY(5, 2) = -A1*F1*DSIN(PHI1 - PHI3)
      DFY(5, 3) = (A2*DCOS(PHI2 - PHI3))
      DFY(5, 4) = -A2*F2*DSIN(PHI2 - PHI3)
      DFY(5, 5) = -1
      DFY(5, 6) = ((A1*F1*DSIN(PHI1 - PHI3) + A2*F2*DSIN(PHI2 - PHI3)))

      DFY(6, 1) = A1*DSIN(PHI1 - PHI3)/F3
      DFY(6, 2) = A1*F1*DCOS(PHI1 - PHI3)/F3
      DFY(6, 3) = A2*DSIN(PHI2 - PHI3)/F3
      DFY(6, 4) = A2*F2*DCOS(PHI2 - PHI3)/F3
      DFY(6, 5) = -(A1*F1*DSIN(PHI1 - PHI3)/F3**2) - A2*F2*DSIN(PHI2 - PHI3)/F3**2
      DFY(6, 6) = -(A1*F1*DCOS(PHI1 - PHI3)/F3) - A2*F2*DCOS(PHI2 - PHI3)/F3
   END SUBROUTINE JAC_RADAU

   SUBROUTINE MAS(N, AM, LMAS, RPAR, IPAR)
!C                    DOUBLE PRECISION AM(LMAS,N)
   END SUBROUTINE MAS

   SUBROUTINE JAC_NAG(NEQ, T, Y, H, D, PP)
      IMPLICIT NONE
!     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION D, H, T
      INTEGER NEQ
!     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION PP(NEQ, NEQ), Y(NEQ)
!     .. LOCAL SCALARS ..
      DOUBLE PRECISION HXD
!     .. EXECUTABLE STATEMENTS ..

      REAL(C_DOUBLE) X1R, X1I, Q31, I1, R1, TH1, DCIR1, DCOS1, DSIN1, &
         X2R, X2I, Q32, I2, R2, TH2, DCIR2, DCOS2, DSIN2, Q333, &
         F1, F2, F3, PHI1, PHI2, PHI3, A1, A2, PH311, PH322, PHI13, PHI23, &
         DSIN311, DSIN322, DCOS311, DCOS322, DCOS13, DSIN13, DCOS23, DSIN23, PHI131, PHI232
      COMPLEX(C_DOUBLE_COMPLEX) X1, X2

      HXD = H*D

      FP(1) = Y(1)*CDEXP(IC*Y(2))
      FP(2) = Y(3)*CDEXP(IC*Y(4))

      CALL SOLVEP_NAG(P(:, 1), P, 'P')

      X1 = XI(P, 1)
      X2 = XI(P, 2)

      X1R = DREAL(X1)
      X1I = DIMAG(X1)
      X2R = DREAL(X2)
      X2I = DIMAG(X2)

      F1 = Y(1)
      PHI1 = Y(2)
      F2 = Y(3)
      PHI2 = Y(4)
      F3 = Y(5)
      PHI3 = Y(6)

      Q31 = Q(3)/Q(1)
      I1 = ICU(1)
      R1 = R(1)
      TH1 = TH(1)
      DCIR1 = DCIR(1)
      DCOS1 = DCOS(PHI1)
      DSIN1 = DSIN(PHI1)

      Q32 = Q(3)/Q(2)
      I2 = ICU(2)
      R2 = R(2)
      TH2 = TH(2)
      DCIR2 = DCIR(2)
      DCOS2 = DCOS(PHI2)
      DSIN2 = DSIN(PHI2)

      PH311 = PHI3 - PHI1 - TH1
      PH322 = PHI3 - PHI2 - TH2
      PHI13 = PHI1 - PHI3
      PHI23 = PHI2 - PHI3
      PHI131 = PHI1 - PHI3 + TH1
      PHI232 = PHI2 - PHI3 + TH2

      DCOS311 = DCOS(PH311)
      DSIN311 = DSIN(PH311)
      DCOS322 = DCOS(PH322)
      DSIN322 = DSIN(PH322)
      DCOS13 = DCOS(PHI13)
      DSIN13 = DSIN(PHI13)
      DCOS23 = DCOS(PHI23)
      DSIN23 = DSIN(PHI23)

      PP(1, 1) = 1 + HXD*NHARM*Q31
      PP(1, 2) = -(HXD*(I1*Q31*X1R*DCOS(PHI1) + I1*Q31*X1I*DSIN(PHI1) - 2*F3*NHARM*Q31*R1*DSIN(PHI1 - PHI3 + TH1)))
      !PP(1,3) = 0
      !PP(1,4) = 0
      PP(1, 5) = -2*HXD*NHARM*Q31*R1*DCOS(PHI1 - PHI3 + TH1)
      PP(1, 6) = -2*F3*HXD*NHARM*Q31*R1*DSIN(PHI1 - PHI3 + TH1)

      PP(2, 1) = -(HXD*(-((I1*Q31*X1R*DCOS(PHI1))/F1**2) - (I1*Q31*X1I*DSIN(PHI1))/F1**2 + (2*F3*NHARM*Q31*R1*DSIN(PHI1 - PHI3 + TH1))/F1**2))
      PP(2, 2) = 1 - HXD*((I1*Q31*X1I*DCOS(PHI1))/F1 - (2*F3*NHARM*Q31*R1*DCOS(PHI1 - PHI3 + TH1))/F1 - (I1*Q31*X1R*DSIN(PHI1))/F1)
      !PP(2,3) = 0
      !PP(2,4) = 0
      PP(2, 5) = (2*HXD*NHARM*Q31*R1*DSIN(PHI1 - PHI3 + TH1))/F1
      PP(2, 6) = (-2*F3*HXD*NHARM*Q31*R1*DCOS(PHI1 - PHI3 + TH1))/F1

      !PP(3,1) = 0
      !PP(3,2) = 0
      PP(3, 3) = 1 + HXD*NHARM*Q32
      PP(3, 4) = -(HXD*(I2*Q32*X2R*DCOS(PHI2) + I2*Q32*X2I*DSIN(PHI2) - 2*F3*NHARM*Q32*R2*DSIN(PHI2 - PHI3 + TH2)))
      PP(3, 5) = -2*HXD*NHARM*Q32*R2*DCOS(PHI2 - PHI3 + TH2)
      PP(3, 6) = -2*F3*HXD*NHARM*Q32*R2*DSIN(PHI2 - PHI3 + TH2)

      !PP(4,1) = 0
      !PP(4,2) = 0
      PP(4, 3) = -(HXD*(-((I2*Q32*X2R*DCOS(PHI2))/F2**2) - (I2*Q32*X2I*DSIN(PHI2))/F2**2 + (2*F3*NHARM*Q32*R2*DSIN(PHI2 - PHI3 + TH2))/F2**2))
      PP(4, 4) = 1 - HXD*((I2*Q32*X2I*DCOS(PHI2))/F2 - (2*F3*NHARM*Q32*R2*DCOS(PHI2 - PHI3 + TH2))/F2 - (I2*Q32*X2R*DSIN(PHI2))/F2)
      PP(4, 5) = (2*HXD*NHARM*Q32*R2*DSIN(PHI2 - PHI3 + TH2))/F2
      PP(4, 6) = (-2*F3*HXD*NHARM*Q32*R2*DCOS(PHI2 - PHI3 + TH2))/F2

      PP(5, 1) = -(A1*HXD*DCOS(PHI1 - PHI3))
      PP(5, 2) = A1*F1*HXD*DSIN(PHI1 - PHI3)
      PP(5, 3) = -(A2*HXD*DCOS(PHI2 - PHI3))
      PP(5, 4) = A2*F2*HXD*DSIN(PHI2 - PHI3)
      PP(5, 5) = 1
      PP(5, 6) = -(HXD*(A1*F1*DSIN(PHI1 - PHI3) + A2*F2*DSIN(PHI2 - PHI3)))

      PP(6, 1) = -((A1*HXD*DSIN(PHI1 - PHI3))/F3)
      PP(6, 2) = -((A1*F1*HXD*DCOS(PHI1 - PHI3))/F3)
      PP(6, 3) = -((A2*HXD*DSIN(PHI2 - PHI3))/F3)
      PP(6, 4) = -((A2*F2*HXD*DCOS(PHI2 - PHI3))/F3)
      PP(6, 5) = -(HXD*(-((A1*F1*DSIN(PHI1 - PHI3))/F3**2) - (A2*F2*DSIN(PHI2 - PHI3))/F3**2))
      PP(6, 6) = 1 - HXD*(-((A1*F1*DCOS(PHI1 - PHI3))/F3) - (A2*F2*DCOS(PHI2 - PHI3))/F3)
      RETURN
   END

   SUBROUTINE MONITR(N, NMAX, T, HLAST, H, Y, YDOT, YSAVE, R, ACOR, IMON, INLN, HMIN, HMXI, NQU)
      USE D02NVF_D02M_N, ONLY: W_MONITR
      IMPLICIT NONE
      !..PARAMETERS..
      INTEGER NOUT, IT, J
      PARAMETER(NOUT=6)
      INTEGER NY2DIM
      PARAMETER(NY2DIM=6)
      !..SCALAR ARGUMENTS..
      DOUBLE PRECISION H, HLAST, HMIN, HMXI, T
      INTEGER IMON, INLN, N, NMAX, NQU
      !..ARRAY ARGUMENTS..
      DOUBLE PRECISION ACOR(NMAX, 2), R(N), Y(N), YDOT(N), YSAVE(NMAX, *)
      !..SCALARS IN COMMON..
      DOUBLE PRECISION XOUT
      !..LOCAL SCALARS..
      INTEGER I, IFAIL
      LOGICAL(4) PRESSED
      CHARACTER(1) KEY
      INTEGER(C_INT), PARAMETER :: ESC = 27

      REAL(C_DOUBLE) PEX(NEQP)

      !..EXTERNAL SUBROUTINES..
      EXTERNAL D02XKF
      !..COMMON BLOCKS..
      COMMON XOUT, IT
      !..EXECUTABLE STATEMENTS..
      IF (IMON .NE. 1) RETURN
20    IF (.NOT. (T - HLAST .LT. XOUT .AND. XOUT .LE. T)) RETURN
      IFAIL = 1
      !C1 INTERPOLATION
      CALL D02XKF(XOUT, R, N, YSAVE, NMAX, NY2DIM, ACOR(1, 2), N, T, NQU, HLAST, H, IFAIL)

      IF (IFAIL .NE. 0) THEN
         IMON = -2
      ELSE
         !WRITE (NOUT, 99999) XOUT, (R(I), I=1, N)
         CALL CALCPEX(R, P, CL1(IT), LHS1(IT), RHS1(IT), CL2(IT), LHS2(IT), RHS2(IT))
         ETA(:, IT) = EFF(P(:, NZ))
         !ETAG(:, IT) = PITCH**2/(PITCH**2 + 1)*ETA(:, IT)
         WRITE (*, '(A,F8.3,A,F8.5,A,F8.5,A,F8.5,A,F8.5,A,F8.5,A,F9.5,A,F9.5,A,F9.5,A,F5.3,A,F5.3,A,\,A)') 'T =', XOUT, &
            '  |F1| =', R(1), '  |F2| =', R(3), '  |F3| =', R(5), '  E1 =', ETA(1, IT), '  E2 =', ETA(2, IT), &
            '  W1 = ', YDOT(2), '  W2 = ', YDOT(4), '  W3 = ', YDOT(6), '  C1 = ', DABS(CL1(IT)/RHS1(IT))*100, '%  C2 = ', DABS(CL2(IT)/RHS2(IT)*100), '%', CHAR(13)
         DO J = 1, N
            F(J, IT) = R(J)
         END DO
         DO J = 1, 3
            W_MONITR(J, IT) = YDOT(2*J)
            W(J, IT) = W_MONITR(J, IT)
         END DO

         XI1(IT) = XI(P, 1)
         XI2(IT) = XI(P, 2)
         XI1_KR(IT) = XI1(IT)/(F(1, IT)*CDEXP(IC*F(2, IT)))
         XI2_KR(IT) = XI2(IT)/(F(3, IT)*CDEXP(IC*F(4, IT)))

         XOUT = XOUT + DT
         NT = IT
         NTEND = NT
         IT = IT + 1
         PRESSED = PEEKCHARQQ()
         IF (PRESSED) THEN
            KEY = GETCHARQQ()
            IF (ICHAR(KEY) .EQ. ESC) THEN
               WRITE (*, '(/,A)') 'QUIT?'
               KEY = GETCHARQQ()
               IF (ICHAR(KEY) .EQ. 121 .OR. ICHAR(KEY) .EQ. 89) THEN
                  NT = IT - 1
                  NTEND = NT
                  CALL WRITE_RESULTS()
                  !IMON = -2
                  !RETURN
               END IF
            END IF
         END IF
         IF (XOUT .LT. TEND) GO TO 20
      END IF

      RETURN

99999 FORMAT(1X, F8.3, 6(F13.5, 2X))
   END SUBROUTINE MONITR

   SUBROUTINE CALCPEX(FPEX, PEX, C1, LHS1, RHS1, C2, LHS2, RHS2)

      IMPLICIT NONE

      REAL(8), INTENT(OUT) :: C1, LHS1, RHS1, C2, LHS2, RHS2

      INTEGER(C_INT) I, IDX(NE)
      REAL(C_DOUBLE) :: ZWANT, ZGOT, PEX(:, :), FPEX(6), P2EX_MEAN(2), P20_MEAN(2)
      COMPLEX(8) PTMP(NE)

      FP(1) = FPEX(1)*CDEXP(IC*FPEX(2))
      FP(2) = FPEX(3)*CDEXP(IC*FPEX(4))

      IF (METH .EQ. 1) THEN
         CALL SOLVEP_DOPR(P(:, 1), PEX, 'P')
      ELSEIF (METH .EQ. 2) THEN
         CALL SOLVEP_NAG(P(:, 1), PEX, 'P')
      ELSEIF (METH .EQ. 3) THEN
         CALL SOLVEP_DOPR(P(:, 1), PEX, 'P')
      ELSEIF (METH .EQ. 4) THEN
         CALL SOLVEP_DOPR(P(:, 1), PEX, 'P')
      END IF

      DO I = 1, 2
         IDX = IDXP(I, :)

         PTMP(:) = DCMPLX(P(IDXRE(I, :), 1), P(IDXIM(I, :), 1))
         P20_MEAN(I) = SUM(CDABS(PTMP(:)*CDABS(PTMP(:))))/NE

         PTMP(:) = DCMPLX(PEX(IDXRE(I, :), NZ), PEX(IDXIM(I, :), NZ))
         P2EX_MEAN(I) = SUM(CDABS(PTMP(:)*CDABS(PTMP(:))))/NE
      END DO

      LHS1 = 2.0*NHARM*FPEX(1)**2 - 4.0*NHARM*R(1)*FPEX(5)*FPEX(1)*DCOS(TH(1) - FPEX(6) + FPEX(2))
      RHS1 = -ICU(1)*(P2EX_MEAN(1) - P20_MEAN(1))
      C1 = LHS1 - RHS1

      LHS2 = 2.0*NHARM*FPEX(3)**2 - 4.0*NHARM*R(2)*FPEX(5)*FPEX(3)*DCOS(TH(2) - FPEX(6) + FPEX(4))
      RHS2 = -ICU(2)*(P2EX_MEAN(2) - P20_MEAN(2))
      C2 = LHS2 - RHS2

   END SUBROUTINE CALCPEX

   SUBROUTINE FREQ()

      IMPLICIT NONE

      INTEGER(C_INT) :: II, IPARF, AIPARP(1), ITP, J
      REAL(C_DOUBLE) T, Z, ZWANT, ZGOT, &
         X1R, X1I, Q31, I1, R1, TH1, DCIR1, DCOS1, DSIN1, &
         X2R, X2I, Q32, I2, R2, TH2, DCIR2, DCOS2, DSIN2, Q3, &
         F1, F2, F3, PHI1, PHI2, PHI3, A1, A2
      COMPLEX(C_DOUBLE_COMPLEX) X1, X2

      DO II = 2, NT

         FP(1) = F(1, II)*EXP(IC*F(2, II))
         FP(2) = F(3, II)*EXP(IC*F(4, II))

         IF (METH .EQ. 1) THEN
            CALL SOLVEP_DOPR(P(:, 1), P, 'P')
         ELSEIF (METH .EQ. 2) THEN
            CALL SOLVEP_NAG(P(:, 1), P, 'P')
         ELSEIF (METH .EQ. 3) THEN
            CALL SOLVEP_RK(P(:, 1), P, 'P')
         ELSEIF (METH .EQ. 4) THEN
            CALL SOLVEP_RK(P(:, 1), P, 'P')
         END IF

         !X1 = XI(P(1:2*NE, :), 1)
         !X2 = XI(P(2*NE + 1:4*NE, :), 1)
         X1 = XI(P, 1)
         X2 = XI(P, 2)

         X1R = DREAL(X1)
         X1I = DIMAG(X1)
         X2R = DREAL(X2)
         X2I = DIMAG(X2)

         F1 = F(1, II)
         PHI1 = F(2, II)
         F2 = F(3, II)
         PHI2 = F(4, II)
         F3 = F(5, II)
         PHI3 = F(6, II)

         Q31 = Q(3)/Q(1)
         I1 = ICU(1)
         R1 = R(1)
         TH1 = TH(1)
         DCIR1 = DCIR(1)
         DCOS1 = DCOS(PHI1)
         DSIN1 = DSIN(PHI1)

         Q32 = Q(3)/Q(2)
         I2 = ICU(2)
         R2 = R(2)
         TH2 = TH(2)
         DCIR2 = DCIR(2)
         DCOS2 = DCOS(PHI2)
         DSIN2 = DSIN(PHI2)

         Q3 = Q(3)
         A1 = A(1)
         A2 = A(2)

         !S(1) = (-NHARM*F1 + I1*(-X1I*DCOS1 + X1R*DSIN1) + 2*R1*NHARM*F3*DCOS(PHI3 - PHI1 - TH1))*Q31
         W(1, II) = -2*DCIR1*Q3 + (I1/F1*(X1R*DCOS1 + X1I*DSIN1) + 2*R1*NHARM*(F3/F1)*DSIN(PHI3 - PHI1 - TH1))*Q31

         !S(3) = (-NHARM*F2 + I2*(-X2I*DCOS2 + X2R*DSIN2) + 2*R2*NHARM*F3*DCOS(PHI3 - PHI2 - TH2))*Q32
         W(2, II) = -2*DCIR2*Q3 + (I2/F2*(X2R*DCOS2 + X2I*DSIN2) + 2*R2*NHARM*(F3/F2)*DSIN(PHI3 - PHI2 - TH2))*Q32

         !S(5) = -F3 + A1*F1*DCOS(PHI1 - PHI3) + A2*F2*DCOS(PHI2 - PHI3)
         W(3, II) = A1*F1/F3*DSIN(PHI1 - PHI3) + A2*F2/F3*DSIN(PHI2 - PHI3)
      END DO
   END SUBROUTINE FREQ

   SUBROUTINE SOLOUT_FICTION
   END SUBROUTINE SOLOUT_FICTION

   SUBROUTINE SOLOUTP(NR, XOLD, X, Y, N, CON, ICOMP, ND, RPARP, IPARP, IRTRN)
      IMPLICIT NONE

      INTERFACE
         FUNCTION CONTD5_P(II, X, CON, ICOMP, ND)
            IMPLICIT DOUBLE PRECISION(A - H, O - Z)
            DIMENSION CON(5*ND), ICOMP(ND)
         END
      END INTERFACE

      INTEGER(C_INT) NR, N, ND, ICOMP(ND), IPARP, IRTRN, J, ITP
      REAL(C_DOUBLE) XOLD, X, CON(5*ND), RPARP, Y(NEQP), XOUTP
      LOGICAL(4) PRESSED
      CHARACTER(1) KEY
      INTEGER(C_INT), PARAMETER :: ESC = 27
      COMMON/INTERNP/XOUTP, ITP

      IF (NR .EQ. 1) THEN
         ITP = 1
         DO J = 1, NEQP
            P(J, ITP) = Y(J)
         END DO
         XOUTP = X + DZ
      ELSE
10       CONTINUE
         IF (X .GE. XOUTP) THEN
            ITP = ITP + 1
            DO J = 1, NEQP
               P(J, ITP) = CONTD5_P(J, XOUTP, CON, ICOMP, ND)
            END DO
            XOUTP = XOUTP + DZ
            GOTO 10
         END IF
      END IF
      RETURN
   END SUBROUTINE SOLOUTP

   SUBROUTINE SOLOUTF(NR, XOLD, X, Y, N, CON, ICOMP, ND, RPARF, IPARF, IRTRN)
      IMPLICIT NONE

      INTERFACE
         FUNCTION CONTD5_F(II, X, CON, ICOMP, ND)
            IMPLICIT DOUBLE PRECISION(A - H, O - Z)
            DIMENSION CON(5*ND), ICOMP(ND)
         END
      END INTERFACE

      INTEGER(C_INT) NR, N, ND, ICOMP(ND), IPARF, IRTRN, J, ITF
      REAL(C_DOUBLE) XOLD, X, CON(5*ND), RPARF, Y(NEQF), XOUTF, PEX(NEQP), YY(NEQF), WCUR(3)
      LOGICAL(4) PRESSED
      CHARACTER(1) KEY
      INTEGER(C_INT), PARAMETER :: ESC = 27
      COMMON/INTERNF/XOUTF, ITF
      !COMMON/INTERNF/ ITF

      IF (NR .EQ. 1) THEN
         ITF = 1
         DO J = 1, NEQF
            F(J, ITF) = Y(J)
         END DO
         CALL CALCPEX(Y, P, CL1(ITF), LHS1(ITF), RHS1(ITF), CL2(ITF), LHS2(ITF), RHS2(ITF))
         !ETA(:, ITF) = EFF(PEX)
         ETA(:, ITF) = EFF(P(:, NZ))
         !ETAG(:, ITF) = PITCH**2/(PITCH**2 + 1)*ETA(:, ITF)
         !WRITE (*, '(A,F10.5,A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F6.3,A,F6.3,\,A,A)') 'TIME = ', XOUTF, &
         WRITE (*, '(A,F8.3,A,F8.5,A,F8.5,A,F8.5,A,F8.5,A,F8.5,A,F9.5,A,F9.5,A,F9.5,A,F5.3,A,F5.3,A,\,A)') 'T =', XOUTF, &
            '  |F1| = ', ABS(F(1, ITF)), '  |F2| = ', ABS(F(3, ITF)), &
            '  |F3| = ', ABS(F(5, ITF)), '  EFF1 = ', ETA(1, ITF), '  EFF2 = ', ETA(2, ITF), &
            '  W1 = ', 0, '  W2 = ', 0, '  W3 = ', 0, &
            !'  PH1 = ', F(2, ITF), '  PH2 = ', F(4, ITF), '  PH3 = ', F(6, ITF), &
            '  C1 = ', ABS(CL1(ITF)/RHS1(ITF))*100, ' %  C2 = ', ABS(CL2(ITF)/RHS2(ITF))*100, ' %', CHAR(13)
         XOUTF = X + DT
         XI1(ITF) = XI(P, 1)
         XI2(ITF) = XI(P, 2)
         XI1_KR(ITF) = XI1(ITF)/(F(1, ITF)*CDEXP(IC*F(2, ITF)))
         XI2_KR(ITF) = XI2(ITF)/(F(3, ITF)*CDEXP(IC*F(4, ITF)))
      ELSE
10       CONTINUE
         IF (X .GE. XOUTF) THEN
            ITF = ITF + 1
            DO J = 1, NEQF
               F(J, ITF) = CONTD5_F(J, XOUTF, CON, ICOMP, ND)
            END DO
            CALL CALCPEX(Y, P, CL1(ITF), LHS1(ITF), RHS1(ITF), CL2(ITF), LHS2(ITF), RHS2(ITF))
            !ETA(:, ITF) = EFF(PEX)
            ETA(:, ITF) = EFF(P(:, NZ))
            !ETAG(:, ITF) = PITCH**2/(PITCH**2 + 1)*ETA(:, ITF)
            DO J = 1, 3
               WCUR(J) = (F(2*J, ITF) - F(2*J, ITF - 1))/DT
            END DO
            XI1(ITF) = XI(P, 1)
            XI2(ITF) = XI(P, 2)
            XI1_KR(ITF) = XI1(ITF)/(F(1, ITF)*CDEXP(IC*F(2, ITF)))
            XI2_KR(ITF) = XI2(ITF)/(F(3, ITF)*CDEXP(IC*F(4, ITF)))
            !WRITE (*, '(A,F10.5,A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F10.6,A,F6.3,A,F6.3,\,A,A)') 'TIME = ', XOUTF, &
            WRITE (*, '(A,F8.3,A,F8.5,A,F8.5,A,F8.5,A,F8.5,A,F8.5,A,F9.5,A,F9.5,A,F9.5,A,F5.3,A,F5.3,A,\,A)') 'T =', XOUTF, &
               '  |F1| = ', ABS(F(1, ITF)), '  |F2| = ', ABS(F(3, ITF)), &
               '  |F3| = ', ABS(F(5, ITF)), '  EFF1 = ', ETA(1, ITF), '  EFF2 = ', ETA(2, ITF), &
               '  W1 = ', WCUR(1), '  W2 = ', WCUR(2), '  W3 = ', WCUR(3), &
               '  C1 = ', DABS(CL1(ITF)/RHS1(ITF))*100, ' %  C2 = ', DABS(CL2(ITF)/RHS2(ITF))*100, ' %', CHAR(13)
            XOUTF = XOUTF + DT
            GOTO 10
         END IF
      END IF

      PRESSED = PEEKCHARQQ()
      IF (PRESSED) THEN
         KEY = GETCHARQQ()
         IF (ICHAR(KEY) .EQ. ESC) THEN
            WRITE (*, '(/,A)') 'QUIT?'
            KEY = GETCHARQQ()
            IF (ICHAR(KEY) .EQ. 121 .OR. ICHAR(KEY) .EQ. 89) THEN
               NT = ITF
               IRTRN = -1
               !RETURN
            END IF
         END IF
      END IF
      RETURN
   END SUBROUTINE SOLOUTF

   SUBROUTINE SOLOUT_RADAU(NR, XOLD, X, Y, CONT, LRC, N, RPAR, IPAR, IRTRN)
!C --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS BY USING "CONTR5"
      IMPLICIT REAL*8(A - H, O - Z)
      DIMENSION Y(N), CONT(LRC)
      COMMON/INTERN/XOUT
      IF (NR .EQ. 1) THEN
         WRITE (6, 99) X, Y(1), Y(3), Y(5), NR - 1
         XOUT = 0.02D0
      ELSE
10       CONTINUE
         IF (X .GE. XOUT) THEN
!C --- CONTINUOUS OUTPUT FOR RADAU5
            WRITE (6, 99) XOUT, CONTR5(1, XOUT, CONT, LRC), &
               CONTR5(3, XOUT, CONT, LRC), CONTR5(5, XOUT, CONT, LRC), NR - 1
            XOUT = XOUT + 0.02D0
            GOTO 10
         END IF
      END IF
99    FORMAT(1X, 'X =', F5.2, '    Y1 =', E18.10, '    Y3 =', E18.10, '    Y5 =', E18.10, '    NSTEP =', I10)
      RETURN
   END SUBROUTINE SOLOUT_RADAU

   SUBROUTINE SOLVEF_RK()
      IMPLICIT NONE

      !INTEGER(C_INT) NT, NEQ
      REAL(C_DOUBLE) H, T0, T

      CALL ODE4F_RK(DFDT_RK, F, NEQF, NT, 0.0D0, DT)
   END SUBROUTINE SOLVEF_RK

   FUNCTION DFDT_RK(T, Y) RESULT(S)
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      REAL(C_DOUBLE) T, Y(:), S(SIZE(Y))

      CALL DFDT_COMMON(NEQF, T, Y, S)
   END FUNCTION DFDT_RK

   SUBROUTINE ODE4F_RK(DYDT, Y, NEQ, NT, T0, H)
      IMPORT
      IMPLICIT NONE

      INTERFACE
         FUNCTION DYDT(T, Y) RESULT(S)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            REAL(C_DOUBLE) T, Y(:), S(SIZE(Y))
         END FUNCTION DYDT
      END INTERFACE

      INTEGER(C_INT) NT, NEQ, I, J, ITF
      REAL(C_DOUBLE) H, T0, T
      REAL(C_DOUBLE), INTENT(INOUT) :: Y(:, :)
      REAL(C_DOUBLE) S1(SIZE(Y, 1)), S2(SIZE(Y, 1)), S3(SIZE(Y, 1)), S4(SIZE(Y, 1)), V(SIZE(Y, 1))
      REAL(C_DOUBLE) PEX(NEQP)
      LOGICAL(4) PRESSED
      CHARACTER(1) KEY
      INTEGER(C_INT), PARAMETER :: ESC = 27

      !!SOLVE EQ. AT T=0
      !FP(1) = Y(1, 1)*CDEXP(IC*Y(2, 1))
      !FP(2) = Y(3, 1)*CDEXP(IC*Y(4, 1))
      !
      !CALL ODE4P(DPDZ_COMMON, P, NE, NZ, 0.0D0, DZ)
      !
      !ETA(:, 1) = EFF(P(:, NZ))
      !ETAG(:, 1) = PITCH**2/(PITCH**2 + 1)*ETA(:, 1)

      !OPEN (1, FILE='TEST.DAT')
      !DO J = 1, NZ
      !    WRITE (1, '(5F17.8)') (J-1)*DZ, P(14, J), P(NE+14, J), P(128, J), P(NE+128, J)
      !END DO
      !CLOSE (1)
      !STOP

      DO I = 1, NT - 1
         V = Y(:, I)
         T = T0 + (I - 1)*H
         S1(:) = DYDT(T, V)
         S2(:) = DYDT(T + H/2, V + H*S1(:)/2)
         S3(:) = DYDT(T + H/2, V + H*S2(:)/2)
         S4(:) = DYDT(T + H, V + H*S3(:))
         Y(:, I + 1) = V + H*(S1(:) + 2.0D0*S2(:) + 2.0D0*S3(:) + S4(:))/6.0D0

         ETA(:, I + 1) = EFF(P(:, NZ))
         !ETAG(:, I + 1) = PITCH**2/(PITCH**2 + 1)*ETA(:, I + 1)

         ITF = I + 1
         CALL CALCPEX(Y(:, I + 1), P, CL1(ITF), LHS1(ITF), RHS1(ITF), CL2(ITF), LHS2(ITF), RHS2(ITF))
         ETA(:, ITF) = EFF(PEX)
         !ETA(:, ITF) = EFF(P(:, NZ))
         !ETAG(:, ITF) = PITCH**2/(PITCH**2 + 1)*ETA(:, ITF)
         W(1, ITF) = (S1(2) + 2.0D0*S2(2) + 2.0D0*S3(2) + S4(2))/6.0D0
         W(2, ITF) = (S1(4) + 2.0D0*S2(4) + 2.0D0*S3(4) + S4(4))/6.0D0
         W(3, ITF) = (S1(6) + 2.0D0*S2(6) + 2.0D0*S3(6) + S4(6))/6.0D0

         !WRITE (*, '(A,F12.7,A,F10.7,A,F10.7,A,F10.7,A,F10.7,A,F10.7,A,F5.3,A,F5.3,A,\,A)') 'TIME = ', T + H, '   |F1| = ', ABS(Y(1, I + 1)), '   |F2| = ', ABS(Y(3, I + 1)), '   |F3| = ', ABS(Y(5, I + 1)), &
         !'   EFF1 = ', ETA(1, I + 1), '   EFF2 = ', ETA(2, I + 1), '  C1 = ', DABS(CL1(ITF)/RHS1(ITF))*100, ' %  C2 = ', DABS(CL2(ITF)/RHS2(ITF))*100, ' %', CHAR(13)
         WRITE (*, '(A,F8.3,A,F8.5,A,F8.5,A,F8.5,A,F8.5,A,F8.5,A,F9.5,A,F9.5,A,F9.5,A,F5.3,A,F5.3,A,\,A)') 'T =', T + H, &
            '  |F1| =', ABS(Y(1, I + 1)), '  |F2| =', ABS(Y(3, I + 1)), '  |F3| =', ABS(Y(5, I + 1)), '  E1 =', ETA(1, I + 1), '  E2 =', ETA(2, I + 1), &
            '  W1 = ', W(1, ITF), '  W2 = ', W(2, ITF), '  W3 = ', W(3, ITF), '  C1 = ', DABS(CL1(ITF)/RHS1(ITF))*100, '%  C2 = ', DABS(CL2(ITF)/RHS2(ITF)*100), '%', CHAR(13)

         PRESSED = PEEKCHARQQ()
         IF (PRESSED) THEN
            KEY = GETCHARQQ()
            IF (ICHAR(KEY) .EQ. ESC) THEN
               WRITE (*, '(/,A)') 'QUIT?'
               KEY = GETCHARQQ()
               IF (ICHAR(KEY) .EQ. 121 .OR. ICHAR(KEY) .EQ. 89) THEN
                  NT = I + 1
                  RETURN
                  !OPEN (1, FILE='F.DAT')
                  !DO J = 1, I + 1
                  !    WRITE (1, '(4E17.8)') (J - 1)*H, ABS(Y(1, J)), ABS(Y(2, J)), ABS(Y(3, J))
                  !END DO
                  !CLOSE (2)
                  !OPEN (2, FILE='E.DAT')
                  !DO J = 1, I + 1
                  !    WRITE (2, '(5E17.8)') (J - 1)*H, ETA(1, J), ETAG(1, J), ETA(2, J), ETAG(2, J)
                  !END DO
                  !CLOSE (2)
                  !STOP
               END IF
            END IF
         END IF
      END DO
   END SUBROUTINE ODE4F_RK

   SUBROUTINE SOLVEP_RK(PIN, PEX, C)
      IMPLICIT NONE

      REAL(C_DOUBLE), INTENT(IN) :: PIN(:)
      REAL(C_DOUBLE), INTENT(INOUT) :: PEX(:, :)
      CHARACTER(C_CHAR), INTENT(IN) :: C

      PEX(:, 1) = PIN

      IF (C .EQ. 'P') THEN
         CALL ODE4P_RK(DPDZ_RK, PEX, NEQP, NZ, 0.0D0, DZ)
      ELSE
         CALL ODE4P_RK(DPDZ_RK, P, NEQP, NZ, 0.0D0, DZ)
         PEX(:, 1) = P(:, NZ)
      END IF
   END SUBROUTINE SOLVEP_RK

   SUBROUTINE ODE4P_RK(DYDT, Y, NEQ, NZ, Z0, H)!, PARAMS)
      IMPORT
      IMPLICIT NONE

      INTERFACE
         SUBROUTINE DYDT(Z, P, PRHS)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            REAL(C_DOUBLE), INTENT(IN) :: Z, P(*)
            REAL(C_DOUBLE), INTENT(OUT) :: PRHS(*)
         END SUBROUTINE DYDT
      END INTERFACE

      INTEGER(C_INT) NZ, NEQ, I
      REAL(C_DOUBLE) H, Z0, Z
      REAL(C_DOUBLE), INTENT(INOUT) :: Y(:, :)
      REAL(C_DOUBLE) S1(SIZE(Y, 1)), S2(SIZE(Y, 1)), S3(SIZE(Y, 1)), S4(SIZE(Y, 1)), V(SIZE(Y, 1))

      DO I = 1, NZ - 1
         V = Y(:, I)
         Z = Z0 + (I - 1)*H
         CALL DYDT(Z, V, S1)
         CALL DYDT(Z + H/2, V + H*S1(:)/2, S2)
         CALL DYDT(Z + H/2, V + H*S2(:)/2, S3)
         CALL DYDT(Z + H, V + H*S3(:), S4)
         Y(:, I + 1) = V + H*(S1(:) + 2*S2(:) + 2*S3(:) + S4(:))/6
      END DO
   END SUBROUTINE ODE4P_RK

   SUBROUTINE DPDZ_RK(Z, P, PRHS)
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      REAL(C_DOUBLE), INTENT(IN) :: Z, P(*)
      REAL(C_DOUBLE), INTENT(OUT) :: PRHS(*)

      CALL DPDZ_COMMON(Z, P, PRHS)

   END SUBROUTINE DPDZ_RK

END MODULE FUN
