MODULE CSPLINE

CONTAINS

   SUBROUTINE SPLINE(N, X, Y, B, C, D)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(N), Y(N), B(N), C(N), D(N)
!
!  THE COEFFICIENTS B(I), C(I) AND D(I) ARE CALCULATED, 1=1,
!  2, ..., N, FOR CUBIC INTERPOLATION SPLINE
!
!  S(X) = Y(I)+B(I)*(X-X(I)) + C(I)*(X-X(I))**2 +
!  -FD(I)*(X - X(I))**3
!
!  FOR X(I) .LE. X .LE. X(I+1)
!
!  INPUT INFORMATION..
!
!  N = NUMBER OF SPECIFIED POINTS OR NODES (N .GE. 2)
!  X = ABSCISSUE OF NODES IN STRICTLY INCREASING ORDER
!  Y = ORDINATES OF NODES
!
!  OUTPUT...
!
!  B, C, D = ARRAYS OF SPLINE COEFFICIENTS DEFINITED ABOVE.
!
!  IF YOU DESIGNATE THE DIFFERENTIATION SYMBOL BY P, THEN
!
!  Y(I)= S(X(I))
!  B(I) = SP(X(I))
!  C(I) = SPP(X(I))/2
!  D(I) = SPPP(X(I))/6 (RIGHT HAND DERIVATIVE)
!
!  USING THE ACCOMPANYING SEVAL FUNCTION SUBROUTINE
!  YOU CAN CALCULATE SPLINE VALUES.
!
      INTEGER NM1, IB, I
      REAL T

      NM1 = N - 1
      IF (N .LT. 2) RETURN
      IF (N .LT. 3) GO TO 50
!
! BUILD A TRIDIAGONAL SYSTEM
! B = DIAGONAL, O = OVERDIAGONAL, C = RIGHT PARTS.
!
      D(1) = X(2) - X(1)
      C(2) = (Y(2) - Y(1))/D(1)
      DO I = 2, NM1
         D(I) = X(I + 1) - X(I)
         B(I) = 2.*(D(I - 1) + D(I))
         C(I + 1) = (Y(I + 1) - Y(I))/D(I)
         C(I) = C(I + 1) - C(I)
      END DO
!
! BOUNDARY CONDITIONS. THIRD DERIVATIVES AT POINTS
! X(1) AND X(N) ARE CALCULATED USING DIVISIONED
! DIFFERENCES
!
      B(1) = -D(1)
      B(N) = -D(N - 1)
      C(1) = 0.
      C(N) = 0.
      IF (N .EQ. 3) GO TO 15
      C(1) = C(3)/(X(4) - X(2)) - C(2)/(X(3) - X(1))
      C(N) = C(N - 1)/(X(N) - X(N - 2)) - C(N - 2)/(X(N - 1) - X(N - 3))
      C(1) = C(1)*D(1)**2/(X(4) - X(1))
      C(N) = -C(N)*D(N - 1)**2/(X(N) - X(N - 3))
!
! STRAIGHT RUN
!
15    DO I = 2, N
         T = D(I - 1)/B(I - 1)
         B(I) = B(I) - T*D(I - 1)
         C(I) = C(I) - T*C(I - 1)
      END DO
!
! REVERSE SUBSTITUSTION
!
      C(N) = C(N)/B(N)
      DO IB = 1, NM1
         I = N - IB
         C(I) = (C(I) - D(I)*C(I + 1))/B(I)
      END DO
!
! C(I) NOW STORES THE VALUE OF SIGMA(I), DEFINED
! IN #4.4.
!
! CALCULATE COEFFICIENTS OF POLYNOMIALS
!
      B(N) = (Y(N) - Y(NM1))/D(NM1) + D(NM1)*(C(NM1) + 2.*C(N))
      DO I = 1, NM1
         B(I) = (Y(I + 1) - Y(I))/D(I) - D(I)*(C(I + 1) + 2.*C(I))
         D(I) = (C(I + 1) - C(I))/D(I)
         C(I) = 3.*C(I)
      END DO
      C(N) = 3.*C(N)
      D(N) = D(N - 1)
      RETURN

50    B(1) = (Y(2) - Y(1))/(X(2) - X(1))
      C(1) = 0.
      D(1) = 0.
      B(2) = B(1)
      C(2) = 0.
      D(2) = 0.
      RETURN
   END SUBROUTINE SPLINE

   DOUBLE PRECISION FUNCTION SEVAL(N, U, X, Y, B, C, D)
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION U, X(N), Y(N), B(N), C(N), D(N)
!
!THIS SUBROUTINE EVALUATES THE CUBIC SPLINE FUNCTION
!
!SEVAL = Y(I)+B(I)*(U-X(I)) + C(I)*(U-X(I)))**2 + D(I)*(U-X(I))**3
!
!WHERE X(I) .LT. U .LT. X(I + 1). USING HORNER'S RULE
!
!IF U .LT. X(1), THEN THE VALUE 1 = 1 IS TAKEN.
!IF U .GE. X(N), THEN THE VALUE I = N IS TAKEN.
!
!INPUT..
!
!N = THE NUMBER OF DATA POINTS
!U = THE ABSCISSA AT WHICH THE SPLINE IS TO BE EVALUATED
!X, Y = THE ARRAYS OF DATA ABSCISSAS AND ORD1NATES
!B, C, D = ARRAYS OF SPLINE COEFFICIENTS, COMPUTED BY SPLINE SUBROUTINE
!
!IF U IS NOT IN THE SAME INTERVAL AS THE PREVIOUS CALL, THEN A
!BINARY SEARCH IS PERFORMED TO DETERMINE THE PROPER INTERVAL.
!
      INTEGER I, J, K
      REAL DX
      DATA I/1/
      IF (I .GE. N) I = 1
      IF (U .LT. X(I)) GO TO 10
      IF (U .LE. X(I + 1)) GO TO 30

!
! BINARY SEARCH
!
10    I = 1
      J = N + 1
20    K = (I + J)/2
      IF (U .LT. X(K)) J = K
      IF (U .GE. X(K)) I = K
      IF (J .GT. I + 1) GO TO 20
!
! EVALUATE SPLINE
!
30    DX = U - X(I)
      SEVAL = Y(I) + DX*(B(I) + DX*(C(I) + DX*D(I)))
      RETURN
   END FUNCTION SEVAL

   SUBROUTINE CMPLXSPLINE(N, ZEX, Z_ANG, RE_ANG, IM_ANG, REB, REC, RED, IMB, IMC, IMD)
      IMPLICIT NONE

      INTEGER(4) N
      REAL(8) C, ZEX, Z_ANG(N), RE_ANG(N), IM_ANG(N), REB(N), REC(N), RED(N), IMB(N), IMC(N), IMD(N)

      C = ZEX/Z_ANG(N)

      Z_ANG = C*Z_ANG

      CALL SPLINE(N, Z_ANG, RE_ANG, REB, REC, RED)
      CALL SPLINE(N, Z_ANG, IM_ANG, IMB, IMC, IMD)

      !
      !
      !DZ=Z_ANG(N)/2999
      !open(1,file='test.dat')
      !do i=1,3000
      !    TEST = CMPLXEVAL(N, (I-1)*DZ, Z_ANG, RE_ANG, IM_ANG, REB, REC, RED, IMB, IMC, IMD)
      !    write(1,'(3f12.6)') (I-1)*DZ, REAL(TEST), IMAG(TEST)
      !enddo
      !close(1)
      !STOP

   END SUBROUTINE CMPLXSPLINE

   COMPLEX(8) FUNCTION CMPLXEVAL(N, U, Z_ANG, RE_ANG, IM_ANG, REB, REC, RED, IMB, IMC, IMD)
      IMPLICIT NONE

      INTEGER(4) N
      REAL(8) U, Z_ANG(N), RE_ANG(N), IM_ANG(N), REB(N), REC(N), RED(N), IMB(N), IMC(N), IMD(N), RE, IM

      RE = SEVAL(N, U, Z_ANG, RE_ANG, REB, REC, RED)
      IM = SEVAL(N, U, Z_ANG, IM_ANG, IMB, IMC, IMD)

      CMPLXEVAL = DCMPLX(RE, IM)

   END FUNCTION CMPLXEVAL

END MODULE CSPLINE
