      SUBROUTINE URAND1(N, X, IR)
************************************************************************
* UNIFORM RANDOM NUMBER GENERATOR (MIXED CONGRUENTIAL METHOD)          *
*     PORTABLE BUT SLOW.  THE PERIOD IS ONLY 1664501.                  *
* PARAMETERS                                                           *
*   (1) N      (I) THE NUMBER OF RANDOM NUMBERS TO BE GENERATED        *
*                  (INPUT)                                             *
*   (2) X      (D) UNIFORM RANDOM NUMBERS (OUTPUT)                     *
*   (3) IR     (I) THE INITIAL SEED  (INPUT)                           *
*                  THE SEED FOR THE NEXT CALL (OUTPUT)                 *
* COPYRIGHT: Y. OYANAGI, JUNE 30, 1989  V.1                            *
************************************************************************
*
       DOUBLE PRECISION X(N), INVM
       PARAMETER (M = 1664501, LAMBDA = 1229, MU = 351750)
       PARAMETER (INVM = 1.0D0 / M)
*PAREMETER CHECK
      IF( N .LE. 0) THEN
       WRITE(6,*) '(SUBR.URAND1) PARAMETER ERROR. N = ', N
       WRITE(6,*) 'RETURN WITH NO FURTHER CALCULATION.'
       RETURN
      END IF
      IF( IR .LT. 0 .OR. IR .GE. M) THEN
       WRITE(6,*) '(SUBR.URAND1) WARNING. IR = ', IR
      END IF
*MAIN LOOP
      DO 10 I = 1, N
       IR = MOD( LAMBDA * IR + MU, M)
       X(I) = IR * INVM
   10 CONTINUE
      RETURN
      END

