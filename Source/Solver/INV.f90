!***************************************************************************************************
!- PURPORSE:
!     CALCULATE THE INVERSE OF MATRIX A.
!  
!- INPUT ARGUMENTS:
!  A: THE N*N MATRIX 
!  N: THE DEMENSION OF MATRIX A
!
!- OUTPUT ARGUMENTS:
!  A: THE INVERSE OF MATRIX A
!
!- CALL PROCEDURES:
!  MODULE PROCEDURE: NUM2STR, ERROR
!
!- CALLED BY
!  STRAIN_NODAL_RECTANGLE, STRESS_NODAL_RECTANGLE
!
!- PROGAMMED BY:
!  ZHIHAI XIANG, DEPARTMENT OF ENGINEERING MECHANICS, TSINGHUA UNIVERSITY, JANUARY 18, 2016
!***************************************************************************************************

SUBROUTINE INV(A,N)

USE BASIC_DATA, ONLY: INT_KIND, REAL_KIND, REAL_EPS, NUM2STR, ERROR

IMPLICIT NONE

INTEGER(INT_KIND), INTENT(IN)::  N           ! THE DEMENSION OF MATRIX A
REAL(REAL_KIND), INTENT(INOUT):: A(N,N)      ! THE N*N MATRIX
INTEGER(INT_KIND)                I, J, K     ! LOOP INDEX
INTEGER(INT_KIND)                ME(N)       ! RECORD THE COLUMN OF THE LARGEST ABSOLUTE VALUE
REAL(REAL_KIND)                  Y           ! THE LARGEST ABSOLUTE VALUE OR A TEMPERARY VARIABLE        
REAL(REAL_KIND)                  B(N), C(N)  ! TEMPERARY ARRAY        

!-- DECOMPOSE THE MATRIX --
DO I = 1, N
   ME(I) = I
ENDDO

DO I = 1, N
   ! SELECT THE LARGEST ABSOLUTE VALUE IN THIS ROW
   Y = 0.0
   DO J = I, N
      IF(ABS(A(I,J)) .LE. ABS(Y)) CYCLE
      K = J       ! THE COLUMN OF THE LARGEST ABSOLUTE VALUE
      Y = A(I,J)
   ENDDO
        
   IF(ABS(Y) .LT. REAL_EPS) CALL ERROR('Singular matrix')
        
   Y = 1.0 / Y
   DO J = 1, N
      C(J)   = A(J,K)
      A(J,K) = A(J,I)
      A(J,I) = -C(J)  * Y
      B(J)   = A(I,J) * Y
      A(I,J) = A(I,J) * Y
   ENDDO
   
   A(I,I) = Y
   ! SWARP COLUMN
   J     = ME(I)
   ME(I) = ME(K)
   ME(K) = J
      
   DO K = 1, N
      IF(K .EQ. I) CYCLE
      DO J = 1, N
         IF(J .EQ. I) CYCLE
         A(K,J) = A(K,J) - B(J) * C(K)
      ENDDO
   ENDDO
ENDDO

!-- BACK SUBSTITUED --
DO I = 1, N
   DO K = 1, N
      IF(ME(K) .EQ. I) EXIT
   ENDDO
   
   IF(K .EQ. I) CYCLE
      
   DO J = 1, N
      Y      = A(I,J)
      A(I,J) = A(K,J)
      A(K,J) = Y
   ENDDO
        
   J     = ME(I)
   ME(I) = ME(K)
   ME(K) = J
ENDDO
     
END SUBROUTINE INV