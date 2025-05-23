!***************************************************************************************************
!- PURPORSE:
!     WRITE RESULTS TO OUTPUT FILE.
!  
!- INPUT ARGUMENTS:
!  OPTION: THE TYPE OF OUTPUT DATA
!
!- OUTPUT ARGUMENTS:
!  NONE
!
!- CALL PROCEDURES:
!  MODULE PROCEDURE : ERROR
!
!- CALLED BY:
!  SOLVER_MANAGER  
!
!- PROGAMMED BY:
!  ZHIHAI XIANG, DEPARTMENT OF ENGINEERING MECHANICS, TSINGHUA UNIVERSITY, JANUARY 18, 2016
!***************************************************************************************************

SUBROUTINE OUTPUT(OPTION)

USE SOLUTION_DATA, ONLY: INT_KIND, REAL_KIND, WORD_KIND, OUT, NUM_OPT, NUM_NODE,  &
                         NODE_DISP_DOF, DISPLACEMENT, NODAL_STRAIN, NODAL_STRESS, ERROR

IMPLICIT NONE

CHARACTER(WORD_KIND), INTENT(IN):: OPTION(NUM_OPT)  ! THE TYPE OF OUTPUT DATA
INTEGER(INT_KIND)                  I, J, K          ! LOOP INDEX
REAL(REAL_KIND)                    B, C, S1, S2     ! PRINCIPAL STRAINS OR STRESSES

!-- OUTPUT RESULTS --
DO I = 1, NUM_OPT
   IF(LEN_TRIM(OPTION(I)) .EQ. 0) EXIT
   
   SELECT CASE(TRIM(OPTION(I)))                
      CASE ('DISPLACEMENT')
         WRITE(OUT,*) '*** DISPLACEMENT ***'
         WRITE(OUT,'(A)') '  NODE         U              V'
         WRITE(OUT,"(35('-'))")
         DO J = 1, NUM_NODE
            WRITE(OUT,'(I7, 7(2X, ES13.6))')&
               J, (DISPLACEMENT(NODE_DISP_DOF(J)%DOFS(K)), K = 1, NODE_DISP_DOF(J)%NUM_DOF)
         ENDDO
         
      CASE ('NODAL_STRAIN')
         WRITE(OUT,*) '*** NODAL STRAIN ***'
         WRITE(OUT,'(A)') '  NODE         EXX            EYY            EXY            E1            E2'
         WRITE(OUT,"(82('-'))")
         DO J = 1, NUM_NODE
            B = NODAL_STRAIN(J)%VALUE(1) + NODAL_STRAIN(J)%VALUE(2)
            C = NODAL_STRAIN(J)%VALUE(1) * NODAL_STRAIN(J)%VALUE(2) - NODAL_STRAIN(J)%VALUE(3)**2
            S1 = (B + SQRT(B**2 - 4.0 * C)) / 2.0
            S2 = (B - SQRT(B**2 - 4.0 * C)) / 2.0
            WRITE(OUT,'(I7, 6(2X, ES13.6))')  J, NODAL_STRAIN(J)%VALUE, S1, S2
         ENDDO
         
      CASE ('NODAL_STRESS')
         WRITE(OUT,*) '*** NODAL STRESS ***'
         WRITE(OUT,'(A)') '  NODE         SXX            SYY            SXY            S1            S2'
         WRITE(OUT,"(82('-'))")
         DO J = 1, NUM_NODE
            B = NODAL_STRESS(J)%VALUE(1) + NODAL_STRESS(J)%VALUE(2)
            C = NODAL_STRESS(J)%VALUE(1) * NODAL_STRESS(J)%VALUE(2) - NODAL_STRESS(J)%VALUE(3)**2
            S1 = (B + SQRT(B**2 - 4.0 * C)) / 2.0
            S2 = (B - SQRT(B**2 - 4.0 * C)) / 2.0
            WRITE(OUT,'(I7, 6(2X, ES13.6))')  J, NODAL_STRESS(J)%VALUE, S1, S2
         ENDDO
         
      CASE DEFAULT
         CALL ERROR("Invalid output item: '" // TRIM(OPTION(I)) // "'!")
   END SELECT
   
   WRITE(OUT,*)
ENDDO

END SUBROUTINE OUTPUT