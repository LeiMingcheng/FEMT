!***************************************************************************************************
!- PURPORSE:
!     BRANCH TO DIFFERENT KIND OF ANALYSIS.
!  
!- INPUT ARGUMENTS:
!  OPTION: SOLUTION OPTIONS
!
!- OUTPUT ARGUMENTS:
!  NONE
!
!- CALL PROCEDURES:
!  MODULE PROCEDURE:   ERROR
!  EXTERNAL PROCEDURE: SOLVE_MECH_STATIC_DIR
!
!- CALLED BY:
!  SOLVER_MANAGER  
!
!- PROGAMMED BY:
!  ZHIHAI XIANG, DEPARTMENT OF ENGINEERING MECHANICS, TSINGHUA UNIVERSITY, JANUARY 17, 2016
!***************************************************************************************************

SUBROUTINE SOLVE(OPTION)

USE SOLUTION_DATA, ONLY: INT_KIND, WORD_KIND, NUM_OPT, ERROR

IMPLICIT NONE

CHARACTER(WORD_KIND), INTENT(IN):: OPTION(NUM_OPT) ! SOLUTION OPTIONS

!-- BRANCH TO DIFFERENT KIND OF ANALYSIS --
SELECT CASE(TRIM(OPTION(1)))
   ! PHYSICAL FIELD
   CASE ('MECHANICAL')
      ! ANALYSIS TYPE
      SELECT CASE (TRIM(OPTION(2)))
         CASE ('STATIC')
            ! SOLVER TYPE
            SELECT CASE(TRIM(OPTION(3)))
               CASE ('DIR')   ! DIRECT ELIMINATION
                  CALL SOLVE_MECH_STATIC_DIR
                  
               CASE ('CGRAD') ! CONJUGATION GRADIENT ITERATION
                  
               CASE DEFAULT
                  CALL ERROR('Invalid equation solver ' // TRIM(OPTION(3)) // '!')
            END SELECT
      
         CASE ('DYNAMIC')
         
         CASE ('MODAL')  

         CASE ('BUCKLING')
        
         CASE DEFAULT
            CALL ERROR('Invalid mechanical problem type ' // TRIM(OPTION(2)) // '!')
      END SELECT

   CASE ('THERMAL')
      
   CASE ('COUPLE')
   
   CASE DEFAULT
      CALL ERROR('Invalid problem class ' // TRIM(OPTION(1)) // '!')
END SELECT

END SUBROUTINE SOLVE