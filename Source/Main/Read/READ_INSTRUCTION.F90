!***************************************************************************************************
!- PURPORSE:
!     READ INSTRUCTIONS IN ONE STEP BLOCK.
!  
!- INPUT ARGUMENTS:
!  NONE
!
!- OUTPUT ARGUMENTS:
!  NONE
!
!- CALL PROCEDURES:
!  MODULE PROCEDURE: GET_MACRO, REMOVE_FIRSTWORD, ERROR
!
!- CALLED BY:
!  SOLVER_MANAGER
!
!- PROGAMMED BY:
!  ZHIHAI XIANG, DEPARTMENT OF ENGINEERING MECHANICS, TSINGHUA UNIVERSITY, NOVEMBER 1, 2015
!***************************************************************************************************

SUBROUTINE READ_INSTRUCTION

USE SOLUTION_DATA, ONLY: INT_KIND, LINE_KIND, WORD_KIND, CONTROL, NUM_OPT, INSTRUCTION_NODE, &
                         INSTRUCTIONS, GET_MACRO, REMOVE_FIRSTWORD, ERROR

IMPLICIT NONE

CHARACTER(LINE_KIND)              COMMAND
CHARACTER(WORD_KIND)              INSTRUCTION
CHARACTER(WORD_KIND)              OPTION(NUM_OPT) ! INSTRUCTION OPTIONS
INTEGER(INT_KIND)                 ERR             ! ERR INDEX
INTEGER(INT_KIND)                 I               ! LOOP INDEX
TYPE(INSTRUCTION_NODE), POINTER:: P_INSTRUCTION

!-- READ EACH INSTRUCTION --
COMMAND = GET_MACRO(CONTROL)
IF(COMMAND(1:8) .NE. 'SOLUTION') CALL ERROR("'SOLUTION' is expected!")

DO
   COMMAND = GET_MACRO(CONTROL) 
   READ(COMMAND,*) INSTRUCTION
   COMMAND = REMOVE_FIRSTWORD(COMMAND)       
   
   ! COULD BE 'END SOLUTION'
   IF(TRIM(INSTRUCTION) .EQ. 'END') THEN
      READ(COMMAND,*) INSTRUCTION
      COMMAND = REMOVE_FIRSTWORD(COMMAND)

      IF(TRIM(INSTRUCTION) .EQ. 'SOLUTION') THEN
         EXIT
      ELSE   
         CALL ERROR("Invalid command '" // TRIM(INSTRUCTION) // "' in SOLUTION block!")      
      END IF
   ENDIF  
   
   ! READ OPTIONS OF EACH INSTRUCTION
   OPTION = ''
   SELECT CASE(TRIM(INSTRUCTION))
      ! STIFFNESS MATRIX
      CASE ('STIFFNESS')
         IF(LEN_TRIM(COMMAND) .GT. 0) THEN   
            READ(COMMAND,*) OPTION(1)
            IF(TRIM(OPTION(1)) .NE. 'PLANE_STRESS' .AND. &
               TRIM(OPTION(1)) .NE. 'PLANE_STRAIN') &
               CALL ERROR("Invalid STIFFNESS option '" // TRIM(OPTION(1)) // "'!")
         ENDIF
            
      CASE ('BOUNDARY')
      
      CASE ('FORCE')
      
      CASE ('SOLVE')
         I = 1
         DO WHILE(LEN_TRIM(COMMAND) .GT. 0)
            READ(COMMAND,*) OPTION(I)
            COMMAND = REMOVE_FIRSTWORD(COMMAND)
            I       = I + 1
         ENDDO
         
      CASE ('STRAIN')
         IF(LEN_TRIM(COMMAND) .GT. 0) THEN   
            READ(COMMAND,*) OPTION(1)
            IF(TRIM(OPTION(1)) .NE. 'NODAL' .AND. &
               TRIM(OPTION(1)) .NE. 'GAUSS') &
               CALL ERROR("Invalid STRAIN option '" // TRIM(OPTION(1)) // "'!")            
         ELSE
            CALL ERROR("Nodal or Gauss STRAIN should be specified!")
         ENDIF
                  
      CASE ('STRESS')
         IF(LEN_TRIM(COMMAND) .GT. 0) THEN   
            READ(COMMAND,*) OPTION(1)
            IF(TRIM(OPTION(1)) .NE. 'NODAL' .AND. &
               TRIM(OPTION(1)) .NE. 'GAUSS') &
               CALL ERROR("Invalid STRESS option '" // TRIM(OPTION(1)) // "'!")
            
            COMMAND = REMOVE_FIRSTWORD(COMMAND)
            IF(LEN_TRIM(COMMAND) .GT. 0) THEN
               READ(COMMAND,*) OPTION(2)
               IF(TRIM(OPTION(2)) .NE. 'PLANE_STRESS' .AND. &
               TRIM(OPTION(2)) .NE. 'PLANE_STRAIN') &
               CALL ERROR("Invalid STRESS option(2) '" // TRIM(OPTION(2)) // "'!")
            ENDIF
         ELSE
            CALL ERROR("Nodal or Gauss STRESS should be specified!")
         ENDIF
         
      CASE ('OUTPUT')
         I = 1
         DO WHILE(LEN_TRIM(COMMAND) .GT. 0)
            READ(COMMAND,*) OPTION(I)
            COMMAND = REMOVE_FIRSTWORD(COMMAND)
            I       = I + 1
         ENDDO 
         
      CASE DEFAULT
         CALL ERROR('Invalid instruction: ' // TRIM(INSTRUCTION) // ' in SOLUTION block!')
   END SELECT
   
   ! CREATE INSTRUCTION LIST
   ALLOCATE(P_INSTRUCTION, STAT = ERR)
   IF(ERR .NE. 0) CALL ERROR('Fail to allocate an instruction node!')
   
   P_INSTRUCTION%INSTRUCTION = INSTRUCTION
   P_INSTRUCTION%OPTION      = OPTION
         
   IF(ASSOCIATED(INSTRUCTIONS%HEAD)) THEN
      INSTRUCTIONS%TAIL%NEXT => P_INSTRUCTION
      P_INSTRUCTION%PRE      => INSTRUCTIONS%TAIL
      INSTRUCTIONS%TAIL      => P_INSTRUCTION
   ELSE
      INSTRUCTIONS%HEAD => P_INSTRUCTION
      INSTRUCTIONS%TAIL => P_INSTRUCTION
   ENDIF
ENDDO

END SUBROUTINE READ_INSTRUCTION