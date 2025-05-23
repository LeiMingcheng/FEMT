!***************************************************************************************************
!- PURPORSE:
!     READ ELEMENT INFORMATION FROM INPUT FILE.
!
!- INPUT ARGUMENTS:
!  NONE
!
!- OUTPUT ARGUMENTS:
!  NONE
!
!- CALL PROCEDURES:
!  MODULE PROCEDURE: GET_MACRO, ERROR, REMOVE_FIRSTWORD, NUM2STR
!
!- CALLED BY:
!  INITIATE
!
!- PROGAMMED BY:
!  ZHIHAI XIANG, DEPARTMENT OF ENGINEERING MECHANICS, TSINGHUA UNIVERSITY, NOVEMBER 1, 2015
!***************************************************************************************************

SUBROUTINE READ_ELEMENTS

USE BASIC_DATA, ONLY: INT_KIND, WORD_KIND, LINE_KIND, LG_KIND, CONTROL, NUM_ELEMENT, ELEMENT_LIB,  &
                      ELEMENTS, REMOVE_FIRSTWORD, NUM2STR, GET_MACRO, ERROR

IMPLICIT NONE

CHARACTER(LINE_KIND) COMMAND
CHARACTER(LINE_KIND) STRING     ! 'TO' OR THE LIST OF ELEMENT NUMBER
CHARACTER(WORD_KIND) ELE_TYPE   ! ELEMENT TYPE
INTEGER(INT_KIND)    ERR        ! ERR INDEX
INTEGER(INT_KIND)    I          ! LOOP INDEX
INTEGER(INT_KIND)    J          ! LOOP INDEX
INTEGER(INT_KIND)    FIRST
INTEGER(INT_KIND)    LAST
INTEGER(INT_KIND)    NUMBER     ! MATERIAL OR SECTION NUMBER OR THE NUMBER OF NODES
INTEGER(INT_KIND)    NO         ! NODE OR ELEMENT NUMBER
LOGICAL(LG_KIND)     EX         ! EXIT MARK

!-- ALLOCATE ELEMENTS --
IF(NUM_ELEMENT .LE. 0) CALL ERROR('NUM_ELEMENT should be greater than 0!')
ALLOCATE(ELEMENTS(NUM_ELEMENT), STAT = ERR)
IF(ERR .NE. 0) CALL ERROR('Fail to allocate ELEMENTS!')

!-- READ ELEMENT INFORMATION --
EX = .FALSE.
DO WHILE(.NOT. EX)
   COMMAND = GET_MACRO(CONTROL)

   SELECT CASE (TRIM(COMMAND))
      CASE ('END ELEMENTS')
         EX = .TRUE.

      CASE ('ELEMENT_TYPE')
         ! READ ELEMENT TYPE
         COMMAND = GET_MACRO(CONTROL)
         DO WHILE(TRIM(COMMAND) .NE. 'END ELEMENT_TYPE')
            ! LOAD THE 4TH CHARACTER AFTER THE 'TYPE' AS THE ELEMENT TYPE
            READ(COMMAND(INDEX(COMMAND, 'TYPE', BACK = .TRUE.) + 4:),*) ELE_TYPE
            SELECT CASE (TRIM(ELE_TYPE))
               CASE ('TRIANGLE3')
                  NUMBER = 1

               CASE ('TRIANGLE6')
                  NUMBER = 2

               CASE ('RECTANGLE4')
                  NUMBER = 3

               CASE ('RECTANGLE8')
                  NUMBER = 4

               CASE ('RECTANGLE9')
                  NUMBER = 5

               CASE DEFAULT
                  CALL ERROR('Invalid element type: ' // TRIM(ELE_TYPE) // ' in READ_ELEMENT!')
            END SELECT

            IF(INDEX(COMMAND, 'TO') .NE. 0) THEN
               ! A GROUP OF ELEMENTS ARE OF THAT TYPE
               READ(COMMAND,*) FIRST, STRING, LAST
               IF(LAST .LT. FIRST) CALL ERROR('Last element No. < first element No. ! (TYPE)')

               DO I = FIRST, LAST
                  ELEMENTS(I)%ELEMENT = NUMBER
               ENDDO
            ELSE
               ! A LIST OF ELEMENTS OF THAT TYPE
               STRING = COMMAND(1:INDEX(COMMAND, 'TYPE') - 1)
               DO WHILE(LEN_TRIM(STRING) .GT. 0)
                  READ(STRING,*) I
                  ELEMENTS(I)%ELEMENT = NUMBER
                  STRING = REMOVE_FIRSTWORD(STRING)
               ENDDO
            ENDIF

            COMMAND = GET_MACRO(CONTROL)
         ENDDO

         ! ALLOCATE MEMORIES FOR ELEMENT NODES
         DO I = 1, NUM_ELEMENT
            ALLOCATE(ELEMENTS(I)%NODE(ELEMENT_LIB(ELEMENTS(I)%ELEMENT)%NUM_NODE), STAT = ERR)
            IF(ERR .NE. 0) CALL ERROR('Fail to allocate node for element ' // TRIM(NUM2STR(I,'I')) //'!')
         ENDDO

      CASE ('ELEMENT_NODES')
         ! READ ELEMENT NODES
         DO I = 1, NUM_ELEMENT
            IF(.NOT. ALLOCATED(ELEMENTS(I)%NODE)) &
               CALL ERROR('The node of element ' // TRIM(NUM2STR(I,'I')) // ' has not been allocated!')
         ENDDO

         DO I = 1, NUM_ELEMENT
            COMMAND = GET_MACRO(CONTROL)
            IF(ICHAR(COMMAND(1:1)) .LT. 49 .OR. ICHAR(COMMAND(1:1)) .GT. 57) THEN
               CALL ERROR('Invalid element number: ' // COMMAND(1:INDEX(COMMAND,' ') - 1))
            ELSE
               READ(COMMAND,*) NO, (ELEMENTS(NO)%NODE(J), J = 1, ELEMENT_LIB(ELEMENTS(NO)%ELEMENT)%NUM_NODE)
            ENDIF
         ENDDO

      CASE ('ELEMENT_MATERIAL')
         ! READ ELEMENT MATERIAL
         COMMAND = GET_MACRO(CONTROL)
         DO WHILE(TRIM(COMMAND) .NE. 'END ELEMENT_MATERIAL')
            READ(COMMAND(INDEX(COMMAND, 'MATERIAL', BACK = .TRUE.) + 8:),*) NUMBER

            IF(INDEX(COMMAND, 'TO') .NE. 0) THEN
               ! A GROUP OF ELEMENTS ARE OF THAT TYPE OF MATERIAL
               READ(COMMAND,*) FIRST, STRING, LAST
               IF(LAST .LT. FIRST) CALL ERROR('Last element No. < first element No. ! (MATERIAL)')

               DO I = FIRST, LAST
                  ELEMENTS(I)%MATERIAL = NUMBER
               ENDDO
            ELSE
               ! A LIST OF ELEMENTS ARE OF THAT TYPE OF MATERIAL
               STRING = COMMAND(1:INDEX(COMMAND, 'MATERIAL') - 1)
               DO WHILE(LEN_TRIM(STRING) .GT. 0)
                  READ(STRING,*) I
                  ELEMENTS(I)%MATERIAL = NUMBER
                  STRING = REMOVE_FIRSTWORD(STRING)
               ENDDO
            ENDIF

            COMMAND = GET_MACRO(CONTROL)
         ENDDO

      CASE ('ELEMENT_GEOMETRY')
         ! READ ELEMENT TYPE
         COMMAND = GET_MACRO(CONTROL)
         DO WHILE(TRIM(COMMAND) .NE. 'END ELEMENT_GEOMETRY')
            READ(COMMAND(INDEX(COMMAND, 'GEOMETRY', BACK = .TRUE.) + 8:),*) NUMBER

            IF(INDEX(COMMAND, 'TO') .NE. 0) THEN
               ! A GROUP OF ELEMENTS ARE IN THAT TYPE OF GEOMETRY
               READ(COMMAND,*) FIRST, STRING, LAST
               IF(LAST .LT. FIRST) CALL ERROR('Last element No. < first element No. ! (GEOMETRY)')

               DO I = FIRST, LAST
                  ELEMENTS(I)%GEOMETRY = NUMBER
               ENDDO
            ELSE
               ! A LIST OF ELEMENTS ARE IN THAT TYPE OF GEOMETRY
               STRING = COMMAND(1:INDEX(COMMAND, 'GEOMETRY') - 1)
               DO WHILE(LEN_TRIM(STRING) .GT. 0)
                  READ(STRING,*) I
                  ELEMENTS(I)%GEOMETRY = NUMBER
                  STRING = REMOVE_FIRSTWORD(STRING)
               ENDDO
            ENDIF

            COMMAND = GET_MACRO(CONTROL)
         ENDDO

      CASE DEFAULT
         CALL ERROR('Invalid command: ' // TRIM(COMMAND) // ' in ELEMENT BLOCK')
   END SELECT
ENDDO

END SUBROUTINE READ_ELEMENTS