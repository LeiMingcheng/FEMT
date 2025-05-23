!***************************************************************************************************
!- PURPORSE:
!     CHECK BASIC DATA IN INPUT FILE.
!  
!- INPUT ARGUMENTS:
!  NONE
!
!- OUTPUT ARGUMENTS:
!  NONE
!
!- CALL PROCEDURES:
!  MODULE PROCEDURE: GET_MACRO, ERROR, NUM2STR
!
!- CALLED BY:
!  INITIATE
!
!- PROGAMMED BY:
!  ZHIHAI XIANG, DEPARTMENT OF ENGINEERING MECHANICS, TSINGHUA UNIVERSITY, NOVEMBER 1, 2015
!***************************************************************************************************

SUBROUTINE CHECK_DATA

USE BASIC_DATA

IMPLICIT NONE

INTEGER(INT_KIND)                 I                ! LOOP INDEX
INTEGER(INT_KIND)                 J                ! LOOP INDEX
INTEGER(INT_KIND)                 K                ! LOOP INDEX
INTEGER(INT_KIND)::               FOUND_ERROR = 0  ! THE COUNTER OF FOUND ERRORS
INTEGER(INT_KIND), ALLOCATABLE :: FLAG(:)          ! ERROR FLAG


!-- ALLOCATE MEMORY FOR FLAG --
ALLOCATE(FLAG(MAX(NUM_NODE, NUM_ELEMENT)))

!-- CHECK IDENTICAL NODES --
FLAG(1:NUM_NODE) = 0

DO I = 1,NUM_NODE
    DO J =1,NUM_NODE
        IF (I == J) CYCLE
        IF (ALL(ABS(XYZ(:,I) - XYZ(:,J)) .LE. SQRT(REAL_EPS)))    THEN
            FLAG(I + J * 0) = J
            FLAG(J + I * 0) = I
        ENDIF
    ENDDO
ENDDO

IF(SUM(FLAG(1:NUM_NODE)) .GT. 0) THEN
   FOUND_ERROR = FOUND_ERROR + 1
   WRITE(*,*) 'Found identical nodes!'
   WRITE(LOG,*)
   WRITE(LOG, '(A, I1, A)') '*ERROR ', FOUND_ERROR, '-- Found the following identical node pairs:'
   DO I = 1, NUM_NODE
      IF(FLAG(I) .NE. 0) THEN
         WRITE(LOG, "( ' (', I0, 1X ',', I0')' )") I, FLAG(I)
         FLAG(FLAG(I)) = 0    ! ONLY LIST ONE TIME
      ENDIF
   ENDDO  
END IF

!-- CHECK USELESS NODES --
FLAG(1:NUM_NODE) = 1
DO I = 1, NUM_NODE
   DO J = 1, NUM_ELEMENT
      IF(ANY(ELEMENTS(J)%NODE .EQ. I)) THEN
         FLAG(I) = 0
         EXIT
      ENDIF
   ENDDO
ENDDO

IF(SUM(FLAG(1:NUM_NODE)) .GT. 0) THEN
   FOUND_ERROR = FOUND_ERROR + 1
   WRITE(*,*) 'Found useless nodes!'
   WRITE(LOG,*)
   WRITE(LOG, '(A, I1, A)') '*ERROR ', FOUND_ERROR, '-- The following nodes do not belong to any element:'
   DO I = 1, NUM_NODE
      IF(FLAG(I) .NE. 0) WRITE(LOG,*) I
   ENDDO  
END IF

!-- CHECK DUPLICATED NODES IN ONE ELEMENT --
FLAG(1:NUM_NODE) = 0
DO K = 1, NUM_ELEMENT
   FORALL(I = 1:ELEMENT_LIB(ELEMENTS(K)%ELEMENT)%NUM_NODE, J = 1:ELEMENT_LIB(ELEMENTS(K)%ELEMENT)%NUM_NODE, &
          I .NE. J .AND. ELEMENTS(K)%NODE(I) .EQ. ELEMENTS(K)%NODE(J))
      FLAG(ELEMENTS(K)%NODE(I) + J * 0) = K
   END FORALL
ENDDO

IF(SUM(FLAG(1:NUM_NODE)) .GT. 0) THEN
   FOUND_ERROR = FOUND_ERROR + 1
   WRITE(*,*) 'Found redundant nodes!'
   WRITE(LOG,*)
   WRITE(LOG, '(A, I1, A)') '*ERROR ', FOUND_ERROR, '-- Found the following duplicated nodes:'
   DO I = 1, NUM_NODE
      IF(FLAG(I) .NE. 0) WRITE(LOG, "(' Element ', I0, ', Node ', I0)") FLAG(I), I
   ENDDO  
END IF

!-- CHECK IDENTICAL ELEMENTS --
FLAG(1:NUM_ELEMENT) = 0

DO I = 1,NUM_ELEMENT
    DO J =1,NUM_ELEMENT
        IF (I == J) CYCLE
        IF (ELEMENTS(I) .EQ. ELEMENTS(J))    THEN
            FLAG(I + J * 0) = J
            FLAG(J + I * 0) = I
        ENDIF
    ENDDO
ENDDO

IF(SUM(FLAG(1:NUM_ELEMENT)) .GT. 0) THEN
   FOUND_ERROR = FOUND_ERROR + 1
   WRITE(*,*) 'Found identical elements!'
   WRITE(LOG,*)   
   WRITE(LOG, '(A, I1, A)') '*ERROR ', FOUND_ERROR, '-- Found the following identical element pairs:'
   DO I = 1, NUM_ELEMENT
      IF(FLAG(I) .NE. 0) THEN
         WRITE(LOG, "( ' (', I0, 1X ',', I0')' )") I, FLAG(I)
         FLAG(FLAG(I)) = 0    ! ONLY LIST ONE TIME
      ENDIF
   ENDDO  
END IF

!-- CHECK THE VALIDITY OF ELEMENT NODES --
FLAG(1:NUM_ELEMENT) = 1
DO K = 1, NUM_ELEMENT
   DO I = 1, ELEMENT_LIB(ELEMENTS(K)%ELEMENT)%NUM_NODE
      IF(ELEMENTS(K)%NODE(I) .LE. 0 .OR. ELEMENTS(K)%NODE(I) .GT. NUM_NODE) THEN
         FLAG(K) = ELEMENTS(K)%NODE(I)
         EXIT
      ENDIF
   ENDDO   
ENDDO

IF(SUM(FLAG(1:NUM_ELEMENT)) .NE. NUM_ELEMENT) THEN
   FOUND_ERROR = FOUND_ERROR + 1
   WRITE(*,*) 'Found invalid element nodes!'
   WRITE(LOG,*)   
   WRITE(LOG, '(A, I1, A)') '*ERROR ', FOUND_ERROR, '-- Found the following invalid element nodes:'
   DO I = 1, NUM_ELEMENT
      IF(FLAG(I) .NE. 1) WRITE(LOG, "(' Element ', I0, ', Node ', I0)") I, FLAG(I)
   ENDDO  
END IF

!-- CHECK ELEMENT MATERIAL NUMBER --
FLAG(1:NUM_ELEMENT) = 0
FORALL(I = 1:NUM_ELEMENT, ELEMENTS(I)%MATERIAL .LT. 1 .OR. ELEMENTS(I)%MATERIAL .GT. NUM_MATERIAL) FLAG(I) = 1
   
IF(SUM(FLAG(1:NUM_ELEMENT)) .GT. 0) THEN
   FOUND_ERROR = FOUND_ERROR + 1
   WRITE(*,*) 'Invalid element material number!'
   WRITE(LOG,*)   
   WRITE(LOG, '(A, I1, A)') '*ERROR ', FOUND_ERROR, '-- Found the following invalid element material numbers:'
   DO I = 1, NUM_ELEMENT
      IF(FLAG(I) .NE. 0) WRITE(LOG, "(' Element ', I0)") I
   ENDDO  
END IF

!-- CHECK GEOMETRY THICKNESS --
FLAG(1:NUM_GEOMETRY) = 0
DO I = 1,NUM_GEOMETRY  
   SELECT CASE(GEOMETRIES(I)% GEOMETRY_TYPE)
      CASE('PLANE')
         IF(GEOMETRIES(I)%THICKNESS .LE. 0.0_REAL_KIND) FLAG(I) = 1
      CASE DEFAULT
   END SELECT
ENDDO 
   
IF(SUM(FLAG(1:NUM_GEOMETRY)) .GT. 0) THEN
   FOUND_ERROR = FOUND_ERROR + 1
   WRITE(*,*) 'None-positive thickness!'
   WRITE(LOG,*)   
   WRITE(LOG, '(A, I1, A)') '*ERROR ', FOUND_ERROR, '-- Found the following geometries with non-positive thickness:'
   DO I = 1, NUM_GEOMETRY
      IF(FLAG(I) .NE. 0) WRITE(LOG, "(' Geometry ', I0)") I
   ENDDO  
END IF

!-- SUMMARY --
IF(FOUND_ERROR .GT. 0) &
   CALL ERROR(TRIM(NUM2STR(FOUND_ERROR, '(I0)')) // ' errors have been found during the data checking!')

DEALLOCATE(FLAG)

END SUBROUTINE CHECK_DATA