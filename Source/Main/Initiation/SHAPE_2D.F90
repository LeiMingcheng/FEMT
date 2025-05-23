!***************************************************************************************************
!- PURPORSE:
!     CALCULATE THE SHAPE FUNCTIONS AND THEIR DERIVATIVES FOR 2D ELEMENTS.
!
!- INPUT ARGUMENTS:
!  INTE_POINT: THE INFORMATION OF INTEGRATION POINTS
!  NUM_POINT:  THE NUMBER OF INTEGRATION POINTS
!  NUM_NODE :  THE NUMBER OF ELEMENT NODES
!
!- OUTPUT ARGUMENTS:
!  SHAPES : THE SHAPE FUNCTION
!  D_SHAPE: THE DERIVATIVE OF SHAPE FUNCTION OVER NATURAL COORDINATES
!
!- CALL PROCEDURES:
!  NONE
!
!- CALLED BY
!  INIT_SHAPE
!
!- PROGAMMED BY:
!  ZHIHAI XIANG, DEPARTMENT OF ENGINEERING MECHANICS, TSINGHUA UNIVERSITY, JANUARY 16, 2016
!***************************************************************************************************

SUBROUTINE SHAPE_2D(INTE_POINT, NUM_POINT, NUM_NODE, SHAPES, D_SHAPE)

USE BASIC_DATA, ONLY: INT_KIND, REAL_KIND, NUM_INTE

IMPLICIT NONE

INTEGER(INT_KIND), INTENT(IN):: NUM_POINT                        ! THE NUMBER OF INTEGRATION POINTS
INTEGER(INT_KIND), INTENT(IN):: NUM_NODE                         ! THE NUMBER OF ELEMENT NODES
TYPE(NUM_INTE), INTENT(IN)::    INTE_POINT(NUM_POINT)            ! THE INFORMATION OF NUMERICAL INTEGRATION POINTS
REAL(REAL_KIND), INTENT(OUT)::  SHAPES(NUM_NODE, NUM_POINT)      ! THE SHAPE FUNCTION
REAL(REAL_KIND), INTENT(OUT)::  D_SHAPE(2, NUM_NODE, NUM_POINT)  ! THE DERIVATIVE OF SHAPE FUNCTION OVER NATURAL COORDINATES

! THE SIGN OF KESI, ETA FOR EACH ELEMENT NODE
REAL(REAL_KIND), PARAMETER::    NODE_SIGN(2,9) = RESHAPE( &
                                (/ &
                                 -1.0, -1.0,  &  ! Node 1: (-1,-1)
                                  1.0, -1.0,  &  ! Node 2: ( 1,-1)
                                  1.0,  1.0,  &  ! Node 3: ( 1, 1)
                                 -1.0,  1.0,  &  ! Node 4: (-1, 1)
                                  0.0, -1.0,  &  ! Node 5: ( 0,-1) (Q8 midside)
                                  1.0,  0.0,  &  ! Node 6: ( 1, 0) (Q8 midside)
                                  0.0,  1.0,  &  ! Node 7: ( 0, 1) (Q8 midside)
                                 -1.0,  0.0,  &  ! Node 8: (-1, 0) (Q8 midside)
                                  0.0,  0.0   &  ! Node 9: ( 0, 0) (Q9 center)
                                 /), &
                                (/2, 9/))

INTEGER(INT_KIND), PARAMETER::  CORNER_EDGE(2,4) = RESHAPE( &
                                (/5, 8,  5, 6,  6, 7,  7, 8/), &
                                (/2, 4/)) !
INTEGER(INT_KIND)               I                                ! LOOP INDEX FOR INTEGRATION POINTS
INTEGER(INT_KIND)               J                                ! LOOP INDEX FOR ELEMENT NODES
REAL(REAL_KIND)                 ksi, eta                         ! Natural coordinates for the current integration point
REAL(REAL_KIND)                 L_neg1_ksi, L_zero_ksi, L_pos1_ksi ! 1D Lagrange polys for ksi
REAL(REAL_KIND)                 L_neg1_eta, L_zero_eta, L_pos1_eta ! 1D Lagrange polys for eta
REAL(REAL_KIND)                 dL_neg1_ksi, dL_zero_ksi, dL_pos1_ksi ! Derivatives of 1D Lagrange polys for ksi
REAL(REAL_KIND)                 dL_neg1_eta, dL_zero_eta, dL_pos1_eta ! Derivatives of 1D Lagrange polys for eta


!-- SET SHAPE FUNCTIONS --
DO I = 1, NUM_POINT
   ksi = INTE_POINT(I)%COORD(1)
   eta = INTE_POINT(I)%COORD(2)

   ! TRIANGULAR ELEMENTS (AREA COORDINATES)
   IF(NUM_NODE .EQ. 3 .OR. NUM_NODE .EQ. 6) THEN   ! ���óɡ����޵�Ԫ����2003��� P108
      DO J = 1, 3
         SHAPES(J,I) = INTE_POINT(I)%COORD(J)
      ENDDO

      IF(NUM_NODE .EQ. 3) CYCLE

      SHAPES(4,I) = 4.0 * INTE_POINT(I)%COORD(1) * INTE_POINT(I)%COORD(2)
      SHAPES(5,I) = 4.0 * INTE_POINT(I)%COORD(2) * INTE_POINT(I)%COORD(3)
      SHAPES(6,I) = 4.0 * INTE_POINT(I)%COORD(3) * INTE_POINT(I)%COORD(1)

      SHAPES(1,I) = SHAPES(1,I) - (SHAPES(4,I) + SHAPES(6,I)) / 2.0
      SHAPES(2,I) = SHAPES(2,I) - (SHAPES(4,I) + SHAPES(5,I)) / 2.0
      SHAPES(3,I) = SHAPES(3,I) - (SHAPES(5,I) + SHAPES(6,I)) / 2.0
   ELSE ! RECTANGULAR ELEMENTS (CARTESIAN COORDINATES) ���óɡ����޵�Ԫ����2003��� P114
      ! NODES: (1,2,3,4) for Q4
      DO J = 1, 4
         SHAPES(J,I) = (1.0 + NODE_SIGN(1,J) * ksi) * &
                       (1.0 + NODE_SIGN(2,J) * eta) / 4.0
      ENDDO

      IF(NUM_NODE .EQ. 4) CYCLE

      ! This section is for 8-node Serendipity element (Q8).
      ! NODE_SIGN array indices 1-8 now correctly correspond to Q8 node numbering (corners then midsides).
      ! Index 9 is for the Q9 center node.
      IF (NUM_NODE .EQ. 8) THEN
          ! Midside nodes 5,7 (on ksi=0 virtual lines, using NODE_SIGN for eta direction)
          DO J = 5, 7, 2 ! Nodes 5 (0,-1), 7 (0,1)
             SHAPES(J,I) = (1.0 + NODE_SIGN(2,J) * eta) / 2.0 * &
                           (1.0 - ksi ** 2)
          ENDDO

          ! Midside nodes 6,8 (on eta=0 virtual lines, using NODE_SIGN for ksi direction)
          DO J = 6, 8, 2 ! Nodes 6 (1,0), 8 (-1,0)
             SHAPES(J,I) = (1.0 + NODE_SIGN(1,J) * ksi) / 2.0 * &
                           (1.0 - eta ** 2)
          ENDDO

          ! ADJUST THE CORNER NODES for Q8
          DO J = 1, 4
             SHAPES(J,I) = SHAPES(J,I) - (SHAPES(CORNER_EDGE(1,J),I) + SHAPES(CORNER_EDGE(2,J),I)) / 2.0
          ENDDO
      ENDIF


      ! FOR 9-NODE LAGRANGIAN ELEMENT (Q9)
      IF(NUM_NODE .EQ. 9) THEN
         ! Define 1D Quadratic Lagrangian Polynomials and their derivatives
         ! For nodes at zeta = -1, 0, 1
         L_neg1_ksi = 0.5 * ksi * (ksi - 1.0)
         L_zero_ksi = (1.0 - ksi**2)
         L_pos1_ksi = 0.5 * ksi * (ksi + 1.0)

         L_neg1_eta = 0.5 * eta * (eta - 1.0)
         L_zero_eta = (1.0 - eta**2)
         L_pos1_eta = 0.5 * eta * (eta + 1.0)

         dL_neg1_ksi = 0.5 * (2.0 * ksi - 1.0)
         dL_zero_ksi = -2.0 * ksi
         dL_pos1_ksi = 0.5 * (2.0 * ksi + 1.0)

         dL_neg1_eta = 0.5 * (2.0 * eta - 1.0)
         dL_zero_eta = -2.0 * eta
         dL_pos1_eta = 0.5 * (2.0 * eta + 1.0)

         ! Corner Nodes (1-4) - Standard Lagrangian
         ! Node 1: (-1,-1)
         SHAPES(1,I) = L_neg1_ksi * L_neg1_eta
         ! Node 2: (1,-1)
         SHAPES(2,I) = L_pos1_ksi * L_neg1_eta
         ! Node 3: (1,1)
         SHAPES(3,I) = L_pos1_ksi * L_pos1_eta
         ! Node 4: (-1,1)
         SHAPES(4,I) = L_neg1_ksi * L_pos1_eta

         ! Mid-side Nodes (5-8) - Standard Lagrangian
         ! Node 5: (0,-1)
         SHAPES(5,I) = L_zero_ksi * L_neg1_eta
         ! Node 6: (1,0)
         SHAPES(6,I) = L_pos1_ksi * L_zero_eta
         ! Node 7: (0,1)
         SHAPES(7,I) = L_zero_ksi * L_pos1_eta
         ! Node 8: (-1,0)
         SHAPES(8,I) = L_neg1_ksi * L_zero_eta

         ! Center node (9) - Standard Lagrangian
         SHAPES(9,I) = L_zero_ksi * L_zero_eta ! Equivalent to (1.0 - ksi**2) * (1.0 - eta**2)
      ENDIF
   ENDIF
ENDDO

!-- SET THE DERIVATIVE OF SHAPE FUNCTIONS --
DO I = 1, NUM_POINT
   ksi = INTE_POINT(I)%COORD(1)
   eta = INTE_POINT(I)%COORD(2)

   ! TRIANGULAR ELEMENTS (AREA COORDINATES)
   IF(NUM_NODE .EQ. 3 .OR. NUM_NODE .EQ. 6) THEN
      D_SHAPE(1,1,I) =  1.0
      D_SHAPE(2,1,I) =  0.0
      D_SHAPE(1,2,I) =  0.0
      D_SHAPE(2,2,I) =  1.0
      D_SHAPE(1,3,I) = -1.0
      D_SHAPE(2,3,I) = -1.0

      IF(NUM_NODE .EQ. 3) CYCLE

      D_SHAPE(1,4,I) =  4.0 * INTE_POINT(I)%COORD(2)
      D_SHAPE(2,4,I) =  4.0 * INTE_POINT(I)%COORD(1)
      D_SHAPE(1,5,I) = -4.0 * INTE_POINT(I)%COORD(2)
      D_SHAPE(2,5,I) =  4.0 * (INTE_POINT(I)%COORD(3) - INTE_POINT(I)%COORD(2))
      D_SHAPE(1,6,I) =  4.0 * (INTE_POINT(I)%COORD(3) - INTE_POINT(I)%COORD(1))
      D_SHAPE(2,6,I) = -4.0 * INTE_POINT(I)%COORD(1)

      D_SHAPE(:,1,I) = D_SHAPE(:,1,I) - (D_SHAPE(:,4,I) + D_SHAPE(:,6,I)) / 2.0
      D_SHAPE(:,2,I) = D_SHAPE(:,2,I) - (D_SHAPE(:,4,I) + D_SHAPE(:,5,I)) / 2.0
      D_SHAPE(:,3,I) = D_SHAPE(:,3,I) - (D_SHAPE(:,5,I) + D_SHAPE(:,6,I)) / 2.0
   ELSE ! RECTANGULAR ELEMENTS (CARTESIAN COORDINATES)
      ! NODES: (1,2,3,4) for Q4
      DO J = 1, 4
         D_SHAPE(1,J,I) = NODE_SIGN(1,J) * (1.0 + NODE_SIGN(2,J) * eta) / 4.0
         D_SHAPE(2,J,I) = NODE_SIGN(2,J) * (1.0 + NODE_SIGN(1,J) * ksi) / 4.0
      ENDDO

      IF(NUM_NODE .EQ. 4) CYCLE

      ! This section is for 8-node Serendipity element (Q8)
      IF (NUM_NODE .EQ. 8) THEN
          ! Midside nodes 5,7 (on ksi=0 virtual lines)
          DO J = 5, 7, 2 ! Nodes 5 (0,-1), 7 (0,1)
             D_SHAPE(1,J,I) = -ksi * (1.0 + NODE_SIGN(2,J) * eta) ! Simplified from -(1.0 + NODE_SIGN(2,J) * eta) * ksi
             D_SHAPE(2,J,I) = NODE_SIGN(2,J) * (1.0 - ksi ** 2)/ 2.0
          ENDDO

          ! Midside nodes 6,8 (on eta=0 virtual lines)
          DO J = 6, 8, 2 ! Nodes 6 (1,0), 8 (-1,0)
             D_SHAPE(1,J,I) = NODE_SIGN(1,J) * (1.0 - eta ** 2)/ 2.0
             D_SHAPE(2,J,I) = -eta * (1.0 + NODE_SIGN(1,J) * ksi) ! Simplified from -(1.0 + NODE_SIGN(1,J) * ksi) * eta
          ENDDO

          ! ADJUST THE CORNER NODES for Q8
          DO J = 1, 4
             D_SHAPE(:,J,I) = D_SHAPE(:,J,I) - (D_SHAPE(:,CORNER_EDGE(1,J),I) + D_SHAPE(:,CORNER_EDGE(2,J),I)) / 2.0
          ENDDO
      ENDIF

      ! FOR 9-NODE LAGRANGIAN ELEMENT (Q9)
      IF(NUM_NODE .EQ. 9) THEN
         ! Define 1D Quadratic Lagrangian Polynomials and their derivatives again for safety / clarity
         ! For nodes at zeta = -1, 0, 1
         L_neg1_ksi = 0.5 * ksi * (ksi - 1.0)
         L_zero_ksi = (1.0 - ksi**2)
         L_pos1_ksi = 0.5 * ksi * (ksi + 1.0)

         L_neg1_eta = 0.5 * eta * (eta - 1.0)
         L_zero_eta = (1.0 - eta**2)
         L_pos1_eta = 0.5 * eta * (eta + 1.0)

         dL_neg1_ksi = 0.5 * (2.0 * ksi - 1.0)
         dL_zero_ksi = -2.0 * ksi
         dL_pos1_ksi = 0.5 * (2.0 * ksi + 1.0)

         dL_neg1_eta = 0.5 * (2.0 * eta - 1.0)
         dL_zero_eta = -2.0 * eta
         dL_pos1_eta = 0.5 * (2.0 * eta + 1.0)

         ! Corner Nodes (1-4) - Standard Lagrangian Derivatives
         ! Node 1: (-1,-1)
         D_SHAPE(1,1,I) = dL_neg1_ksi * L_neg1_eta
         D_SHAPE(2,1,I) = L_neg1_ksi * dL_neg1_eta
         ! Node 2: (1,-1)
         D_SHAPE(1,2,I) = dL_pos1_ksi * L_neg1_eta
         D_SHAPE(2,2,I) = L_pos1_ksi * dL_neg1_eta
         ! Node 3: (1,1)
         D_SHAPE(1,3,I) = dL_pos1_ksi * L_pos1_eta
         D_SHAPE(2,3,I) = L_pos1_ksi * dL_pos1_eta
         ! Node 4: (-1,1)
         D_SHAPE(1,4,I) = dL_neg1_ksi * L_pos1_eta
         D_SHAPE(2,4,I) = L_neg1_ksi * dL_pos1_eta

         ! Mid-side Nodes (5-8) - Standard Lagrangian Derivatives
         ! Node 5: (0,-1)
         D_SHAPE(1,5,I) = dL_zero_ksi * L_neg1_eta
         D_SHAPE(2,5,I) = L_zero_ksi * dL_neg1_eta
         ! Node 6: (1,0)
         D_SHAPE(1,6,I) = dL_pos1_ksi * L_zero_eta
         D_SHAPE(2,6,I) = L_pos1_ksi * dL_zero_eta
         ! Node 7: (0,1)
         D_SHAPE(1,7,I) = dL_zero_ksi * L_pos1_eta
         D_SHAPE(2,7,I) = L_zero_ksi * dL_pos1_eta
         ! Node 8: (-1,0)
         D_SHAPE(1,8,I) = dL_neg1_ksi * L_zero_eta
         D_SHAPE(2,8,I) = L_neg1_ksi * dL_zero_eta

         ! Center node (9) - Standard Lagrangian Derivatives
         D_SHAPE(1,9,I) = dL_zero_ksi * L_zero_eta ! Equivalent to -2.0 * ksi * (1.0 - eta**2)
         D_SHAPE(2,9,I) = L_zero_ksi * dL_zero_eta ! Equivalent to -2.0 * eta * (1.0 - ksi**2)
      ENDIF
   ENDIF
ENDDO

END SUBROUTINE SHAPE_2D