MODULE Quadrangle_Class

  USE Element_Class

  IMPLICIT NONE

  !=========================================
  TYPE, PUBLIC, EXTENDS(element) :: quadrangle

     ! < >

   CONTAINS

     PROCEDURE, PUBLIC :: initilisize => initilize_sub

  END TYPE quadrangle
  !=========================================

  PRIVATE :: initilize_sub

CONTAINS

  !=========================================
  SUBROUTINE initilize_sub(e, Nodes, Coords)
  !=========================================

    IMPLICIT NONE

    CLASS(quadrangle)                          :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    !-------------------------------------------------

    INTEGER :: id
    !-------------------------------------------------

    SELECT CASE( SIZE(Nodes) )

    CASE(4)

       e%Type     = QUA_Q1
       e%N_verts  = 4
       e%N_points = 4
       e%N_faces  = 5
       e%Type_f   = SEG_P1

    CASE(9)

       e%Type     = QUA_Q2
       e%N_verts  = 4
       e%N_points = 9
       e%N_faces  = 5
       e%Type_f   = SEG_P2

     CASE DEFAULT

        WRITE(*,*) 'ERROR: Unsupported Quadrangle type'
        STOP

     END SELECT

     ALLOCATE( e%Coords(SIZE(Coords,1), SIZE(Coords,2)) )
     ALLOCATE( e%NU(SIZE(Nodes)) )
     
     DO id = 1, SIZE(Coords,1)
        e%Coords(id, :) = Coords(id, :)
     ENDDO

     e%NU = Nodes
     
  END SUBROUTINE initilize_sub
  !===========================
  
END MODULE Quadrangle_Class

