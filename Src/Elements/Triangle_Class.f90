MODULE Trianlge_Class

  USE Element_Class

  IMPLICIT NONE

  !=========================================
  TYPE, PUBLIC, EXTENDS(element) :: triangle

     ! < >

   CONTAINS

     PROCEDURE, PUBLIC :: initilisize => initilize_sub

  END TYPE triangle
  !=========================================

  PRIVATE :: initilize_sub

CONTAINS

  !=========================================
  SUBROUTINE initilize_sub(e, Nodes, Coords)
  !=========================================

    IMPLICIT NONE

    CLASS(triangle)                          :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    !-------------------------------------------------

    INTEGER :: id
    !-------------------------------------------------

    SELECT CASE( SIZE(Nodes) )

    CASE(3)

       e%Type     = TRI_P1
       e%N_verts  = 3
       e%N_points = 3
       e%N_faces  = 3
       e%Type_f   = SEG_P1

    CASE(6)

       e%Type     = TRI_P2
       e%N_verts  = 3
       e%N_points = 6
       e%N_faces  = 3
       e%Type_f   = SEG_P2

     CASE DEFAULT

        WRITE(*,*) 'ERROR: Unsupported Triangle type'
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
  
END MODULE Trianlge_Class
