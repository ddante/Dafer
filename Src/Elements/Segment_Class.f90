MODULE Segment_Class

  USE Face_Class

  IMPLICIT NONE
  
  !=========================================
  TYPE, PUBLIC, EXTENDS(face) :: segment

     ! < >

   CONTAINS

     PROCEDURE, PUBLIC :: initialisize => initialize_sub

  END TYPE segment
  !=========================================

  PRIVATE :: initialize_sub

CONTAINS

  !=========================================
  SUBROUTINE initialize_sub(e, Nodes, Coords)
  !=========================================

    IMPLICIT NONE

    CLASS(segment)                           :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    !-------------------------------------------------

    INTEGER :: id
    !-------------------------------------------------

    SELECT CASE( SIZE(Nodes) )

    CASE(2)

       e%N_verts  = 2
       e%N_points = 2

    CASE(3)

       e%N_verts  = 2
       e%N_points = 3

     CASE DEFAULT

        WRITE(*,*) 'ERROR: Unsupported Segment type'
        STOP

     END SELECT

     !---------------------------
     ! Attach data to the element
     !---------------------------------------------------
     ALLOCATE( e%Coords(SIZE(Coords,1), SIZE(Coords,2)) )
     ALLOCATE( e%NU(SIZE(Nodes)) )
     
     DO id = 1, SIZE(Coords,1)
        e%Coords(id, :) = Coords(id, :)
     ENDDO

     e%NU = Nodes
     
  END SUBROUTINE initialize_sub
  !===========================

END MODULE Segment_Class

