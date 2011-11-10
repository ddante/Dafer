MODULE Quadrangle_Class

  USE Element_Class
  USE Segment_Class

  IMPLICIT NONE

  !=========================================
  TYPE, PUBLIC, EXTENDS(element) :: quadrangle

     ! < >

   CONTAINS

     PROCEDURE, PUBLIC :: initialize => initialize_sub

  END TYPE quadrangle
  !=========================================

  PRIVATE :: initialize_sub

CONTAINS

  !==========================================
  SUBROUTINE initialize_sub(e, Nodes, Coords)
  !==========================================

    IMPLICIT NONE

    CLASS(quadrangle)                        :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    !-------------------------------------------------

    TYPE(segment), POINTER :: seg

    INTEGER,      DIMENSION(:),   ALLOCATABLE :: VV
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: RR

    INTEGER :: id, i, j, k, istat
    !-------------------------------------------------

    SELECT CASE( SIZE(Nodes) )

    CASE(4)

       !---------
       ! ELELMENT
       !----------------------------------
       e%Type     = QUA_Q1 ! Element type
       e%N_verts  = 4      ! # Vertices
       e%N_points = 4      ! # DoFs
       e%N_faces  = 4      ! # Faces
       e%Type_f   = SEG_P1 ! Faces type

       !------
       ! FACES
       !----------------------------------
       ALLOCATE( e%faces(e%N_faces) )

       ! Dofs on the faces
       ALLOCATE( VV(2) )

       ! Coordinates on the faces
       ALLOCATE( RR(SIZE(Coords, 1), 2) )

       DO i = 1, e%N_faces

          j = MOD(i, MIN(i+1, e%N_verts)) + 1

          VV = (/ Nodes(i), Nodes(j) /)
          RR = Coords( :, (/ i, j /) )

          ALLOCATE(seg, STAT=istat)          
          IF(istat /=0) THEN
             WRITE(*,*) 'ERROR: failed trianlge allocation'
          ENDIF

          CALL seg%initialisize( VV, RR )

          e%faces(i)%f => seg
          
       ENDDO

       DEALLOCATE( VV, RR )

    CASE(9)

       !---------
       ! ELELMENT
       !----------------------------------       
       e%Type     = QUA_Q2 ! Element type
       e%N_verts  = 4      ! # Vertices
       e%N_points = 9      ! # DoFs
       e%N_faces  = 4      ! # Faces
       e%Type_f   = SEG_P2 ! Faces type

       !------
       ! FACES
       !----------------------------------
       ALLOCATE( e%faces(e%N_faces) )

       ! Dofs on the faces
       ALLOCATE( VV(3) )

       ! Coordinates on the faces
       ALLOCATE( RR(SIZE(Coords, 1), 3) )

       DO i = 1, e%N_faces

          j = MOD(i, MIN(i+1, e%N_verts)) + 1

          k = i + 4

          VV = (/ Nodes(i), Nodes(j), Nodes(k) /)
          RR = Coords( :, (/ i, j, k /) )

          ALLOCATE(seg, STAT=istat)          
          IF(istat /=0) THEN
             WRITE(*,*) 'ERROR: failed trianlge allocation'
          ENDIF

          CALL seg%initialisize( VV, RR )

          e%faces(i)%f => seg
          
       ENDDO

       DEALLOCATE( VV, RR )

     CASE DEFAULT

        WRITE(*,*) 'ERROR: Unsupported Quadrangle type'
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
  !============================
  
END MODULE Quadrangle_Class

