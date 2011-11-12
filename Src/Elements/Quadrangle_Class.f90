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
  PRIVATE :: init_faces_QUA_Q1
  PRIVATE :: init_faces_QUA_Q2

  PRIVATE :: rd_normal_fun
  PRIVATE :: rd_normal_QUA_Q1
  PRIVATE :: rd_normal_QUA_Q2

CONTAINS

  !================================================
  SUBROUTINE initialize_sub(e, mode, Nodes, Coords)
  !================================================

    IMPLICIT NONE

    CLASS(quadrangle)                        :: e
    CHARACTER(*),                 INTENT(IN) :: mode
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    !-------------------------------------------------

    INTEGER :: i, id
    !-------------------------------------------------

    SELECT CASE( SIZE(Nodes) )

    CASE(4)

       !---------
       ! ELELMENT
       !----------------------------------
       e%Type     = QUA_Q1 ! Element type
       e%N_verts  = 4      ! # Vertices
       e%N_points = 4      ! # DoFs

       IF( mode == "element") THEN
          ! FACES
          CALL init_faces_QUA_Q1(e, Nodes, Coords)
       ENDIF

    CASE(9)

       !---------
       ! ELELMENT
       !----------------------------------       
       e%Type     = QUA_Q2 ! Element type
       e%N_verts  = 4      ! # Vertices
       e%N_points = 9      ! # DoFs

       IF( mode == "element") THEN
          ! FACES
          CALL init_faces_QUA_Q2(e, Nodes, Coords)
       ENDIF

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

     !--------------------------
     ! Normals for the RD scheme
     !--------------------------------------
     ALLOCATE( e%rd_n(e%N_dim, e%N_points) )

     DO i = 1, e%N_points
        e%rd_n(:, i) = rd_normal_fun(e, i)
     ENDDO
     
  END SUBROUTINE initialize_sub
  !============================

  !========================================
  FUNCTION rd_normal_fun(e, i) RESULT(nn_i)
  !========================================

    IMPLICIT NONE

    CLASS(quadrangle)   :: e
    INTEGER, INTENT(IN) :: i

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_i
    !---------------------------------------    

    SELECT CASE(e%Type)

    CASE(QUA_Q1)

       nn_i = rd_normal_QUA_Q1(e, i)

    CASE(QUA_Q2)

       nn_i = rd_normal_QUA_Q2(e, i)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: Unsupported Quadrangle type'
       STOP

    END SELECT
    
  END FUNCTION rd_normal_fun
  !=========================  

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFIC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%    QUADRANGLE Q1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !=============================================
  SUBROUTINE init_faces_QUA_Q1(e, Nodes, Coords)
  !=============================================

    IMPLICIT NONE

    CLASS(quadrangle)                        :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    !-------------------------------------------------

    TYPE(segment), POINTER :: seg

    INTEGER,      DIMENSION(:),   ALLOCATABLE :: loc
    INTEGER,      DIMENSION(:),   ALLOCATABLE :: VV
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: RR

    INTEGER :: id, if, i, j, k, istat
    !-------------------------------------------------

    e%N_faces  = 4 ! # Faces

    ALLOCATE( e%faces(e%N_faces) )

    ! Dofs on the faces
    ALLOCATE( VV(2), loc(2) )

    ! Coordinates on the faces
    ALLOCATE( RR(SIZE(Coords, 1), 2) )

    DO if = 1, e%N_faces

       i = MOD(if, MIN(if+1, e%N_faces)) + 1
       j = MOD(i,  MIN(if+2, e%N_faces)) + 1

       loc = (/ i, j /)

       VV = Nodes(loc)
       RR = Coords(:, loc)

       ALLOCATE(seg, STAT=istat)          
       IF(istat /=0) THEN
          WRITE(*,*) 'ERROR: failed quadrangle allocation'
       ENDIF

       CALL seg%initialisize( loc, VV, RR )
       CALL seg%init_quadrature()
          
       e%faces(if)%f => seg
          
    ENDDO

    DEALLOCATE( VV, RR, loc )

  END SUBROUTINE init_faces_QUA_Q1
  !===============================

  !===========================================
  FUNCTION rd_normal_QUA_Q1(e, k) RESULT(nn_k)
  !===========================================

    IMPLICIT NONE

    CLASS(quadrangle)   :: e
    INTEGER, INTENT(IN) :: k

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_k
    !---------------------------------------

    INTEGER :: i, j
    !---------------------------------------

    SELECT CASE(k)
    CASE(1)
       i = 2; j = 4
    CASE(2)
       i = 3; j = 1
    CASE(3)
       i = 4; j = 2
    CASE(4)
       i = 1; j = 3
    CASE DEFAULT
       WRITE(*,*) 'ERROR: wrong nodes for rd_normal'
       STOP       
    END SELECT    

    nn_k(1) = e%Coords(2,i) - e%Coords(2,j)
    nn_k(2) = e%Coords(1,j) - e%Coords(1,i)

  END FUNCTION rd_normal_QUA_Q1
  !============================


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFIC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%    QUADRANGLE Q2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !=============================================
  SUBROUTINE init_faces_QUA_Q2(e, Nodes, Coords)
  !=============================================

    IMPLICIT NONE

    CLASS(quadrangle)                        :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    !-------------------------------------------------

    TYPE(segment), POINTER :: seg

    INTEGER,      DIMENSION(:),   ALLOCATABLE :: loc
    INTEGER,      DIMENSION(:),   ALLOCATABLE :: VV
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: RR

    INTEGER :: id, if, i, j, k, istat
    !-------------------------------------------------

    e%N_faces  = 4 ! # Faces

    ALLOCATE( e%faces(e%N_faces) )

    ! Dofs on the faces
    ALLOCATE( VV(3), loc(3) )

    ! Coordinates on the faces
    ALLOCATE( RR(SIZE(Coords, 1), 3) )

    DO if = 1, e%N_faces

       i = MOD(if, MIN(if+1, e%N_faces)) + 1
       j = MOD(i,  MIN(if+2, e%N_faces)) + 1

       k = i + 4

       loc = (/ i, j, k /)

       VV = Nodes(loc)
       RR = Coords(:, loc)

       ALLOCATE(seg, STAT=istat)          
       IF(istat /=0) THEN
          WRITE(*,*) 'ERROR: failed quadrangle allocation'
       ENDIF

       CALL seg%initialisize( loc, VV, RR )
       CALL seg%init_quadrature()

       e%faces(if)%f => seg
          
    ENDDO

    DEALLOCATE( VV, RR, loc )

  END SUBROUTINE init_faces_QUA_Q2
  !===============================

  !===========================================
  FUNCTION rd_normal_QUA_Q2(e, k) RESULT(nn_k)
  !===========================================

    IMPLICIT NONE

    CLASS(quadrangle)   :: e
    INTEGER, INTENT(IN) :: k

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_k
    !---------------------------------------

    INTEGER :: i, j
    !---------------------------------------

    SELECT CASE(k)
    CASE(1)
       i = 2; j = 4
    CASE(2)
       i = 3; j = 1
    CASE(3)
       i = 4; j = 2
    CASE(4)
       i = 1; j = 3
    CASE(5)
       i = 2; j = 1
    CASE(6)
       i = 3; j = 2
    CASE(7)
       i = 4; j = 3
    CASE(8)
       i = 1; j = 4
    CASE(9)
       i = 1; j = 1        
    CASE DEFAULT
       WRITE(*,*) 'ERROR: wrong nodes for rd_normal'
       STOP       
    END SELECT    

    nn_k(1) = e%Coords(2,i) - e%Coords(2,j)
    nn_k(2) = e%Coords(1,j) - e%Coords(1,i)

  END FUNCTION rd_normal_QUA_Q2
  !============================

END MODULE Quadrangle_Class

