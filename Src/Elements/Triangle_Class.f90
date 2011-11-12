MODULE Trianlge_Class

  USE Element_Class
  USE Segment_Class

  IMPLICIT NONE

  !=========================================
  TYPE, PUBLIC, EXTENDS(element) :: triangle

     ! < >

   CONTAINS

     PROCEDURE, PUBLIC :: initialize => initialize_sub
    
  END TYPE triangle
  !=========================================

  PRIVATE :: initialize_sub
  PRIVATE :: init_faces_TRI_P1
  PRIVATE :: init_faces_TRI_P2

  PRIVATE :: rd_normal_fun
  PRIVATE :: rd_normal_TRI_P1
  PRIVATE :: rd_normal_TRI_P2

CONTAINS

  !================================================
  SUBROUTINE initialize_sub(e, mode, Nodes, Coords)
  !================================================

    IMPLICIT NONE

    CLASS(triangle)                          :: e
    CHARACTER(*),                 INTENT(IN) :: mode
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    !-------------------------------------------------

    INTEGER :: i, id
    !-------------------------------------------------

    SELECT CASE( SIZE(Nodes) )

    CASE(3)

       !---------
       ! ELELMENT
       !----------------------------------
       e%Type     = TRI_P1 ! Element type
       e%N_verts  = 3      ! # Vertices 
       e%N_points = 3      ! # DoFs       

       IF( mode == "element") THEN
          ! FACES
          CALL init_faces_TRI_P1(e, Nodes, Coords)
       ENDIF
       
    CASE(6)

       !---------
       ! ELELMENT
       !----------------------------------
       e%Type     = TRI_P2 ! Element type
       e%N_verts  = 3      ! # Vertices 
       e%N_points = 6      ! # DoFs

       IF( mode == "element") THEN
          ! FACES
          CALL init_faces_TRI_P2(e, Nodes, Coords)
       ENDIF
   
     CASE DEFAULT

        WRITE(*,*) 'ERROR: Unsupported Triangle type'
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

    CLASS(triangle)     :: e
    INTEGER, INTENT(IN) :: i

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_i
    !---------------------------------------    

    SELECT CASE(e%Type)

    CASE(TRI_P1)

       nn_i = rd_normal_TRI_P1(e, i)

    CASE(TRI_P2)

       nn_i = rd_normal_TRI_P2(e, i)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: Unsupported Triangle type'
       STOP

    END SELECT
    
  END FUNCTION rd_normal_fun
  !=========================  

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFIC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%     TRIANGLE P1    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !=============================================
  SUBROUTINE init_faces_TRI_P1(e, Nodes, Coords)
  !=============================================

    IMPLICIT NONE

    CLASS(triangle)                          :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    !-------------------------------------------------

     TYPE(segment), POINTER :: seg

    INTEGER,      DIMENSION(:),   ALLOCATABLE :: loc
    INTEGER,      DIMENSION(:),   ALLOCATABLE :: VV
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: RR

    INTEGER :: id, if, i, j, k, istat
    !-------------------------------------------------

    e%N_faces  = 3 ! # Faces
       
    ALLOCATE( e%faces(e%N_faces) )
    
    ! DoFs on the faces
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
          WRITE(*,*) 'ERROR: failed trianlge faces allocation'
       ENDIF

       CALL seg%initialisize( loc, VV, RR )
       CALL seg%init_quadrature()

       e%faces(if)%f => seg
          
    ENDDO

    DEALLOCATE( VV, RR, loc )

  END SUBROUTINE init_faces_TRI_P1
  !===============================

  !===========================================
  FUNCTION rd_normal_TRI_P1(e, k) RESULT(nn_k)
  !===========================================

    IMPLICIT NONE

    CLASS(triangle)     :: e
    INTEGER, INTENT(IN) :: k

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_k
    !---------------------------------------

    INTEGER :: i, j
    !---------------------------------------

    i = e%faces(k)%f%l_NU(1)
    j = e%faces(k)%f%l_NU(2)

    nn_k(1) = e%Coords(2,i) - e%Coords(2,j)
    nn_k(2) = e%Coords(1,j) - e%Coords(1,i)

  END FUNCTION rd_normal_TRI_P1
  !============================

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFIC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%     TRIANGLE P2    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  !=============================================
  SUBROUTINE init_faces_TRI_P2(e, Nodes, Coords)
  !=============================================

    IMPLICIT NONE

    CLASS(triangle)                          :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    !-------------------------------------------------

     TYPE(segment), POINTER :: seg

    INTEGER,      DIMENSION(:),   ALLOCATABLE :: loc
    INTEGER,      DIMENSION(:),   ALLOCATABLE :: VV
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: RR

    INTEGER :: id, if, i, j, k, istat
    !-------------------------------------------------

    e%N_faces  = 3 ! # Faces

    ALLOCATE( e%faces(e%N_faces) )

    ! Dofs on the faces
    ALLOCATE( VV(3), loc(3) )

    ! Coordinates on the faces
    ALLOCATE( RR(SIZE(Coords, 1), 3) )

    DO if = 1, e%N_faces

       i = MOD(if, MIN(if+1, e%N_faces)) + 1          
       j = MOD(i,  MIN(if+2, e%N_faces)) + 1

       k = i + 3

       loc = (/ i, j, k /)

       VV = Nodes(loc)
       RR = Coords(:, loc)

       ALLOCATE(seg, STAT=istat)          
       IF(istat /=0) THEN
          WRITE(*,*) 'ERROR: failed trianlge allocation'
       ENDIF

       CALL seg%initialisize( loc, VV, RR )
       CALL seg%init_quadrature()
         
       e%faces(if)%f => seg
          
    ENDDO

    DEALLOCATE( VV, RR, loc)

  END SUBROUTINE init_faces_TRI_P2
  !===============================

  !===========================================
  FUNCTION rd_normal_TRI_P2(e, k) RESULT(nn_k)
  !===========================================

    IMPLICIT NONE

    CLASS(triangle)     :: e
    INTEGER, INTENT(IN) :: k

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_k
    !---------------------------------------

    INTEGER :: i, j, kk
    !---------------------------------------

    SELECT CASE(k)

    CASE(1, 2, 3)

       i = e%faces(k)%f%l_NU(1)
       j = e%faces(k)%f%l_NU(2)

       nn_k(1) = e%Coords(2,i) - e%Coords(2,j)
       nn_k(2) = e%Coords(1,j) - e%Coords(1,i)       

    CASE(4, 5, 6)

       kk = MOD(k+1, MIN(k+2, e%N_faces)) + 1

       i = e%faces(kk)%f%l_NU(1)      
       j = e%faces(kk)%f%l_NU(2)
       
       nn_k(1) = e%Coords(2,i) - e%Coords(2,j)       
       nn_k(2) = e%Coords(1,j) - e%Coords(1,i)

       nn_k = -nn_k       

    CASE DEFAULT

       WRITE(*,*) 'ERROR: wrong nodes for rd_normal'
       STOP

    END SELECT
       
  END FUNCTION rd_normal_TRI_P2
  !============================

END MODULE Trianlge_Class
