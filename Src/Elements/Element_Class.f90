MODULE Element_Class

  USE Face_Class
  
  IMPLICIT NONE

  !==================================  
  INTEGER, PARAMETER :: TRI_P1 = 10
  INTEGER, PARAMETER :: TRI_P2 = 11
  
  INTEGER, PARAMETER :: QUA_Q1 = 20
  INTEGER, PARAMETER :: QUA_Q2 = 21
  !==================================

  TYPE :: faces_ptr     
     CLASS(face), POINTER :: f     
  END TYPE faces_ptr
       
  !==================================
  TYPE, PUBLIC :: element

     INTEGER :: Type      ! Element type
     INTEGER :: N_dim = 2 ! # spatial dimension
     INTEGER :: N_verts   ! # vertices
     INTEGER :: N_points  ! # DoFs
     INTEGER :: N_Faces   ! # faces
     
     !REAL(KIND=8) :: Volume

     !------------------
     ! Element geometry
     !-------------------------------------------------
     ! Vertices coordinates
     ! Coords(id, j): id = 1, N_dim, j = 1, N_verts
     REAL(KIND=8), DIMENSION(:,:), POINTER :: Coords
     !-------------------------------------------------
     ! Nodes (global numeration of the mesh)
     ! NU(j), j = 1, N_points
     INTEGER,      DIMENSION(:),   POINTER :: NU
     ! Normals for the RD scheme
     ! rd_n(id, j), id = 1, N_dim, j = 1, N_points
     REAL(KIND=8), DIMENSION(:,:), POINTER :: rd_n

     !---------------
     !Faces structure
     !---------------------------------------------------
     TYPE(faces_ptr), DIMENSION(:), ALLOCATABLE :: faces
     
!!$   CONTAINS

  END TYPE element  
  !==================================
    
CONTAINS
     
END MODULE Element_Class
