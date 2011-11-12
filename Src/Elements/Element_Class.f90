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

     INTEGER :: Type
     INTEGER :: N_dim = 2
     INTEGER :: N_verts
     INTEGER :: N_points
     INTEGER :: N_Faces
     
     REAL(KIND=8) :: Volume

     !-------------------------------
     
     REAL(KIND=8), DIMENSION(:,:), POINTER :: Coords
     INTEGER,      DIMENSION(:),   POINTER :: NU

     REAL(KIND=8), DIMENSION(:,:), POINTER :: rd_n

     !-------------------------------
     
     TYPE(faces_ptr), DIMENSION(:), ALLOCATABLE :: faces
     
!!$   CONTAINS

  END TYPE element  
  !==================================
    
CONTAINS
     
END MODULE Element_Class
