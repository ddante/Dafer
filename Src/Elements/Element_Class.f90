MODULE Element_Class

  IMPLICIT NONE

  !==================================  
  INTEGER, PARAMETER :: SEG_P1 = 1
  INTEGER, PARAMETER :: SEG_P2 = 2
  
  INTEGER, PARAMETER :: TRI_P1 = 10
  INTEGER, PARAMETER :: TRI_P2 = 11
  
  INTEGER, PARAMETER :: QUA_Q1 = 20
  INTEGER, PARAMETER :: QUA_Q2 = 21
  !==================================

  !==================================
  TYPE, PUBLIC :: element

     INTEGER :: Type
     INTEGER :: N_dim = 2
     INTEGER :: N_verts
     INTEGER :: N_points
     INTEGER :: N_Faces
     INTEGER :: Type_f

     REAL(KIND=8), DIMENSION(:,:), POINTER :: Coords
     INTEGER,      DIMENSION(:),   POINTER :: NU

!!$   CONTAINS

  END type element     
  !==================================
    
CONTAINS
     
END MODULE Element_Class
