MODULE Face_Class

  IMPLICIT NONE

  !==================================
  TYPE, PUBLIC :: face

     INTEGER :: Type
     INTEGER :: N_verts
     INTEGER :: N_points

     REAL(KIND=8), DIMENSION(:,:), POINTER :: Coords
     INTEGER,      DIMENSION(:),   POINTER :: NU
integer :: ax
  END TYPE face
  !=================================

CONTAINS

END MODULE Face_Class
