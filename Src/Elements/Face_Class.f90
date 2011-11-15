MODULE Face_Class

  IMPLICIT NONE

  !==================================
  INTEGER, PARAMETER :: SEG_P1 = 1
  INTEGER, PARAMETER :: SEG_P2 = 2
  !==================================

  !==================================
  TYPE, PUBLIC :: face

     INTEGER :: Type
     INTEGER :: N_dim = 2
     INTEGER :: N_verts
     INTEGER :: N_points

     REAL(KIND=8) :: Area

     !-------------------------------

     REAL(KIND=8), DIMENSION(:,:), POINTER :: Coords
     INTEGER,      DIMENSION(:),   POINTER :: NU
     INTEGER,      DIMENSION(:),   POINTER :: l_nu
     INTEGER                               :: g_seg
     INTEGER,      DIMENSION(:),   POINTER :: c_ele
     
     !-------------------------------
     
     INTEGER :: N_quad

     REAL(KIND=8), DIMENSION(:),     POINTER :: w_q
     REAL(KIND=8), DIMENSION(:,:),   POINTER :: n_q
     REAL(KIND=8), DIMENSION(:,:),   POINTER :: phi_q
     REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_q
          

!   CONTAINS

  END TYPE face
  !=================================

CONTAINS

END MODULE Face_Class
