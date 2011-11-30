MODULE Quadrature_rules

  USE Element_Class
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: oInt_n, Int_d

CONTAINS

  !====================================
  FUNCTION oInt_n(ele, ff) RESULT(i_ff)
  !====================================
  !
  ! \oint_{\partial E} {ff * n}
  !
    IMPLICIT NONE

    TYPE(element),                INTENT(IN) :: ele
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ff

    REAL(KIND=8) :: i_ff
    !----------------------------------------------

    INTEGER,                      POINTER :: N_quad
    INTEGER,                      POINTER :: N_points
    INTEGER,      DIMENSION(:),   POINTER :: loc
    REAL(KIND=8), DIMENSION(:,:), POINTER :: p
    REAL(KIND=8), DIMENSION(:,:), POINTER :: n
    REAL(KIND=8), DIMENSION(:),   POINTER :: w
    
    
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ff_q

    INTEGER :: j, iq, k
    !----------------------------------------------

    i_ff = 0.d0

    ALLOCATE( ff_q(SIZE(ff,1)) )

    DO j = 1, ele%N_faces

       N_quad   => ele%faces(j)%f%N_quad
       N_points => ele%faces(j)%f%N_points
       loc      => ele%faces(j)%f%l_nu
       p        => ele%faces(j)%f%phi_q
       w        => ele%faces(j)%f%w_q
       n        => ele%faces(j)%f%n_q
       
       DO iq = 1, N_quad

          ff_q = 0.d0
          DO k = 1, N_points 

             ff_q = ff_q + ff(:, loc(k)) * p(k, iq)

          ENDDO 

          i_ff = i_ff + w(iq) * DOT_PRODUCT(ff_q, n(:, iq) )

       ENDDO
       
       NULLIFY(N_quad, N_points, loc, p, w, n)
       
    ENDDO

    DEALLOCATE( ff_q )

  END FUNCTION oInt_n
  !==================  

  !===================================
  FUNCTION Int_d(ele, ff) RESULT(i_ff)
  !===================================
  !
  ! \int_{E} {ff}
  !
    IMPLICIT NONE

    TYPE(element),              INTENT(IN) :: ele
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: ff

    REAL(KIND=8) :: i_ff
    !----------------------------------------------

    INTEGER :: N_quad,  N_points
   
    REAL(KIND=8), DIMENSION(:,:), POINTER :: p
    REAL(KIND=8), DIMENSION(:),   POINTER :: w
        
    REAL(KIND=8) :: ff_q

    INTEGER :: iq
    !----------------------------------------------

    i_ff = 0.d0

    N_quad   =  ele%N_quad
    N_points =  ele%N_points
    p        => ele%phi_q
    w        => ele%w_q

    DO iq = 1, N_quad

       ff_q = SUM( ff * p(:, iq) )

       i_ff = i_ff + w(iq) * ff_q

    ENDDO
       
    NULLIFY( p, w )
       
  END FUNCTION Int_d
  !=================

END MODULE Quadrature_rules
