MODULE LLxFS_method

  USE Element_class
  USE geometry,       ONLY: N_dim, elements
  USE init_problem,   ONLY: pb_type
  USE models,         ONLY: advection_speed
  
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: LLxFS_scheme

CONTAINS

  !=====================================================
  SUBROUTINE LLxFS_scheme(ele, Phi_tot, u, Phi_i, alpha)
  !=====================================================

    IMPLICIT NONE

    TYPE(element),              INTENT(IN)  :: ele
    REAL(KIND=8),               INTENT(IN)  :: Phi_tot
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: u
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: Phi_i
    REAL(KIND=8),               INTENT(OUT) :: alpha
    !--------------------------------------------------

    INTEGER :: Ns, k, j, Nu

    REAL(KIND=8) :: theta_e, xi, diff_u, k_s
    
    REAL(KIND=8)                   :: u_m, x_m, y_m
    REAL(KIND=8), DIMENSION(N_dim) :: a_m

    REAL(KIND=8), DIMENSION(SIZE(Phi_i)) :: Stab_i
    !--------------------------------------------------

    Ns = ele%N_points

    !-----------
    ! Mean state
    !----------------------------------------------
    u_m = SUM(u) / REAL(Ns)

    x_m = SUM(ele%Coords(1, :)) / REAL(Ns)
    y_m = SUM(ele%Coords(2, :)) / REAL(Ns)

    a_m = advection_speed(pb_type, u_m, x_m, y_m)

    !--------------------
    ! Diffusive parameter
    !-------------------------------------------------
    alpha = 0.d0

    DO k = 1, Ns
       alpha = 0.5d0 * MAX( alpha, ABS( DOT_PRODUCT(a_m, ele%rd_n(:, k)) ) )
    ENDDO

    alpha = 3.d0 * alpha

    ! LxF Scheme
    Phi_i = Phi_tot/REAL(Ns) + alpha*(u - u_m)

    ! Limiting
    Phi_i = limitation(Phi_i)

 stab_i = CIP_stabilization(ele, u) 
    
  END SUBROUTINE LLxFS_scheme
  !==========================

  !=============================================
  FUNCTION limitation(Phi_i__L) RESULT(Phi_i__H)
  !=============================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: Phi_i__L
 
    REAL(KIND=8), DIMENSION(SIZE(Phi_i__L)) :: Phi_i__H
    !---------------------------------------------------

    REAL(KIND=8), DIMENSION(SIZE(Phi_i__L)) :: beta

    REAL(KIND=8) :: Phi_tot

    REAL(KIND=8), PARAMETER :: eps = 1.0E-15
    !---------------------------------------------------

    Phi_tot = SUM(Phi_i__L)

    IF( ABS(Phi_tot) > eps ) THEN

       beta = Phi_i__L / Phi_tot 

       beta = MAX(0.d0, beta) + eps

       beta = beta / SUM(beta)

    ELSE

       beta = 0.d0

    ENDIF

    Phi_i__H = beta * Phi_tot
    
  END FUNCTION limitation
  !======================

  !==============================================
  FUNCTION CIP_stabilization(ele, u) RESULT(Stab)
  !==============================================

    IMPLICIT NONE

    TYPE(element),              INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: u

    REAL(KIND=8), DIMENSION(SIZE(u)) :: Stab
    !---------------------------------------------

    INTEGER :: if, iq, k
    !---------------------------------------------

    DO if = 1, ele%N_faces

write(*,*) 'if', if

       DO iq = 1, ele%faces(if)%f%N_quad

write(*,*) 'iq', iq

do k = 1, ele%N_points
          write(*,*) ele%faces(if)%f%p_Dphi_q(:, k, iq)
enddo
       ENDDO       

    ENDDO
    
write(*,*) '--------------'

  END FUNCTION CIP_stabilization
  !=============================
    
END MODULE LLxFS_method
