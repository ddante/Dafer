MODULE LLxFS_method

  USE Element_class
  USE geometry,       ONLY: N_dim, elements
  USE init_problem,   ONLY: pb_type
  USE models,         ONLY: advection_speed, source_term
  
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

    !-----------
    ! LxF Scheme
    !----------------------------------------
    Phi_i = Phi_tot/REAL(Ns) + alpha*(u - u_m)

    !---------
    ! Limiting
    !------------------------
    Phi_i = limitation(Phi_i)

    !--------------
    ! Stabilization
    !---------------------------------
    Stab_i = GLS_stabilization(ele, u)

    !Stab_i = CIP_stabilization(ele, u)

    Phi_i = Phi_i + Stab_i

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
  FUNCTION GLS_stabilization(ele, u) RESULT(Stab)
  !==============================================

    IMPLICIT NONE

    TYPE(element),              INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: u

    REAL(KIND=8), DIMENSION(SIZE(u)) :: Stab
    !---------------------------------------------

    REAL(KIND=8) :: x_i, y_i
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: a
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: SS
    
    REAL(KIND=8), DIMENSION(N_dim) :: D_u_q, a_q
    REAL(KIND=8)                   :: s_q

    REAL(KIND=8), DIMENSION(:),     POINTER :: w
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: D_phi_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: phi_q
    REAL(KIND=8) :: tau, Stab_D_u, Stab_D_phi
    
    INTEGER :: iq, i, N_s, N_v 
    !---------------------------------------------

    Stab = 0.d0

    N_s = ele%N_points
    N_v = ele%N_verts

    ALLOCATE( a(N_dim, N_s), SS(N_s) )
    
          w => ele%w_q
    D_phi_q => ele%D_phi_q
      phi_q => ele%phi_q

    SS = 0.d0

    DO i = 1, N_s
       
       x_i = ele%coords(1, i)
       y_i = ele%coords(2, i)

       a(:, i) = advection_speed(pb_type, u(i), x_i, y_i)
       SS(i)   = source_term(    pb_type, u(i), x_i, y_i)
       
    ENDDO

    !------------------
    ! Scaling parameter
    !---------------------------------------------
    tau = 0.d0
    
    DO i = 1, N_v

       tau = tau + &
             0.5d0 * ABS( DOT_PRODUCT(a(:, i), ele%rd_n(:, i)) ) / ele%volume

    ENDDO

    tau = 1.d0 / tau

    !-------------------
    ! Stabilization term
    !--------------------------------------------------
    DO iq = 1, ele%N_quad

       D_u_q = 0.d0
         a_q = 0.d0
         s_q = 0.d0

       DO i = 1, N_s
          D_u_q = D_u_q +  D_phi_q(:, i, iq) * u(i)
            a_q = a_q   +  phi_q(i, iq) * a(:, i)
            s_q = s_q   +  phi_q(i, iq) * SS(i)
       ENDDO

       Stab_D_u = DOT_PRODUCT(a_q, D_u_q) - s_q
       Stab_D_u = Stab_D_u * tau

       DO i = 1, N_s

         Stab_D_phi = DOT_PRODUCT(a_q, D_phi_q(:, i, iq))

         Stab(i) = Stab(i) + Stab_D_u * Stab_D_phi * w(iq)
               
      ENDDO

   ENDDO

   NULLIFY( w, D_phi_q, phi_q )
   DEALLOCATE( a, SS )

  END FUNCTION GLS_stabilization
  !=============================
  
  !==============================================
  FUNCTION CIP_stabilization(ele, u) RESULT(Stab)
  !==============================================

    IMPLICIT NONE

    TYPE(element),              INTENT(IN)  :: ele
    REAL(KIND=8), DIMENSION(:), INTENT(IN)  :: u

    REAL(KIND=8), DIMENSION(SIZE(u)) :: Stab
    !---------------------------------------------   
    
    INTEGER,                        POINTER :: N_quad
    INTEGER,                        POINTER :: N_points_f
    INTEGER,      DIMENSION(:),     POINTER :: loc
    REAL(KIND=8),                   POINTER :: h_f
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: n_1
    REAL(KIND=8), DIMENSION(:),     POINTER :: w
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: p
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_1_q
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_2_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: p_Du_1_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: p_Du_2_q
    INTEGER,      DIMENSION(:),     POINTER :: loc_con

    REAL(KIND=8), DIMENSION(N_dim) :: n_2
    REAL(KIND=8), DIMENSION(N_dim) :: a_q

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: a    
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: p_Dphi_l_2_q 
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: Jump_grad_phi
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: Jump_grad_u

    REAL(KIND=8) :: g_1, b_n, x_i, y_i

    INTEGER :: i, if, iq, k
    !---------------------------------------------

    Stab = 0.d0

    g_1 = 0.01
    b_n = 1.d0

    ALLOCATE ( a(N_dim, ele%N_points) )

    DO i = 1, ele%N_points

       x_i = ele%Coords(1, i)
       y_i = ele%Coords(2, i)

       a(:, i) = advection_speed(pb_type, u(i), x_i, y_i)

    ENDDO
    
    ! Construct the stabilization term for each node
    DO i = 1, SIZE(Stab)

       ! Split the boundary integral on the element
       ! as the sum of integral on each face
       DO if = 1, ele%N_faces

          N_quad     => ele%faces(if)%f%N_quad
          N_points_f => ele%faces(if)%f%N_points
          loc        => ele%faces(if)%f%l_nu
          w          => ele%faces(if)%f%w_q
          h_f        => ele%faces(if)%f%area
          n_1        => ele%faces(if)%f%n_q
          p          => ele%faces(if)%f%phi_q
          p_Dphi_1_q => ele%faces(if)%f%p_Dphi_1_q
          p_Dphi_2_q => ele%faces(if)%f%p_Dphi_2_q
          p_Du_1_q   => ele%faces(if)%f%p_Du_1_q
          p_Du_2_q   => ele%faces(if)%f%p_Du_2_q
          loc_con    => ele%faces(if)%f%loc_con

          ALLOCATE( Jump_grad_phi(N_quad), &
                    Jump_grad_u(N_quad)    )

          ALLOCATE( p_Dphi_l_2_q(N_dim, N_quad) )

          IF(loc_con(i) /= 0) THEN
             p_Dphi_l_2_q = p_Dphi_2_q(:, loc_con(i), :)             
          ELSE
             p_Dphi_l_2_q = 0.d0             
          ENDIF

b_n = 0.d0
DO iq = 1, N_quad
a_q = 0.0
DO k = 1, N_points_f
a_q = a_q + a(:, loc(k)) * p(k, iq)
ENDDO
b_n = MAX(b_n, ABS(DOT_PRODUCT(a_q, n_1(:, iq))) )
ENDDO

          ! Perform the integration on the single face
          DO iq = 1, N_quad

             n_2 = -n_1(:, iq)
          
             ! Jump of the gradient of the i-th shape function
             ! at the quadrature points
             Jump_grad_phi(iq) = DOT_PRODUCT(p_Dphi_1_q(:, i,iq), n_1(:, iq) ) + &
                                 DOT_PRODUCT(p_Dphi_l_2_q(:, iq), n_2 )

             ! Jump of the gradient of the sulution
             ! at the quadrature points
             Jump_grad_u(iq) = DOT_PRODUCT(p_Du_1_q(:, iq), n_1(:, iq) ) + &
                               DOT_PRODUCT(p_Du_2_q(:, iq), n_2)

             Stab(i) = Stab(i) + &
                       g_1 * b_n * h_f*h_f * &
                       w(iq) * Jump_grad_phi(iq) * Jump_grad_u(iq)

          ENDDO

          NULLIFY( N_quad, N_points_f, loc,    &
                   w, h_f, n_1, p,             &
                   p_Dphi_1_q, p_Dphi_2_q,     &
                   p_Du_1_q, p_Du_2_q, loc_con )

          DEALLOCATE( Jump_grad_phi, &
                      Jump_grad_u,   &
                      p_Dphi_l_2_q   )

       ENDDO

    ENDDO

    DEALLOCATE( a )

  END FUNCTION CIP_stabilization
  !=============================
    
END MODULE LLxFS_method
