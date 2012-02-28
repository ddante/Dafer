MODULE LW_method

  USE Element_class
  USE geometry,         ONLY: N_dim, elements
  USE init_problem,     ONLY: pb_type
  USE models,           ONLY: advection_speed, source_term
  USE init_problem,     ONLY: is_visc, visc
  USE Quadrature_rules, ONLY: int_d_G

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: LW_scheme
  
CONTAINS

  !========================================================
  SUBROUTINE LW_scheme(ele, Phi_tot, u, D_u, Phi_i, inv_dt)
  !========================================================

    IMPLICIT NONE

    TYPE(element),                INTENT(IN) :: ele
    REAL(KIND=8),                 INTENT(IN)  :: Phi_tot
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: u
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: D_u
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: Phi_i
    REAL(KIND=8),                 INTENT(OUT) :: inv_dt
    !--------------------------------------------------
    INTEGER :: Ns, j, iq, i

    REAL(KIND=8) :: u_m, x_m, y_m, s_i, s_q, &
                    N_m, mod_am, T_par, K_plus

    REAL(KIND=8), DIMENSION(N_dim) :: a_m, int_G_i, a_i, a_q

    REAL(KIND=8), DIMENSION(SIZE(Phi_i)) :: K_i, beta_i, &
                                            omega_i

    REAL(KIND=8) :: alpha, tau, Lr, Tr, Re_l, l_min, l_max
    !---------------------------------------------------

    REAL(KIND=8), DIMENSION(:),     POINTER :: w
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: D_phi_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: phi_q

    REAL(KIND=8), DIMENSION(N_dim) :: D_u_q, &
                                      Dr_u_q
    !--------------------------------------------------

    Ns = ele%N_points
    
    !-----------
    ! Mean state
    !-------------------------------------------------
    u_m = SUM(u) / REAL(Ns)

    x_m = SUM(ele%Coords(1, :)) / REAL(Ns)
    y_m = SUM(ele%Coords(2, :)) / REAL(Ns)

    a_m = advection_speed(pb_type, u_m, x_m, y_m)

    Lr = reference_length(a_m, visc)

    Tr = Lr / ( SQRT( SUM(a_m**2) ) + (visc/Lr) )

    Re_l = SQRT( SUM(a_m**2) ) * Lr / visc

    alpha = 1.0d0

    tau = alpha * Tr
    
    !-----------
    ! LW wcheme
    !-------------------------------------------------
    l_max = 0.0; l_min = HUGE(1.d0); K_plus = 0.0

    DO j = 1, ele%N_verts !Ns

       K_i(j) = 0.5d0 * DOT_PRODUCT(a_m, ele%rd_n(:, j))

       l_max = MAX( l_max, SQRT( SUM(ele%rd_n(:, j)**2)) )

       l_min = MIN( l_min, SQRT( SUM(ele%rd_n(:, j)**2)) )

      !K_plus = K_plus + MAX( K_i(j), 0.d0 )
       K_plus = K_plus + ABS( K_i(j) )
       
    ENDDO

    K_i    = K_i    / ele%volume
    K_plus = K_plus / ele%volume

    mod_am = SQRT( SUM( a_m**2 ) )

!   T_par = 0.25d0 * ele%volume / MAXVAL(ABS(K_i))
!   T_par = 0.5d0 * ele%volume / K_plus
    T_par = 2.d0 / ( K_plus * ele%N_verts )

    omega_i = 1.d0/REAL(Ns, 8)

!!$    DO i = 1, Ns
!!$
!!$       int_G_i = int_d_G(ele, i)
!!$
!!$       beta_i(i) = omega_i(i) + &
!!$                   omega_i(i) * DOT_PRODUCT(a_m, int_G_i) / K_plus
!!$
!!$    !   beta_i(i) = omega_i(i) + &
!!$    !              (T_par / ele%volume) * DOT_PRODUCT(a_m, int_G_i)
!!$                 
!!$    ENDDO
!!$    
!!$    phi_i = beta_i * Phi_tot

    phi_i = omega_i * Phi_tot

          w => ele%w_q
    D_phi_q => ele%D_phi_q
      phi_q => ele%phi_q


    DO iq = 1, ele%N_quad

       Dr_u_q = 0.d0 ! Reconstructed gradient      
        D_u_q = 0.d0 ! Continuous gradient
          a_q = 0.d0 
          s_q = 0.d0

       DO i = 1, Ns

          Dr_u_q = Dr_u_q + phi_q(i, iq) * D_u(:, i)
           D_u_q = D_u_q  + D_phi_q(:, i, iq) * u(i)

          a_i = advection_speed( pb_type, u(i),    &
                                 ele%Coords(1, i), &
                                 ele%Coords(2, i)  )

          s_i = source_term( pb_type, u(i), visc, &
                             ele%Coords(1, i),    &
                             ele%Coords(2, i)     )

          a_q = a_q + phi_q(i, iq) * a_i

          s_q = s_q + phi_q(i, iq) * s_i

       ENDDO

       DO i = 1, Ns

          !----------------
          ! Advective part
          !--------------------------------------------------------
          phi_i(i) = phi_i(i) + &
                     T_par * ( DOT_PRODUCT(a_q, D_u_q) - s_q ) * &
                     DOT_PRODUCT(a_q, D_phi_q(:, i, iq)) * w(iq)

          !--------------
          ! Viscous part
          !--------------------------------------------------------
          IF(is_visc) THEN
             
             phi_i(i) = phi_i(i) + &
                        visc * DOT_PRODUCT(D_phi_q(:, i, iq), (D_u_q - Dr_u_q)) * w(iq) 

          ENDIF

       ENDDO

    ENDDO


!!$    !--------------
!!$    ! Viscous part
!!$    !--------------------------------------------------------
!!$    IF(is_visc) THEN
!!$
!!$              w => ele%w_q
!!$        D_phi_q => ele%D_phi_q
!!$          phi_q => ele%phi_q
!!$
!!$       DO iq = 1, ele%N_quad
!!$
!!$          Dr_u_q = 0.d0 ! Reconstructed gradient
!!$           D_u_q = 0.d0 ! Continuous gradient
!!$
!!$          DO i = 1, Ns
!!$
!!$             Dr_u_q = Dr_u_q + phi_q(i, iq) * D_u(:, i)
!!$              D_u_q = D_u_q  + D_phi_q(:, i, iq) * u(i)
!!$
!!$          ENDDO
!!$
!!$          DO i = 1, Ns
!!$
!!$             phi_i(i) = phi_i(i) + &
!!$                        visc * DOT_PRODUCT(D_phi_q(:, i, iq), (D_u_q - Dr_u_q)) * w(iq) 
!!$
!!$          ENDDO
!!$
!!$       ENDDO
!!$
!!$    ENDIF

    ! Time step
    !--------------
    inv_dt = 0.0

    DO j = 1, Ns

       inv_dt = MAX( inv_dt, &
                     DOT_PRODUCT(a_m, ele%rd_n(:, j)) * (1.d0 + 1.d0/Re_l) )

    ENDDO

    inv_dt = 2.d0 * inv_dt

    NULLIFY( w, D_phi_q, phi_q )

  END SUBROUTINE LW_scheme
  !=======================

  !==========================================
  FUNCTION reference_length(a, nu) RESULT(Lr)
  !==========================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: a
    REAL(KIND=8),               INTENT(IN) :: nu

    REAL(KIND=8) :: Lr
    !-------------------------------------------

    REAL(KIND=8), PARAMETER :: Pi = DACOS(-1.d0)
    REAL(KIND=8) :: Re_pi
    !-------------------------------------------

    IF( is_visc ) THEN
       Re_pi = ( SQRT( SUM(a**2) ) / Pi ) / (nu)
    ELSE
       Re_pi = HUGE(1.d0)
    ENDIF

    Lr = (0.5d0/Pi) * &
         (  ( Re_pi/(SQRT(1.d0 + Re_pi*Re_pi) + 1.d0) ) + &
         SQRT( 1.d0 + 2.d0/(SQRT(1.d0 + Re_pi*Re_pi) + 1.d0) ) )

  END FUNCTION reference_length
  !============================  

  !==============================================
  FUNCTION CIP_stabilization2(ele, u) RESULT(Stab)
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
    REAL(KIND=8), DIMENSION(:),     POINTER :: w
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: n_1
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: p
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_1_q
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_2_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: p_Du_1_q
    REAL(KIND=8), DIMENSION(:,:),   POINTER :: p_Du_2_q
    INTEGER,      DIMENSION(:),     POINTER :: loc_con

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: a_q
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: a
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: p_Dphi_l_2_q 

    REAL(KIND=8) :: Jump_D_phi, Jump_D_u
    REAL(KIND=8) :: b_n, b_t, g_1, x_i, y_i, mod_a_q

    INTEGER :: i, if, iq, k
    !---------------------------------------------

    Stab = 0.d0

    g_1 = 0.05

    ALLOCATE ( a(N_dim, ele%N_points) )
    DO i = 1, ele%N_points

       x_i = ele%Coords(1, i)
       y_i = ele%Coords(2, i)

       a(:, i) = advection_speed(pb_type, u(i), x_i, y_i)

    ENDDO
   
    DO i = 1, SIZE(Stab)

       ! Split the boundary integral on the element
       ! as the sum of integral on each face
       DO if = 1, ele%N_faces

          N_quad     => ele%faces(if)%f%N_quad      ! ok
          N_points_f => ele%faces(if)%f%N_points    ! ok
          loc        => ele%faces(if)%f%l_nu        ! ok
          n_1        => ele%faces(if)%f%n_q         ! ok
          w          => ele%faces(if)%f%w_q         ! ok
          h_f        => ele%faces(if)%f%area        ! ok
          p          => ele%faces(if)%f%phi_q       ! ok
          p_Dphi_1_q => ele%faces(if)%f%p_Dphi_1_q  ! ok
          p_Dphi_2_q => ele%faces(if)%f%p_Dphi_2_q  ! ok
          p_Du_1_q   => ele%faces(if)%f%p_Du_1_q    ! ok
          p_Du_2_q   => ele%faces(if)%f%p_Du_2_q    ! ok
          loc_con    => ele%faces(if)%f%loc_con     ! ok

          ALLOCATE( p_Dphi_l_2_q(N_dim, N_quad), & 
                             a_q(N_dim, N_quad) )

          IF(loc_con(i) /= 0) THEN
             p_Dphi_l_2_q = p_Dphi_2_q(:, loc_con(i), :)
          ELSE
             p_Dphi_l_2_q = 0.d0
          ENDIF

          DO iq = 1, N_quad

             a_q(:, iq) = 0.d0
             DO k = 1, N_points_f
                a_q(:, iq) = a_q(:, iq) + a(:, loc(k)) * p(k, iq) 
             ENDDO              
             
          ENDDO
          

          DO iq = 1, N_quad

             !---------------
             ! Gradient Jump
             !-----------------------------------------------------------
             Jump_D_phi = DOT_PRODUCT( p_Dphi_1_q(:, i, iq), n_1(:, iq) ) + &
                          DOT_PRODUCT( p_Dphi_l_2_q(:, iq), -n_1(:, iq) )
             
             Jump_D_u   = DOT_PRODUCT( p_Du_1_q(:, iq),  n_1(:, iq) ) + &
                          DOT_PRODUCT( p_Du_2_q(:, iq), -n_1(:, iq) )

             b_n = ABS( DOT_PRODUCT(a_q(:, iq), n_1(:,iq)) )
             
             !--------------------
             ! Stabilization term
             !-----------------------------------------------------------
             Stab(i) = Stab(i) + g_1 * b_n * &
                       w(iq) * (h_f**2) * (Jump_D_phi * Jump_D_u)

          ENDDO

          !--------------------
          ! Destruc face stuff
          !----------------------------------------
          DEALLOCATE( p_Dphi_l_2_q, a_q )

          NULLIFY( N_quad, N_points_f, loc, w, n_1, &
                   h_f, p, p_Dphi_1_q, p_Dphi_2_q,  &
                   p_Du_1_q, p_Du_2_q, loc_con      )


       ENDDO

    ENDDO

    DEALLOCATE ( a )

  END FUNCTION CIP_stabilization2
  !=============================

END MODULE LW_method
