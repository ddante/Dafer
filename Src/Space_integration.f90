MODULE space_integration

  USE Element_Class

  USE geometry,       ONLY: N_dim, N_dofs, N_elements, elements

  USE init_problem,   ONLY: order, pb_type, visc, CFL, &
                            with_source

  USE models,         ONLY: advection_flux, strong_bc, &
                            source_term
  USE Num_scheme
  USE Gradient_Reconstruction

  USE Quadrature_rules, ONLY: oInt_n, Int_d

  IMPLICIT NONE

  !================================================
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: D_uu
  !================================================

  PRIVATE  
  PUBLIC :: compute_rhs

CONTAINS

  !====================================
  SUBROUTINE compute_rhs(uu, rhs, Dt_V)
  !====================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: uu
    REAL(KIND=8), DIMENSION(:), INTENT(OUT)   :: rhs
    REAL(KIND=8), DIMENSION(:), INTENT(OUT)   :: Dt_V
    !-------------------------------------------------

    TYPE(element) :: loc_ele, adh_ele

    INTEGER,      DIMENSION(:), ALLOCATABLE :: Nu
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: u1, Phi_i
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: u2

    REAL(KIND=8) :: Phi_tot, inv_dt

    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_1_q
    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_2_q

    REAL(KIND=8), DIMENSION(N_dim) :: p_Du_1_q, p_Du_2_q
    
    INTEGER :: Ns, je, i_f, k, iq, n_ele, istat
    !----------------------------------------------

    rhs = 0.d0

    Dt_v = 0.d0

    IF (.NOT. ALLOCATED(D_uu) ) THEN
       ALLOCATE( D_uu(N_dim, SIZE(uu)) )
    ENDIF
    
    D_uu = Compute_gradient(uu)

    DO je = 1, N_elements

       loc_ele = elements(je)%p

       Ns = loc_ele%N_points

       ALLOCATE( Nu(Ns), u1(Ns), Phi_i(Ns) )

       Nu = loc_ele%NU

       u1 = uu(Nu)

       !--------------------------
       ! Compute total fluctuation
       !------------------------------------
       Phi_tot = total_residual(loc_ele, u1)

       !--------------------------
       ! Update the face gradients
       !------------------------------------------------------
       DO i_f = 1, loc_ele%N_faces

          p_Dphi_1_q => loc_ele%faces(i_f)%f%p_Dphi_1_q
          p_Dphi_2_q => loc_ele%faces(i_f)%f%p_Dphi_2_q
          
          n_ele = loc_ele%faces(i_f)%f%c_ele

          IF( n_ele /= 0 ) THEN

             adh_ele = elements(n_ele)%p

             ALLOCATE( u2(adh_ele%N_points) )
             
             u2 = uu(adh_ele%NU)

             DO iq = 1, loc_ele%faces(i_f)%f%N_quad

                p_Du_1_q = 0.d0
                DO k = 1, loc_ele%N_points
                   p_Du_1_q = p_Du_1_q + p_Dphi_1_q(:, k,iq)*u1(k)
                ENDDO

                p_Du_2_q = 0.d0
                DO k = 1, adh_ele%N_points
                   p_Du_2_q = p_Du_2_q + p_Dphi_2_q(:, k,iq)*u2(k)
                ENDDO

                loc_ele%faces(i_f)%f%p_Du_1_q(:,iq) = p_Du_1_q
                loc_ele%faces(i_f)%f%p_Du_2_q(:,iq) = p_Du_2_q

             ENDDO

             DEALLOCATE( u2 )

             NULLIFY( p_Dphi_1_q, p_Dphi_2_q )

          ELSE

             loc_ele%faces(i_f)%f%p_Du_1_q = 0.d0
             loc_ele%faces(i_f)%f%p_Du_2_q = 0.d0

!!$             DO iq = 1, loc_ele%faces(i_f)%f%N_quad                !* 
!!$                p_Du_1_q = 0.d0                                    !*
!!$                DO k = 1, loc_ele%N_points                         !*
!!$                   p_Du_1_q = p_Du_1_q + p_Dphi_1_q(:, k,iq)*u1(k) !*
!!$                ENDDO
!!$             ENDDO

          ENDIF

       ENDDO
       !------------------------------------------------------

       !---------------------------
       ! Distribute the fluctuation
       !------------------------------------------------------------
       CALL distribute_residual(loc_ele, Phi_tot, u1, Phi_i, inv_dt)

       ! Gather the nadal residual
       rhs(Nu) = rhs(Nu) + Phi_i
       
       Dt_V(Nu) = Dt_V(Nu) + inv_dt
       
       DEALLOCATE( Nu, u1, Phi_i )

    ENDDO

    Dt_V = CFL/DT_V

    !---------------------
    ! Impose strong bc 
    ! RHS = 0 on the inlet
    !-------------------------------------
    CALL strong_bc(pb_type, visc, uu, rhs)

  END SUBROUTINE compute_rhs
  !=========================

  !==============================================
  FUNCTION total_residual(ele, u) RESULT(Phi_tot)
  !==============================================

    IMPLICIT NONE

    TYPE(element),              INTENT(IN) :: ele
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: u

    REAL(KIND=8) :: Phi_tot
    !---------------------------------------------

    REAL(KIND=8) :: x, y
    REAL(KIND=8) :: Phi_advec, Phi_diff, Source

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ff
    REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: SS

    INTEGER :: Ns, i, j
    !---------------------------------------------

    Phi_advec = 0.d0
    Phi_diff  = 0.d0
    Source    = 0.d0

    Ns = ele%N_points

    ALLOCATE( ff(N_dim, Ns), &
              SS(Ns)         )

    DO i = 1, Ns

       x = ele%Coords(1, i)
       y = ele%Coords(2, i)

       ff(:, i) = advection_flux(pb_type, u(i), x, y)
       SS(i)    = source_term(   pb_type, u(i), x, y)

    ENDDO

    Phi_advec = oInt_n(ele, ff)

    IF( with_source ) THEN
       Source = Int_d(ele, SS)
    ENDIF
    
    DEALLOCATE( ff, SS )

    Phi_tot = Phi_advec - Phi_diff - Source
                      !^^^       !^^^

  END FUNCTION total_residual
  !==========================

END MODULE space_integration
