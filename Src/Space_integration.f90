MODULE space_integration

  USE Element_Class

  USE geometry,       ONLY: N_dim, N_dofs, N_elements, elements

  USE init_problem,   ONLY: order, pb_type, visc, CFL, &
                            theta, theta_t, with_source

  USE models,         ONLY: advection_flux, strong_bc, &
                            source_term
  USE Num_scheme

  USE Quadrature_rules, ONLY: oInt_n

  IMPLICIT NONE
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

    TYPE(element) :: loc_ele

    INTEGER,      DIMENSION(:), ALLOCATABLE :: Nu
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: u, Phi_i

    REAL(KIND=8) :: Phi_tot, inv_dt

    INTEGER :: Ns, je, istat
    !----------------------------------------------

    rhs = 0.d0;  theta_t = 0.d0

    Dt_v = 0.d0

    DO je = 1, N_elements
                   
       loc_ele = elements(je)%p

       Ns = loc_ele%N_points

       ALLOCATE( Nu(Ns), u(Ns), Phi_i(Ns) )

       Nu = loc_ele%NU

       u = uu(Nu)

       !--------------------------
       ! Compute total fluctuation
       !-----------------------------------
       Phi_tot = total_residual(loc_ele, u)

       !----------------------------
       ! Distribute the fluctuation
       !-----------------------------------------------------------
       CALL distribute_residual(loc_ele, Phi_tot, u, Phi_i, inv_dt)

       ! Gather the nadal residual       
       rhs(Nu) = rhs(Nu) + Phi_i
       
       Dt_V(Nu) = Dt_V(Nu) + inv_dt
       
       DEALLOCATE( Nu, u, Phi_i )
      
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
    
    INTEGER :: Ns, i
    !---------------------------------------------

    Phi_advec = 0.d0
    Phi_diff  = 0.d0
    Source    = 0.d0

    Ns = ele%N_points

    ALLOCATE( ff(N_dim, Ns) )

    DO i = 1, Ns

       x = ele%Coords(1, i)
       y = ele%Coords(2, i)

       ff(:, i) = advection_flux(pb_type, u(i), x, y)

    ENDDO

    Phi_advec = oInt_n(ele, ff)

    DEALLOCATE( ff )

    Phi_tot = Phi_advec - Phi_diff - Source
                      !^^^

  END FUNCTION total_residual
  !==========================

END MODULE space_integration
