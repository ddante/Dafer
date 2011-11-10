MODULE space_integration

  USE Element_Class

  USE geometry,       ONLY: N_dim, N_dofs, N_elements, elements

  USE init_problem,   ONLY: order, pb_type, visc, CFL, &
                            theta, theta_t, with_source

  USE models,         ONLY: flux_function, strong_bc, &
                            source_term
  !USE Num_scheme

  !USE Quadrature

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
       !------------------------------------
       Phi_tot = total_residual(loc_ele, u)

       DEALLOCATE( Nu, u, Phi_i )
      
    ENDDO

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

    REAL(KIND=8), DIMENSION(N_dim) :: ff
    
    INTEGER :: Ns, i, j, k, kk
    !---------------------------------------------

    Phi_advec = 0.d0
    Phi_diff  = 0.d0
    Source    = 0.d0


    Ns = ele%N_points

! --- TEST --- OK
write(*,*) ele%Type, '***'
do k = 1, ele%N_verts
j = MOD(k, MIN(k+1, ele%N_verts)) + 1
kk = k + ele%N_verts
write(*,*) ele%NU(k), ele%NU(j), ele%NU(kk)
enddo
do j = 1, ele%N_faces
write(*,*) ele%faces(j)%f%NU(:), 'f'
enddo
write(*,*)
! --- TEST ---

    DO i = 1, Ns

       x = ele%Coords(1, i)
       y = ele%Coords(2, i)

       ff = flux_function(pb_type, u(i), x, y)

       !Phi_advec = Int_AllEdges(ele, ff, Order)

!!$       Phi_advec = Phi_advec + &
!!$        DOT_PRODUCT( ff, ele%normale(i) ) * ele%weights(i)
       
    ENDDO

    ! IF visc, compute Phi_diff
!!$
!!$    Source = compute_source_int(ele, u)

    Phi_tot = Phi_advec - Phi_diff - Source
                      !^^^
         
  END FUNCTION total_residual
  !==========================

END MODULE space_integration
