MODULE Gradient_Reconstruction

  USE Element_class  
  USE geometry,       ONLY: N_dim, N_elements, elements, &
                            N_dofs, inv_A, A_t, nn_NU

  USE Quadrature_rules

  USE init_problem, ONLY: order, pb_type, visc, grad_recovery
  USE Models,       ONLY: exact_grad, advection_speed

  IMPLICIT NONE

#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"

  !===================
  Mat  :: MM_LL
  KSP  :: ksp
  PC   :: pc
  Vec  :: rhs_x, rhs_y
  Vec  :: sol
  !===================
  LOGICAL :: is_factorized = .FALSE.

  !=========================================== 
  INTEGER, PARAMETER :: EXACT         = 0, &
                        GREEN_GAUSS   = 1, &
                        L2_PROJECTION = 2, &
                        LEAST_SQUARE  = 3, & 
                        SPR_ZZ        = 4
                        
  !===========================================
 
  PRIVATE
  PUBLIC :: Compute_gradient, DestroyPETSc

CONTAINS

  !=========================================
  FUNCTION Compute_gradient(uu) RESULT(D_uu)
  !=========================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: uu

    REAL(KIND=8), DIMENSION(N_dim, SIZE(uu)) :: D_uu 
    !-----------------------------------------------

    SELECT CASE( grad_recovery )

    CASE( EXACT )

       D_uu = Exact_gradient(uu)

    CASE( GREEN_GAUSS )

       D_uu = GreenGauss_gradient(uu)

    CASE( L2_PROJECTION )

       D_uu = L2projection_gradient(uu)
              
    CASE( LEAST_SQUARE )
       
       IF( Order == 2) THEN
          D_uu = LSQ_gradient_O2(uu)
       ELSE
          D_uu = LSQ_gradient_O3(uu)
       ENDIF

    CASE( SPR_ZZ )

       D_uu = SPRZZ_gradient(uu)

    CASE DEFAULT

       WRITE(*,*) 'Unknown gradient reconstruction type'
       STOP
       
    END SELECT    

  END FUNCTION Compute_gradient
  !============================

  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!

  !=======================================
  FUNCTION Exact_gradient(uu) RESULT(D_uu)
  !=======================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu

    REAL(KIND=8), DIMENSION(N_dim, SIZE(uu)) :: D_uu
    !------------------------------------------------

    TYPE(element) :: loc_ele

    REAL(KIND=8), DIMENSION(N_dim) :: xy

    INTEGER :: Nu, je, k
    !------------------------------------------------

    DO je = 1, N_elements

       loc_ele = elements(je)%p

       DO k = 1, loc_ele%N_points

          Nu = loc_ele%Nu(k)
          xy = loc_ele%Coords(:, k)

          D_uu(:, Nu) = exact_grad(pb_type, xy, visc)

       ENDDO

    ENDDO

  END FUNCTION Exact_gradient
  !==========================  

  !============================================
  FUNCTION GreenGauss_Gradient(uu) RESULT(D_uu)
  !============================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu

    REAL(KIND=8), DIMENSION(N_dim, SIZE(uu)) :: D_uu
    !------------------------------------------------

    TYPE(element) :: loc_ele
    
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: u
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Den

    REAL(KIND=8), DIMENSION(N_dim) :: D_u, xy

    INTEGER :: k, i, je, jf
    !------------------------------------------------

    ALLOCATE( Den(SIZE(uu)) )

    D_uu = 0.d0;  Den = 0.d0
    
    DO je = 1, N_elements

       loc_ele = elements(je)%p

       ALLOCATE( u(loc_ele%N_points) )

       u = uu( loc_ele%NU )

       DO k = 1, loc_ele%N_points

          DO i = 1, N_dim
             D_u(i) = SUM( u * loc_ele%D_phi_k(i, :, k) )
          ENDDO

          D_uu(:, loc_ele%NU(k)) = D_uu(:, loc_ele%NU(k)) + &
                                   D_u * loc_ele%volume

       ENDDO

       Den(loc_ele%NU) = Den(loc_ele%NU) + loc_ele%volume

       DEALLOCATE( u )

    ENDDO

    DO i = 1, N_dim
       D_uu(i, :) = D_uu(i, :)/Den
    ENDDO
    
    DEALLOCATE( Den )

  END FUNCTION GreenGauss_Gradient
  !===============================

  !========================================
  FUNCTION LSQ_Gradient_O2(uu) RESULT(D_uu)
  !========================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: uu

    REAL(KIND=8), DIMENSION(N_dim, SIZE(uu)) :: D_uu
    !-------------------------------------------

    TYPE(element) :: ele

    TYPE :: vector
       REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: v
    END TYPE vector
    TYPE(vector), DIMENSION(:), ALLOCATABLE :: b

    INTEGER, DIMENSION(:), ALLOCATABLE :: Nu
    !-------------------------------------------

    REAL(KIND=8) :: x_i, y_i, x_k, y_k, u_i, u_k, w_ik

    INTEGER :: je, i, j, k, Ni, Nk
    !-------------------------------------------

    ALLOCATE( b(N_dofs) )
    DO i = 1, N_dofs
       ALLOCATE( b(i)%v(2) )
       b(i)%v = 0.d0
    ENDDO

    DO je = 1, N_elements

       ele = elements(je)%p

       ALLOCATE( NU(ele%N_points) )

       Nu = ele%NU
            
       DO i = 1, ele%N_points

          Ni = Nu(i)

          DO k = 1, ele%N_points

             Nk = Nu(k)

             IF( k == i ) CYCLE

             DO j = 1, SIZE(nn_NU(Ni)%con)

                IF( Nk == nn_NU(Ni)%con(j) ) THEN

                   nn_NU(Ni)%con(j) = -nn_NU(Ni)%con(j)

                   x_i = ele%Coords(1, i)
                   y_i = ele%Coords(2, i)

                   x_k = ele%Coords(1, k)
                   y_k = ele%Coords(2, k)

                   w_ik = 1.d0 / ( (x_i - x_k)**2 + (y_i - y_k)**2 )

                   u_i = uu(Ni)
                   u_k = uu(Nk)

                   ! b --> i
                   b(Ni)%v(1) = b(Ni)%v(1) + 2.d0*w_ik*(x_k - x_i)*(u_k - u_i)
                   b(Ni)%v(2) = b(Ni)%v(2) + 2.d0*w_ik*(y_k - y_i)*(u_k - u_i)          

                ENDIF

             ENDDO

          ENDDO
      
       ENDDO

       DEALLOCATE(Nu)

    ENDDO

    DO i = 1, N_dofs      
       D_uu(:, i) = MATMUL(inv_A(i)%MM, b(i)%v)
       nn_NU(i)%con = -nn_NU(i)%con
    ENDDO

    DEALLOCATE( b )

  END FUNCTION LSQ_GRADIENT_O2  
  !===========================

  !========================================
  FUNCTION LSQ_Gradient_O3(uu) RESULT(D_uu)
  !========================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: uu

    REAL(KIND=8), DIMENSION(N_dim, SIZE(uu)) :: D_uu
    !-------------------------------------------

    TYPE(element) :: ele

    TYPE :: vector
       REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: v
    END TYPE vector
    TYPE(vector), DIMENSION(:), ALLOCATABLE :: b

    INTEGER, DIMENSION(:), ALLOCATABLE :: Nu
    !-------------------------------------------

    REAL(KIND=8), DIMENSION(5) :: x_v

    REAL(KIND=8) :: x_i, y_i, x_k, y_k, u_i, u_k, w_ik

    INTEGER :: je, i, j, k, Ni, Nk
    !-------------------------------------------

    ALLOCATE( b(N_dofs) )
    DO i = 1, N_dofs
       ALLOCATE( b(i)%v(5) )
       b(i)%v = 0.d0
    ENDDO

    DO je = 1, N_elements

       ele = elements(je)%p

       ALLOCATE( NU(ele%N_points) )

       Nu = ele%NU
            
       DO i = 1, ele%N_points

          Ni = Nu(i)

          DO k = 1, ele%N_points

             Nk = Nu(k)

             IF( k == i ) CYCLE

             DO j = 1, SIZE(nn_NU(Ni)%con)

                IF( Nk == nn_NU(Ni)%con(j) ) THEN

                   nn_NU(Ni)%con(j) = -nn_NU(Ni)%con(j)

                   x_i = ele%Coords(1, i)
                   y_i = ele%Coords(2, i)

                   x_k = ele%Coords(1, k)
                   y_k = ele%Coords(2, k)

                   w_ik = 1.d0 / ( (x_i - x_k)**2 + (y_i - y_k)**2 )

                   u_i = uu(Ni)
                   u_k = uu(Nk)

                   ! b --> i
                   b(Ni)%v(1) = b(Ni)%v(1) + 2.d0*w_ik*(u_k - u_i)*(x_k - x_i)
                   b(Ni)%v(2) = b(Ni)%v(2) + 2.d0*w_ik*(u_k - u_i)*(y_k - y_i)
                   b(Ni)%v(3) = b(Ni)%v(3) +      w_ik*(u_k - u_i)*(x_k - x_i)**2
                   b(Ni)%v(4) = b(Ni)%v(4) +      w_ik*(u_k - u_i)*(y_k - y_i)**2
                   b(Ni)%v(5) = b(Ni)%v(5) + 2.d0*w_ik*(u_k - u_i)*(x_k - x_i)*(y_k - y_i)

                ENDIF

             ENDDO

          ENDDO

       ENDDO
       
       DEALLOCATE(Nu)

    ENDDO
    
    DO i = 1, N_dofs

       x_v = MATMUL(inv_A(i)%MM, b(i)%v)

       D_uu(:, i) = x_v(1:2)

       nn_NU(i)%con = -nn_NU(i)%con

    ENDDO

    DEALLOCATE( b )

  END FUNCTION LSQ_GRADIENT_O3
  !===========================

  !=======================================
  FUNCTION SPRZZ_gradient(uu) RESULT(D_uu)
  !=======================================
  ! 
  ! WARNING: there is a problem in the case
  !          of a triangle with all the vertices
  !          wich belong to the boundaies
  !
    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu

    REAL(KIND=8), DIMENSION(N_dim, SIZE(uu)) :: D_uu
    !------------------------------------------------

    TYPE(element) :: ele

    TYPE :: vector
       REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: v
    END TYPE vector
    TYPE(vector), DIMENSION(:), ALLOCATABLE :: b, ss
    !------------------------------------------------

    REAL(KIND=8), DIMENSION(N_dim) :: GG

    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhs

    INTEGER, DIMENSION(:), ALLOCATABLE :: Nk

    REAL(KIND=8) :: x_0, y_0, x_k, y_k

    INTEGER :: je, i, j, jf, k, iq, n, Nu, Nu_i, Nu_b
    !------------------------------------------------

    D_uu = 0.d0

    ALLOCATE(  b(SIZE(A_t)), &
              ss(SIZE(A_t)) )

    DO i = 1, SIZE(b)

       n = SIZE(A_t(i)%MM, 2)
       ALLOCATE( b(i)%v(n, N_dim) )
       b(i)%v = 0.d0

       n = SIZE(inv_A(i)%MM, 1)
       ALLOCATE( ss(i)%v(n, N_dim) )
       ss(i)%v = 0.d0

    ENDDO

    ALLOCATE( Nk(SIZE(A_t)) )
    Nk = 0

    ! Construction of the RHS terms
    ! for each mesh vertex
    !--------------------------------
    DO je = 1, N_elements

       ele = elements(je)%p

       DO iq = 1, ele%N_rcv

          DO k = 1, N_dim
             GG(k) = SUM( uu(ele%Nu) * ele%D_phi_R(k, :, iq) )
          ENDDO
     
          DO i = 1, ele%N_verts

             Nu = ele%NU(i)

             Nk(Nu) = Nk(Nu) + 1

             b(Nu)%v(Nk(Nu),:) = GG

          ENDDO

       ENDDO

    ENDDO
    !--------------------------------

    DEALLOCATE(Nk)

    ! Polynomial coeff. for the recoverd
    ! gradient at each mesh vertex
    !--------------------------------
    DO i = 1, SIZE(b)
       
       ALLOCATE( rhs(SIZE(A_t(i)%MM, 1)) )

       DO k = 1, N_dim

          rhs = MATMUL( A_t(i)%MM, b(i)%v(:,k) )

          ss(i)%v(:, k) = MATMUL( inv_A(i)%MM, rhs )
          
          D_uu(k, i) = ss(i)%v(1, k)

       ENDDO

       DEALLOCATE( rhs )
       
    ENDDO
    !--------------------------------

    !/////////////////////////////////////////
    !// For the other DOFs the gradient     //
    !// is computed extrapolating the value //
    !// from the vertces                    //
    !/////////////////////////////////////////
    
    ! Dofs on the faces
    DO je = 1, N_elements
       
       ele = elements(je)%p

       DO j = 1, ele%N_faces
          
          ! boundary nodes
          IF( ele%faces(j)%f%c_ele == 0 ) THEN

             DO k = 1, ele%faces(j)%f%N_points
             
                Nu_b = ele%faces(j)%f%Nu(k)

                x_k = ele%faces(j)%f%Coords(1, k)
                y_k = ele%faces(j)%f%Coords(2, k)
                
                ! *** Find the domain vertex ***
                DO i = 1, ele%N_verts

                   n = 0

                   DO jf = 1, ele%N_faces

                      IF( ele%faces(jf)%f%c_ele == 0 .AND. &
                           (ele%Nu(i) == ele%faces(jf)%f%Nu(1) .OR. &
                            ele%Nu(i) == ele%faces(jf)%f%Nu(2)) ) THEN

                         n = n+1
                      ENDIF

                   ENDDO

                   IF(n == 0) THEN                        

                      Nu_i = ele%Nu(i)

                      x_0 = ele%Coords(1, i)
                      y_0 = ele%Coords(2, i)

                      EXIT

                   ENDIF

                ENDDO
                ! *** done ***
                
                SELECT CASE(Order) 

                CASE(2)

                   D_uu(1, Nu_b) = &
                        ss(Nu_i)%v(1, 1) + &
                        ss(Nu_i)%v(2, 1)*(x_k - x_0) + &
                        ss(Nu_i)%v(3, 1)*(y_k - y_0)

                   D_uu(2, Nu_b) = &
                        ss(Nu_i)%v(1, 2) + &
                        ss(Nu_i)%v(2, 2)*(x_k - x_0) + &
                        ss(Nu_i)%v(3, 2)*(y_k - y_0)

                CASE(3)

                   D_uu(1, Nu_b) = &
                        ss(Nu_i)%v(1, 1) + &
                        ss(Nu_i)%v(2, 1)*(x_k - x_0) + &
                        ss(Nu_i)%v(3, 1)*(y_k - y_0) + &
                        ss(Nu_i)%v(4, 1)*(x_k - x_0)**2 + &
                        ss(Nu_i)%v(5, 1)*(y_k - y_0)**2 + &
                        ss(Nu_i)%v(6, 1)*(x_k - x_0)*(y_k - y_0)

                   D_uu(2, Nu_b) = &
                        ss(Nu_i)%v(1, 2) + &
                        ss(Nu_i)%v(2, 2)*(x_k - x_0) + &
                        ss(Nu_i)%v(3, 2)*(y_k - y_0) + &
                        ss(Nu_i)%v(4, 2)*(x_k - x_0)**2 + &
                        ss(Nu_i)%v(5, 2)*(y_k - y_0)**2 + &
                        ss(Nu_i)%v(6, 2)*(x_k - x_0)*(y_k - y_0) 

                CASE DEFAULT

                   WRITE(*,*) 'ERROR, order nont supported'
                   STOP

                END SELECT
                   
             ENDDO
             
          ! Internal nodes
          ELSEIF( ele%faces(j)%f%c_ele /=0 .AND. Order > 2 ) THEN

             
             DO k = ele%faces(j)%f%N_Verts+1, ele%faces(j)%f%N_points

                Nu_b = ele%faces(j)%f%Nu(k)

                x_k = ele%faces(j)%f%Coords(1, k)
                y_k = ele%faces(j)%f%Coords(2, k)

                ! *** Find the domain vertex ***
                DO i = 1, ele%faces(j)%f%N_Verts

                   n = 0
                   
                   DO jf = 1, ele%N_faces

                      IF( ele%faces(jf)%f%c_ele == 0 .AND. &
                           (ele%faces(j)%f%Nu(i) == ele%faces(jf)%f%Nu(1) .OR. &
                            ele%faces(j)%f%Nu(i) == ele%faces(jf)%f%Nu(2)) ) THEN

                         n = n+1
                      ENDIF

                   ENDDO
                   
                   IF(n == 0) THEN
                      
                      Nu_i = ele%faces(j)%f%Nu(i)

                      x_0 = ele%faces(j)%f%Coords(1, i)
                      y_0 = ele%faces(j)%f%Coords(2, i)

                      EXIT

                   ENDIF

                ENDDO
                ! *** done ***

                D_uu(1, Nu_b) = &
                     ss(Nu_i)%v(1, 1) + &
                     ss(Nu_i)%v(2, 1)*(x_k - x_0) + &
                     ss(Nu_i)%v(3, 1)*(y_k - y_0) + &
                     ss(Nu_i)%v(4, 1)*(x_k - x_0)**2 + &
                     ss(Nu_i)%v(5, 1)*(y_k - y_0)**2 + &
                     ss(Nu_i)%v(6, 1)*(x_k - x_0)*(y_k - y_0)

                D_uu(2, Nu_b) = &
                     ss(Nu_i)%v(1, 2) + &
                     ss(Nu_i)%v(2, 2)*(x_k - x_0) + &
                     ss(Nu_i)%v(3, 2)*(y_k - y_0) + &
                     ss(Nu_i)%v(4, 2)*(x_k - x_0)**2 + &
                     ss(Nu_i)%v(5, 2)*(y_k - y_0)**2 + &
                     ss(Nu_i)%v(6, 2)*(x_k - x_0)*(y_k - y_0)

                ENDDO

             ENDIF

          ENDDO

          ! Internal dof
          IF(ele%N_verts+ele%N_faces < ele%N_points) THEN

             DO k = ele%N_verts+ele%N_faces+1, ele%N_points

                x_k = ele%Coords(1, k)
                y_k = ele%Coords(2, k)

                ! *** Find the domain vertex ***
                DO i = 1, ele%N_verts

                   n = 0

                   DO jf = 1, ele%N_faces

                      IF( ele%faces(jf)%f%c_ele == 0 .AND. &
                           (ele%Nu(i) == ele%faces(jf)%f%Nu(1) .OR. &
                            ele%Nu(i) == ele%faces(jf)%f%Nu(2)) ) THEN

                         n = n+1
                      ENDIF

                   ENDDO

                   IF(n == 0) THEN          

                      Nu_i = ele%Nu(i)

                      x_0 = ele%Coords(1, i)
                      y_0 = ele%Coords(2, i)

                      EXIT

                   ENDIF

                ENDDO
                ! *** done ***

                D_uu(1, ele%NU(k)) = &
                     ss(Nu_i)%v(1, 1) + &
                     ss(Nu_i)%v(2, 1)*(x_k - x_0) + &
                     ss(Nu_i)%v(3, 1)*(y_k - y_0) + &
                     ss(Nu_i)%v(4, 1)*(x_k - x_0)**2 + &
                     ss(Nu_i)%v(5, 1)*(y_k - y_0)**2 + &
                     ss(Nu_i)%v(6, 1)*(x_k - x_0)*(y_k - y_0)

                D_uu(2, ele%NU(k)) = &
                     ss(Nu_i)%v(1, 2) + &
                     ss(Nu_i)%v(2, 2)*(x_k - x_0) + &
                     ss(Nu_i)%v(3, 2)*(y_k - y_0) + &
                     ss(Nu_i)%v(4, 2)*(x_k - x_0)**2 + &
                     ss(Nu_i)%v(5, 2)*(y_k - y_0)**2 + &
                     ss(Nu_i)%v(6, 2)*(x_k - x_0)*(y_k - y_0)

             ENDDO

          ENDIF

       ENDDO

  END FUNCTION SPRZZ_gradient
  !==========================  

  !==============================================
  FUNCTION L2projection_gradient(uu) RESULT(D_uu)
  !==============================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu

    REAL(KIND=8), DIMENSION(N_dim, SIZE(uu)) :: D_uu
    !------------------------------------------------

    TYPE(element) :: ele
    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: b_xy

    INTEGER, ALLOCATABLE, DIMENSION(:) :: rows

    INTEGER :: k, i, j, je, N_s

    PetscScalar, DIMENSION(:), POINTER :: sol_v

    PetscReal :: zero
    PetscErrorCode :: ierr
    !------------------------------------------------

    zero = 0.d0

    IF( .NOT. is_factorized ) THEN

       CALL Assembly_MASSMatrix()
       
       is_factorized = .TRUE.

    ENDIF

    ! Construct the rhs for the two
    ! components of the gradients
    !------------------------------
    DO je = 1, N_elements

       ele = elements(je)%p

       N_s = ele%N_points

       ALLOCATE( b_xy(ele%N_dim, N_s) )
       ALLOCATE( rows(N_s) )

       DO i = 1, N_s

          b_xy(:, i) = int_d_Gu_i( ele, uu(ele%NU), i )
                
       ENDDO

       rows = ele%NU - 1

       CALL VecSetValues(rhs_x, N_s, rows, b_xy(1, :), ADD_VALUES, ierr)
       CALL VecSetValues(rhs_y, N_s, rows, b_xy(2, :), ADD_VALUES, ierr)

       DEALLOCATE( b_xy, rows )

    ENDDO

    CALL VecAssemblyBegin(rhs_x, ierr)
    CALL VecAssemblyBegin(rhs_y, ierr)

    CALL VecAssemblyEnd(rhs_x, ierr)
    CALL VecAssemblyEnd(rhs_y, ierr)

    ! Solve Grad_x
    !-----------------------------------
    CALL KSPSolve(ksp, rhs_x, sol, ierr)

    CALL VecGetArrayF90(sol, sol_v, ierr)
    D_uu(1, :) = sol_v
    CALL VecRestoreArrayF90(sol, sol_v, ierr)

    ! Solve Grad_y
    !-----------------------------------
    CALL KSPSolve(ksp, rhs_y, sol, ierr)

    CALL VecGetArrayF90(sol, sol_v, ierr)
    D_uu(2, :) = sol_v
    CALL VecRestoreArrayF90(sol, sol_v, ierr)

    ! Reset rhs
    !----------------------------
    CALL VecSet(rhs_x, zero, ierr)
    CALL VecSet(rhs_y, zero, ierr)    

  END FUNCTION L2projection_gradient  
  !=================================

  !===============================
  SUBROUTINE Assembly_MassMatrix()
  !===============================
  !
  ! Mass Matrix of the L2 projection.
  ! To be construced and factorized once for all
  !
    IMPLICIT NONE

    INTEGER :: je, i, j, k, N_s, n

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: M_ij

    INTEGER, DIMENSION(:), ALLOCATABLE :: rows, cols
    !-----------------------------------------------

    INTEGER, DIMENSION(:), ALLOCATABLE :: nnz

    TYPE :: con_Nu
       INTEGER, DIMENSION(:), ALLOCATABLE :: con
    END TYPE con_Nu
    TYPE(con_NU), DIMENSION(N_dofs) :: nn_Nu

    INTEGER, DIMENSION(:), POINTER :: Nu

    INTEGER, DIMENSION(:), ALLOCATABLE :: temp

    LOGICAL :: find

    TYPE(element) :: ele
    
    PetscErrorCode :: ierr
    !------------------------------------------------
    
    ! Count the #non-zero elements for each row
    ALLOCATE( nnz(N_dofs) )

    DO je = 1, N_elements
       
       N_s = elements(je)%p%N_points
       Nu => elements(je)%p%Nu

       DO i = 1, N_s

          IF( .NOT. ALLOCATED(nn_Nu(Nu(i))%con) ) THEN

             ALLOCATE(nn_Nu(Nu(i))%con(N_s))

             nn_Nu(Nu(i))%con = Nu

          ELSE

             DO k = 1, N_s

                n = SIZE(nn_Nu(Nu(i))%con)

                find = .FALSE.

                DO j = 1, n

                   IF(nn_Nu(Nu(i))%con(j) == Nu(k)) THEN
                      find = .TRUE.
                      EXIT
                   ENDIF

                ENDDO

                IF(.NOT. find) THEN

                   ALLOCATE( temp(SIZE(nn_Nu(Nu(i))%con)) )
                   temp = nn_Nu(Nu(i))%con

                   DEALLOCATE( nn_Nu(Nu(i))%con )
                     ALLOCATE( nn_Nu(Nu(i))%con(SIZE(temp)+1) )

                   nn_Nu(Nu(i))%con = (/ temp, Nu(k) /)
                      
                   DEALLOCATE(temp)                      

                ENDIF

             ENDDO ! k <= #N_points

          ENDIF ! Not allocated

       ENDDO ! i <= #N_points

    ENDDO ! je <= #elements


    DO i = 1, N_dofs
       nnz(i) = SIZE(nn_Nu(i)%con)
    ENDDO
    

    CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)

    CALL MatCreateSeqAIJ(PETSC_COMM_WORLD, N_dofs,N_dofs, &
                         PETSC_NULL_INTEGER, nnz, MM_LL, ierr)

    CALL MatSetFromOptions(MM_LL, ierr)

    DEALLOCATE( nnz )

    DO je = 1, N_elements

       ele = elements(je)%p

       N_s = ele%N_points

       ALLOCATE( M_ij(N_s, N_s) )
       ALLOCATE( rows(N_s), cols(N_s) )

       DO i = 1, ele%N_points

          DO j = 1, ele%N_points

             M_ij(i, j) = int_d_Mij(ele, i, j)

          ENDDO

       ENDDO

       rows = ele%NU - 1
       cols = ele%NU - 1
       
       DO i = 1,  ele%N_points
          CALL MatSetValues(MM_LL, 1,rows(i), N_s,cols, M_ij(i, :), ADD_VALUES, ierr)
       ENDDO

       DEALLOCATE( M_ij, rows, cols )

    ENDDO

    CALL MatAssemblyBegin(MM_LL, MAT_FINAL_ASSEMBLY, ierr)
    CALL MatAssemblyEnd(  MM_LL, MAT_FINAL_ASSEMBLY, ierr)

    CALL KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    CALL KSPSetOperators(ksp, MM_LL, MM_LL, SAME_PRECONDITIONER, ierr)
    
    CALL KSPGetPC(ksp, pc, ierr)
    CALL PCSetType(pc, PCLU, ierr)
    CALL PCFactorSetMatSolverPackage(pc, MATSOLVERMUMPS, ierr)

    CALL KSPSetFromOptions(ksp, ierr)

    ! Init RHS: G_x, G_y
    CALL VecCreate(PETSC_COMM_WORLD, rhs_x, ierr)
    CALL VecSetSizes(rhs_x, PETSC_DECIDE, N_dofs, ierr)
    CALL VecSetFromOptions(rhs_x, ierr)

    CALL VecDuplicate(rhs_x, rhs_y, ierr)
    CALL VecDuplicate(rhs_x, sol,   ierr)

  END SUBROUTINE Assembly_MassMatrix
  !=================================

  !========================
  SUBROUTINE DestroyPETSc()
  !========================

    IMPLICIT NONE

    PetscErrorCode :: ierr
    !----------------------

    IF( is_factorized ) THEN

       CALL MatDestroy(MM_LL, ierr)

       CALL VecDestroy(rhs_x, ierr)
       CALL VecDestroy(rhs_y, ierr)
       CALL VecDestroy(sol,   ierr)
       
       CALL KSPDestroy(ksp, ierr)

       CALL PetscFinalize(ierr)

    ENDIF

  END SUBROUTINE DestroyPETSc
  !==========================  

END MODULE Gradient_Reconstruction
