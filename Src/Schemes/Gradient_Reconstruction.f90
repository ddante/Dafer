MODULE Gradient_Reconstruction

  USE Element_class  
  USE geometry,       ONLY: N_dim, N_elements, elements, N_dofs, inv_A

  USE Quadrature_rules

  USE init_problem, ONLY: order, pb_type, visc
  USE Models,       ONLY: exact_grad

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
  INTEGER, PARAMETER :: RESERVED      = 0, &
                        GREEN_GAUSS   = 1, &
                        L2_PROJECTION = 2, &
                        LEAST_SQUARE  = 3
                        
  !===========================================
  integer :: type_reconstruction = 3

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

    SELECT CASE( type_reconstruction )

    CASE( GREEN_GAUSS )

       D_uu = GreenGauss_gradient(uu)

    CASE( L2_PROJECTION )

       D_uu = L2projection_gradient(uu)
              
    CASE( LEAST_SQUARE )

       D_uu = LSQ_gradient(uu)

    CASE DEFAULT

       WRITE(*,*) 'Unknown gradient reconstruction type'
       STOP
       
    END SELECT    

  END FUNCTION Compute_gradient
  !============================

  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!

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

!!$    ! Exact gradients
!!$    DO je = 1, N_elements
!!$       
!!$       loc_ele = elements(je)%p
!!$
!!$       DO k = 1, loc_ele%N_points
!!$
!!$          xy = loc_ele%Coords(:, k)
!!$          
!!$          D_uu(:, loc_ele%NU(k)) = exact_grad(pb_type, xy, visc)
!!$
!!$       ENDDO
!!$
!!$    ENDDO
    
    DEALLOCATE( Den )

  END FUNCTION GreenGauss_Gradient
  !===============================

  !=====================================
  FUNCTION LSQ_Gradient(uu) RESULT(D_uu)
  !=====================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: uu

    REAL(KIND=8), DIMENSION(N_dim, SIZE(uu)) :: D_uu
    !-------------------------------------------

    TYPE(element) :: ele

    TYPE :: vector
       REAL(KIND=8), DIMENSION(2) :: v
    END TYPE vector
    TYPE(vector), DIMENSION(:), ALLOCATABLE :: b

    INTEGER, DIMENSION(:), ALLOCATABLE :: Nu
    !-------------------------------------------

    REAL(KIND=8) :: x_i, y_i, x_k, y_k, u_i, u_k

    INTEGER :: je, i, k, jf, N1, N2
    !-------------------------------------------

    ALLOCATE( b(N_dofs) )
    DO i = 1, N_dofs
       b(i)%v = 0.d0
    ENDDO

    DO je = 1, N_elements

       ele = elements(je)%p

       ALLOCATE( NU(ele%N_points) )

       Nu = ele%NU
            
       DO jf = 1, ele%N_faces  

          N1 = Nu(ele%faces(jf)%f%l_nu(1))
          N2 = Nu(ele%faces(jf)%f%l_nu(2))         

          x_i = ele%Coords(1, ele%faces(jf)%f%l_nu(1))
          y_i = ele%Coords(2, ele%faces(jf)%f%l_nu(1))
          
          x_k = ele%Coords(1, ele%faces(jf)%f%l_nu(2))
          y_k = ele%Coords(2, ele%faces(jf)%f%l_nu(2))

          u_i = uu(N1)
          u_k = uu(N2)

          ! b --> i
          b(N1)%v(1) = b(N1)%v(1) + 2.d0*(x_k - x_i)*(u_k - u_i)
          b(N1)%v(2) = b(N1)%v(2) + 2.d0*(y_k - y_i)*(u_k - u_i)

          ! b --> k
          b(N2)%v(1) = b(N2)%v(1) + 2.d0*(x_i - x_k)*(u_i - u_k)
          b(N2)%v(2) = b(N2)%v(2) + 2.d0*(y_i - y_k)*(u_i - u_k)

          
       ENDDO

       DEALLOCATE(Nu)

    ENDDO

    DO i = 1, N_dofs      
       D_uu(:, i) = MATMUL(inv_A(i)%MM, b(i)%v)
    ENDDO

    DEALLOCATE( b )

  END FUNCTION LSQ_Gradient
  !========================  

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

    INTEGER :: je, i, j, N_s

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: M_ij

    INTEGER, DIMENSION(:), ALLOCATABLE :: rows, cols

    TYPE(element) :: ele
    
    PetscErrorCode :: ierr
    PetscInt :: nz
    !----------------------
    
    CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)

    ! approx non-0-elements for each row
    nz = 9*(order-1)

    CALL MatCreateSeqAIJ(PETSC_COMM_WORLD, N_dofs,N_dofs, &
                         nz, PETSC_NULL_INTEGER, MM_LL, ierr)

    CALL MatSetFromOptions(MM_LL, ierr)

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

       CALL MatSetValues(MM_LL, N_s,rows, N_s,cols, M_ij, ADD_VALUES, ierr)

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
