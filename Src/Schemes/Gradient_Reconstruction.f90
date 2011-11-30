MODULE Gradient_Reconstruction

  USE Element_class  
  USE geometry,       ONLY: N_dim, N_elements, elements

  IMPLICIT NONE

  !=========================================== 
  INTEGER, PARAMETER :: RESERVED      = 0, &
                        LEAST_SQUARE  = 1
  !===========================================
  integer :: type_reconstruction = 1

  PRIVATE
  PUBLIC :: Compute_gradient

CONTAINS

  !=========================================
  FUNCTION Compute_gradient(uu) RESULT(D_uu)
  !=========================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: uu

    REAL(KIND=8), DIMENSION(N_dim, SIZE(uu)) :: D_uu
    !-----------------------------------------------

    SELECT CASE( type_reconstruction )

    CASE( LEAST_SQUARE )

       D_uu = Least_Square_gradient(uu)

    CASE DEFAULT

       WRITE(*,*) 'Unknown gradient reconstruction type'
       STOP
       
    END SELECT    

  END FUNCTION Compute_gradient
  !============================



  !==============================================
  FUNCTION Least_Square_Gradient(uu) RESULT(D_uu)
  !==============================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: uu

    REAL(KIND=8), DIMENSION(N_dim, SIZE(uu)) :: D_uu
    !------------------------------------------------

    TYPE(element) :: loc_ele
    
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: u
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Den

    REAL(KIND=8), DIMENSION(N_dim) :: D_u

    INTEGER :: k, i, je
    !------------------------------------------------

    D_uu = 0.d0

    ALLOCATE( Den(SIZE(uu)) )
    
    DO je = 1, N_elements

       loc_ele = elements(je)%p

       ALLOCATE( u(loc_ele%N_points) )

       u = uu( loc_ele%NU )

       DO k = 1, loc_ele%N_points

          DO i = 1, N_dim
             D_u(i) = SUM( u * loc_ele%D_phi_k(i, :, k) )
          ENDDO

          D_uu(:, loc_ele%NU(k)) = D_uu(:, loc_ele%NU(k)) + &
                                   D_u*loc_ele%volume

       ENDDO

       Den(loc_ele%NU) = Den(loc_ele%NU) + loc_ele%volume

       DEALLOCATE( u )

    ENDDO

    DO i = 1, N_dim
       D_uu(i, :) = D_uu(i, :)/Den
    ENDDO
    
    DEALLOCATE( Den )

  END FUNCTION Least_Square_Gradient
  !=================================

  
!!$  !=================================
!!$  SUBROUTINE Face_Average_Gradient()
!!$  !=================================
!!$
!!$    IMPLICIT NONE
!!$
!!$    TYPE(element) :: ele    
!!$    
!!$    REAL(KIND=8), DIMENSION(N_dim) :: Mean_grad_u
!!$
!!$    REAL(KIND=8), DIMENSION(:,:), POINTER :: p_Du_1_q
!!$    REAL(KIND=8), DIMENSION(:,:), POINTER :: p_Du_2_q
!!$
!!$    INTEGER :: je, if, iq
!!$    !-------------------------------------------------
!!$
!!$    DO je = 1, N_elements
!!$
!!$       ele = elements(je)%p
!!$
!!$       DO if = 1, ele%N_faces
!!$
!!$          N_quad   => ele%faces(if)%f%N_quad
!!$          p_Du_1_q => ele%faces(if)%f%p_Du_1_q
!!$          p_Du_2_q => ele%faces(if)%f%p_Du_2_q
!!$          
!!$          DO iq = 1, N_quad
!!$
!!$             Mean_grad_u = 0.5d0 * &
!!$                          ( p_Du_1_q(:, iq) + p_Du_2_q(:, iq) )
!!$
!!$             !ele%faces(i_f)%f%G_u(:, iq) = Mean_grad_u
!!$
!!$          ENDDO
!!$          
!!$
!!$       ENDDO
!!$
!!$       NULLIFY( N_quad, p_Du_1_q, p_Du_2_q ) 
!!$
!!$    ENDDO
!!$    
!!$  END SUBROUTINE Face_Average_Gradient
!!$  !===================================  

END MODULE Gradient_Reconstruction
