MODULE Quadrangle_Class

  USE Element_Class
  USE Segment_Class
  USE Lin_Algebra

  IMPLICIT NONE

  !=========================================
  TYPE, PUBLIC, EXTENDS(element) :: quadrangle

     ! < >

   CONTAINS

     PROCEDURE, PUBLIC :: initialize => initialize_sub

  END TYPE quadrangle
  !=========================================

  PRIVATE :: initialize_sub

  PRIVATE :: init_faces
  PRIVATE :: init_faces_QUA_Q1
  PRIVATE :: init_faces_QUA_Q2

  PRIVATE :: volume_quadrature
  PRIVATE :: init_quadrature_QUA_Q1
  PRIVATE :: init_quadrature_QUA_Q2

  PRIVATE :: basis_function
  PRIVATE :: basis_function_QUA_Q1
  PRIVATE :: basis_function_QUA_Q2

  PRIVATE :: gradient

  PRIVATE :: gradient_ref
  PRIVATE :: gradient_ref_QUA_Q1
  PRIVATE :: gradient_ref_QUA_Q2

  PRIVATE :: DOFs_ref
  PRIVATE :: DOFs_ref_QUA_Q1
  PRIVATE :: DOFs_ref_QUA_Q2

  PRIVATE :: gradient_trace

  PRIVATE :: face_trace
  PRIVATE :: face_trace_QUA_Q1
  PRIVATE :: face_trace_QUA_Q2

  PRIVATE :: rd_normal
  PRIVATE :: rd_normal_QUA_Q1
  PRIVATE :: rd_normal_QUA_Q2

  PRIVATE :: Compute_Jacobian

  PRIVATE :: recovery_procedure
  PRIVATE :: init_Gradiet_Recovery_QUA_Q1
  PRIVATE :: init_Gradiet_Recovery_QUA_Q2

CONTAINS

  !===============================================================
  SUBROUTINE initialize_sub(e, mode, Nodes, Coords, Nu_seg, n_ele)
  !===============================================================

    IMPLICIT NONE

    CLASS(quadrangle)                        :: e
    CHARACTER(*),                 INTENT(IN) :: mode
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords    
    INTEGER,      DIMENSION(:),   INTENT(IN) :: NU_seg
    INTEGER,      DIMENSION(:),   INTENT(IN) :: n_ele
    !-------------------------------------------------

    INTEGER :: i, id
    !-------------------------------------------------

    !---------------------------
    ! Attach data to the element
    !---------------------------------------------------
    ALLOCATE( e%Coords(SIZE(Coords,1), SIZE(Coords,2)) )
    ALLOCATE( e%NU(SIZE(Nodes)) )
     
    DO id = 1, SIZE(Coords,1)
       e%Coords(id, :) = Coords(id, :)
    ENDDO

    e%NU = Nodes

    SELECT CASE( SIZE(Nodes) )

    CASE(4)

       !---------
       ! ELELMENT
       !----------------------------------
       e%Type     = QUA_Q1 ! Element type
       e%N_verts  = 4      ! # Vertices
       e%N_points = 4      ! # DoFs

    CASE(9)

       !---------
       ! ELELMENT
       !----------------------------------       
       e%Type     = QUA_Q2 ! Element type
       e%N_verts  = 4      ! # Vertices
       e%N_points = 9      ! # DoFs

     CASE DEFAULT

        WRITE(*,*) 'ERROR: Unsupported Quadrangle type'
        STOP

     END SELECT

     IF( mode == "element") THEN

        CALL init_faces(e, Nodes, Coords, NU_seg, n_ele)

        ! Store the informations at the quadrature points
        CALL volume_quadrature(e)

        !--------------------------
        ! Normals for the RD scheme
        !--------------------------------------
        ALLOCATE( e%rd_n(e%N_dim, e%N_points) )
        
        DO i = 1, e%N_points
           e%rd_n(:, i) = rd_normal(e, i)
        ENDDO

        CALL nodal_gradients(e)

        ! Store the information at the recovery points
        CALL recovery_procedure(e)

     ENDIF
   
  END SUBROUTINE initialize_sub
  !============================

  !=====================================================
  SUBROUTINE init_faces(e, Nodes, Coords, NU_seg, n_ele)
  !=====================================================

    IMPLICIT NONE
    
    CLASS(quadrangle)                        :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    INTEGER,      DIMENSION(:),   INTENT(IN) :: NU_seg
    INTEGER,      DIMENSION(:),   INTENT(IN) :: n_ele
    !-------------------------------------------------

    SELECT CASE(e%Type)
              
    CASE(QUA_Q1)

       CALL init_faces_QUA_Q1(e, Nodes, Coords, NU_seg, n_ele)
       
    CASE(QUA_Q2)

       CALL init_faces_QUA_Q2(e, Nodes, Coords, NU_seg, n_ele)
              
    CASE DEFAULT

       WRITE(*,*) 'Unknown Quadrangle type for face initialization'
       STOP

    END SELECT 

  END SUBROUTINE init_faces
  !=========================

  !==============================
  SUBROUTINE volume_quadrature(e)
  !==============================
  !
  ! Store the following  quantities at the quadrature points
  !    - value of the basis function
  !    - normal versor
  !    - weight of the quadrature formula (multipplied by 
  !      the jacobian of the transformation)
  !    - value of the trace of the gradient of basis functions
  !    - physical coordiantes of the quadrature point
  !
  ! Compute the lenght of the segment
  !
    IMPLICIT NONE

    CLASS(quadrangle) :: e
    !-----------------------------------------------

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: JJ

    INTEGER :: iq, k
    !-----------------------------------------------

    SELECT CASE(e%Type)

    CASE(QUA_Q1)

       CALL init_quadrature_QUA_Q1(e)
       
    CASE(QUA_Q2)

       CALL init_quadrature_QUA_Q2(e)
           
    CASE DEFAULT

       WRITE(*,*) 'Unknown Quadrangle type for quadrature'
       WRITE(*,*) 'STOP'

    END SELECT

    !-------------------------------------    
    ! Attach data to the quadrature points
    !---------------------------------------------------
    ALLOCATE( JJ(e%N_dim, e%N_dim) )

    DO iq = 1, e%N_quad

       DO k = 1, e%N_points

          e%phi_q(k, iq) = basis_function( e, k, e%x_q(iq, :) )

          e%xx_q(:, iq) = e%xx_q(:, iq) + &               
                         basis_function( e, k, e%x_q(iq, :) ) * e%Coords(:, k)

          e%D_phi_q(:, k, iq) = gradient( e, k, e%x_q(iq, :) )

       ENDDO

       ! warning: not optimized, jacobian computed twice       
       JJ = Compute_Jacobian( e, e%x_q(iq, :) )

       e%w_q(iq) = e%w_q(iq) * determinant(JJ)
      
    ENDDO

    DEALLOCATE( JJ )

    ! Area of the element
    e%volume = SUM( e%w_q )    

  END SUBROUTINE volume_quadrature
  !===============================

  !============================
  SUBROUTINE nodal_gradients(e)
  !============================
  !
  ! Compute the gradient of the basis functions
  ! at the DOFs of the element
  !
    IMPLICIT NONE

    CLASS(quadrangle) :: e    
    !-----------------------------------------------

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: x_dof

    INTEGER :: i, k
    !-----------------------------------------------

    ALLOCATE( x_dof(e%N_points, 2))
    
    x_dof = DOFs_ref(e)
    
    ALLOCATE( e%D_phi_k(e%N_dim, e%N_points, e%N_points) )

    DO i = 1, e%N_points

       DO k = 1, e%N_points

          e%D_phi_k(:, i, k) = gradient( e, i, x_dof(k, :) )

       ENDDO

    ENDDO
    
    DEALLOCATE(x_dof)

  END SUBROUTINE nodal_gradients
  !=============================

  !===============================
  SUBROUTINE recovery_procedure(e)
  !===============================

    IMPLICIT NONE

    CLASS(quadrangle) :: e    
    !-----------------------------------------------

    INTEGER :: iq, k
    !-----------------------------------------------

    SELECT CASE(e%Type)

    CASE(QUA_Q1)

       CALL init_Gradiet_Recovery_QUA_Q1(e)
       
    CASE(QUA_Q2)

       CALL init_Gradiet_Recovery_QUA_Q2(e)
       
    CASE DEFAULT

       WRITE(*,*) 'Unknown Quadrangle type for recovery'
       WRITE(*,*) 'STOP'

    END SELECT

    !-------------------------------------    
    ! Attach data to the recovey points
    !---------------------------------------------------
    DO iq = 1, e%N_rcv

       DO k = 1, e%N_points
         
          e%xx_R(:, iq) = e%xx_R(:, iq) + &               
                          basis_function( e, k, e%x_R(iq, :) ) * e%Coords(:, k)

          e%D_phi_R(:, k, iq) = gradient( e, k, e%x_R(iq, :) )

       ENDDO
     
    ENDDO

  END SUBROUTINE recovery_procedure
  !================================

!*******************************************************************************
!*******************************************************************************
!                         COMPLEMENTARY FUNCTIONS                              !
!*******************************************************************************
!*******************************************************************************

  !==============================================
  FUNCTION basis_function(e, i, xi) RESULT(psi_i)
  !==============================================

    IMPLICIT NONE

    CLASS(quadrangle)                      :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8) :: psi_i
    !-------------------------------------------
      
    SELECT CASE(e%Type)
       
    CASE(QUA_Q1)

       psi_i = basis_function_QUA_Q1(e, i, xi)
       
    CASE(QUA_Q2)

       psi_i = basis_function_QUA_Q2(e, i, xi)
       
    CASE DEFAULT

       WRITE(*,*) 'Unknown Quadrangle type for basis functions'
       WRITE(*,*) 'STOP'

    END SELECT

  END FUNCTION basis_function
  !========================== 

  !==============================================
  FUNCTION gradient_ref(e, i, xi) RESULT(D_psi_i)
  !==============================================

    IMPLICIT NONE

    CLASS(quadrangle)                       :: e
    INTEGER,                     INTENT(IN) :: i    
    REAL(KIND=8),  DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(2) :: D_psi_i
    !-------------------------------------------
      
    SELECT CASE(e%Type)
       
    CASE(QUA_Q1)

       D_psi_i = gradient_ref_QUA_Q1(e, i, xi)
       
    CASE(QUA_Q2)

       D_psi_i = gradient_ref_QUA_Q2(e, i, xi)
     
    CASE DEFAULT

       WRITE(*,*) 'Unknown Quadrangle type for gradients'
       WRITE(*,*) 'STOP'

    END SELECT    

  END FUNCTION gradient_ref
  !========================

  !==========================================
  FUNCTION gradient(e, i, xi) RESULT(D_psi_i)
  !==========================================
  !
  ! Compute the gradient of the shape function i
  ! on the actual quadrangle at the point of baricentric
  ! coordinates xi
  !  
    IMPLICIT NONE

    CLASS(quadrangle)                       :: e
    INTEGER,                     INTENT(IN) :: i    
    REAL(KIND=8),  DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(e%N_dim) :: D_psi_i
    !-------------------------------------------

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: JJ, inv_J

    INTEGER :: k, l, m
    !---------------------------------------------

    ALLOCATE(    JJ(e%N_dim, e%N_dim) )
    ALLOCATE( inv_J(e%N_dim, e%N_dim) )
        
    JJ = compute_Jacobian( e, xi )
    
    inv_J = inverse(JJ)

    D_Psi_i = MATMUL( inv_J, gradient_ref(e, i, xi) )

    DEALLOCATE( JJ, inv_J )
 
  END FUNCTION gradient
  !====================

  !====================================
  FUNCTION rd_normal(e, i) RESULT(nn_i)
  !====================================

    IMPLICIT NONE

    CLASS(quadrangle)   :: e
    INTEGER, INTENT(IN) :: i

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_i
    !---------------------------------------    

    SELECT CASE(e%Type)

    CASE(QUA_Q1)

       nn_i = rd_normal_QUA_Q1(e, i)

    CASE(QUA_Q2)

       nn_i = rd_normal_QUA_Q2(e, i)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: Unsupported Quadrangle type'
       STOP

    END SELECT
    
  END FUNCTION rd_normal
  !=====================  

  !========================================
  SUBROUTINE gradient_trace(e, if, p_D_phi)
  !========================================
  !
  ! Compute the trace of the gradients of all shape
  ! functions on the face if
  !
    IMPLICIT NONE

    CLASS(quadrangle)                           :: e
    INTEGER,                        INTENT(IN)  :: if
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: p_D_phi
    !-----------------------------------------------------

    REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: xi_f
    INTEGER :: N_points_ele, N_points_face, k, i
    !-----------------------------------------------------

    N_points_ele  = e%N_points
    N_points_face = SIZE(p_D_phi, 2)

    ALLOCATE( xi_f(N_points_face, 2) )

    xi_f = face_trace(e, if)    

    DO k = 1, N_points_face      

       DO i = 1, N_points_ele

          p_D_phi(:, k, i) = gradient( e, i, xi_f(k,:) )

       ENDDO

    ENDDO
    
    DEALLOCATE(xi_f)    

  END SUBROUTINE gradient_trace
  !============================

  !======================================
  FUNCTION face_trace(e, if) RESULT(xi_f)
  !=======================================
  !
  ! Compute the coordinates on the face if
  !
    IMPLICIT NONE

    CLASS(quadrangle)        :: e
    INTEGER,      INTENT(IN) :: if
    
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xi_f
    !---------------------------------------------
     
    SELECT CASE(e%Type)
       
    CASE(QUA_Q1)

       ALLOCATE( xi_f(2, 2) )

       xi_f = face_trace_QUA_Q1(e, if)
       
    CASE(QUA_Q2)

       ALLOCATE( xi_f(3, 2) )

       xi_f = face_trace_QUA_Q2(e, if)
       
    CASE DEFAULT

       WRITE(*,*) 'Unknown Quadrangle type for face trace'
       WRITE(*,*) 'STOP'

    END SELECT 

  END FUNCTION face_trace
  !======================

  !=================================
  FUNCTION DOFs_ref(e) RESULT(x_dof)
  !=================================
  !
  ! Compute the coordinates of the DOFs
  ! on the reference element
  !
    IMPLICIT NONE

    CLASS(quadrangle) :: e
    
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: x_dof
    !-------------------------------------------------
     
    SELECT CASE(e%Type)
       
    CASE(QUA_Q1)

       ALLOCATE( x_dof(4, 2) )

       x_dof = DOFs_ref_QUA_Q1(e)
       
    CASE(QUA_Q2)

       ALLOCATE( x_dof(9, 2) )

       x_dof = DOFs_ref_QUA_Q2(e)
          
    CASE DEFAULT

       WRITE(*,*) 'Unknown Quadrangle type for DOFs'
       STOP

    END SELECT 

  END FUNCTION DOFs_ref
  !====================

  !==========================================
  FUNCTION Compute_Jacobian(e, xi) RESULT(JJ)
  !==========================================

    IMPLICIT NONE

    CLASS(quadrangle)                       :: e
    REAL(KIND=8),  DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8),  DIMENSION(e%N_dim, e%N_dim) :: JJ
    !------------------------------------------------

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: d

    INTEGER :: k, l, m
    !------------------------------------------------

    ALLOCATE( d(e%N_dim, e%N_points) )

    DO m = 1, e%N_points
       d(:, m) = gradient_ref(e, m, xi)
    ENDDO

    ! Construction of the Jacobian matrix
    DO k = 1, e%N_dim

       DO l = 1, e%N_dim

          JJ(l, k) = SUM( e%Coords(k, :) * d(l, :) )

       ENDDO

    ENDDO

    DEALLOCATE( d )
    
  END FUNCTION Compute_Jacobian
  !============================

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFIC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%    QUADRANGLE Q1   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !============================================================
  SUBROUTINE init_faces_QUA_Q1(e, Nodes, Coords, NU_seg, n_ele)
  !============================================================

    IMPLICIT NONE

    CLASS(quadrangle)                        :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    INTEGER,      DIMENSION(:),   INTENT(IN) :: NU_seg
    INTEGER,      DIMENSION(:),   INTENT(IN) :: n_ele
    !-------------------------------------------------

    TYPE(segment), POINTER :: seg

    INTEGER,      DIMENSION(:),   ALLOCATABLE :: loc
    INTEGER,      DIMENSION(:),   ALLOCATABLE :: VV
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: RR

    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: p_D_phi

    INTEGER :: id, if, i, j, k, N_ln, istat
    !-------------------------------------------------

    e%N_faces  = 4 ! # Faces

    ALLOCATE( e%faces(e%N_faces) )

    ! # of local nodes
    N_ln = 2

    ! Dofs on the faces
    ALLOCATE( VV(N_ln), loc(N_ln) )

    ! Coordinates on the faces
    ALLOCATE( RR(SIZE(Coords, 1), N_ln) )

    ALLOCATE( p_D_phi(e%N_dim, N_ln, e%N_points) )

    DO if = 1, e%N_faces

       i = MOD(if, MIN(if+1, e%N_faces)) + 1
       j = MOD(i,  MIN(if+2, e%N_faces)) + 1

       loc = (/ i, j /)

       VV = Nodes(loc)
       RR = Coords(:, loc)

       ALLOCATE(seg, STAT=istat)          
       IF(istat /=0) THEN
          WRITE(*,*) 'ERROR: failed quadrangle allocation'
       ENDIF

       CALL seg%initialize( loc, VV, RR, NU_seg(if), n_ele(if) )

       CALL gradient_trace(e, if, p_D_phi)
       
       CALL seg%face_quadrature(p_D_phi)
          
       e%faces(if)%f => seg
          
    ENDDO

    DEALLOCATE( VV, RR, loc, p_D_phi )

  END SUBROUTINE init_faces_QUA_Q1
  !===============================

  !=====================================================
  FUNCTION basis_function_QUA_Q1(e, i, xi) RESULT(psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(quadrangle)                      :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8) :: psi_i
    !-------------------------------------------

    REAL(KIND=8) :: x, y
    !-------------------------------------------

    x = xi(1); y = xi(2)

    SELECT CASE (i)
       
    CASE(1)

       psi_i = (1.d0 - x)*(1.d0 - y)/4.d0

    CASE(2)

       psi_i = (1.d0 + x)*(1.d0 - y)/4.d0

    CASE(3)

       psi_i = (1.d0 + x)*(1.d0 + y)/4.d0
       
    CASE(4)       
       
       psi_i = (1.d0 - x)*(1.d0 + y)/4.d0

    CASE DEFAULT

       WRITE(*,*) 'ERROR: not supported Dof in Quadrangle  baisis function'
       STOP
       
    END SELECT

  END FUNCTION basis_function_QUA_Q1
  !=================================

  !=====================================================
  FUNCTION gradient_ref_QUA_Q1(e, i, xi) RESULT(D_psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(quadrangle)                      :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(2) :: D_psi_i
    !-------------------------------------------

    REAL(KIND=8) :: x, y
    !-------------------------------------------

    x = xi(1); y = xi(2)

    SELECT CASE(i)

    CASE(1)

       D_psi_i(1) = (-1.d0 + y)/4.d0
       D_psi_i(2) = (-1.d0 + x)/4.d0

    CASE(2)

       D_psi_i(1) = ( 1.d0 - y)/4.d0       
       D_psi_i(2) = (-1.d0 - x)/4.d0

    CASE(3)

       D_psi_i(1) = ( 1.d0 + y)/4.d0       
       D_psi_i(2) = ( 1.d0 + x)/4.d0
 
    CASE(4)

       D_psi_i(1) = (-1.d0 - y)/4.d0       
       D_psi_i(2) = ( 1.d0 - x)/4.d0

    CASE DEFAULT

       WRITE(*,*) 'ERROR: not supported Dof in Quadrangle gradient'
       STOP

    END SELECT  
    
  END FUNCTION gradient_ref_QUA_Q1
  !===============================  
  
  !========================================
  FUNCTION DOFs_ref_QUA_Q1(e) RESULT(x_dof)
  !========================================

    IMPLICIT NONE

    CLASS(quadrangle) :: e

    REAL(KIND=8), DIMENSION(4,2) :: x_dof
    !------------------------------------

    x_dof(1, :) = (/ -1.d0,  -1.d0  /)
    x_dof(2, :) = (/  1.d0,  -1.d0  /)
    x_dof(3, :) = (/  1.d0,   1.d0  /)
    x_dof(4, :) = (/ -1.d0,   1.d0  /)
    
  END FUNCTION DOFs_ref_QUA_Q1
  !===========================

  !=============================================
  FUNCTION face_trace_QUA_Q1(e, if) RESULT(xi_f)
  !=============================================

    IMPLICIT NONE

    CLASS(quadrangle)        :: e
    INTEGER,      INTENT(IN) :: if

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xi_f
    !------------------------------------------------

    ALLOCATE( xi_f(2, 2) )

    SELECT CASE(if)

    CASE(1)

       xi_f(1, :) = (/ 1.d0, -1.d0 /)
       xi_f(2, :) = (/ 1.d0,  1.d0 /)

    CASE(2)

       xi_f(1, :) = (/  1.d0, 1.d0 /)
       xi_f(2, :) = (/ -1.d0, 1.d0 /)
       
    CASE(3)
  
       xi_f(1, :) = (/ -1.d0,  1.d0 /)
       xi_f(2, :) = (/ -1.d0, -1.d0 /)

    CASE(4)
  
       xi_f(1, :) = (/ -1.d0, -1.d0 /)
       xi_f(2, :) = (/  1.d0, -1.d0 /)

    CASE DEFAULT

       WRITE(*,*) 'ERROR: wrong Quadrangle face in face trace Q1'
       STOP

    END SELECT    
    
  END FUNCTION face_trace_QUA_Q1
  !=============================

  !===================================
  SUBROUTINE init_quadrature_QUA_Q1(e)
  !===================================

    IMPLICIT NONE

    CLASS(quadrangle) :: e
    !----------------------------------------------

    e%N_quad = 4

    ALLOCATE( e%phi_q(e%N_points, e%N_quad) )

    ALLOCATE( e%D_phi_q(e%N_dim, e%N_points, e%N_quad) )

    ALLOCATE( e%w_q(e%N_quad   ) )
    ALLOCATE( e%x_q(e%N_quad, 2) )

    ALLOCATE( e%xx_q(e%N_dim, e%N_quad) )
    
    !-------------------
    ! Quadrature formula
    !--------------------------------------
    e%x_q(1,:) = (/ -1.d0/DSQRT(3.d0), -1.d0/DSQRT(3.d0) /)
    e%x_q(2,:) = (/  1.d0/DSQRT(3.d0), -1.d0/DSQRT(3.d0) /)
    e%x_q(3,:) = (/  1.d0/DSQRT(3.d0),  1.d0/DSQRT(3.d0) /)
    e%x_q(4,:) = (/ -1.d0/DSQRT(3.d0),  1.d0/DSQRT(3.d0) /)
    
    e%w_q = 1.d0
    !--------------------------------------
        
    e%xx_q = 0.d0
    
  END SUBROUTINE init_quadrature_QUA_Q1
  !====================================

  !===========================================
  FUNCTION rd_normal_QUA_Q1(e, k) RESULT(nn_k)
  !===========================================

    IMPLICIT NONE

    CLASS(quadrangle)   :: e
    INTEGER, INTENT(IN) :: k

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_k
    !---------------------------------------

    INTEGER :: i, j
    !---------------------------------------

    SELECT CASE(k)
    CASE(1)
       i = 2; j = 4
    CASE(2)
       i = 3; j = 1
    CASE(3)
       i = 4; j = 2
    CASE(4)
       i = 1; j = 3
    CASE DEFAULT
       WRITE(*,*) 'ERROR: wrong nodes for rd_normal'
       STOP       
    END SELECT    

    nn_k(1) = e%Coords(2,i) - e%Coords(2,j)
    nn_k(2) = e%Coords(1,j) - e%Coords(1,i)

  END FUNCTION rd_normal_QUA_Q1
  !============================

  !=========================================
  SUBROUTINE init_Gradiet_Recovery_QUA_Q1(e)
  !=========================================

    IMPLICIT NONE

    CLASS(quadrangle) :: e
    !----------------------

    e%N_rcv = 1

    ALLOCATE( e%D_phi_R(e%N_dim, e%N_points, e%N_rcv) )
    
    ALLOCATE( e%x_R(e%N_rcv, 2) )

    ALLOCATE( e%xx_R(e%N_dim, e%N_rcv) )

    !----------------
    ! Recovery points
    !---------------------------------------------------------
    e%x_R(1,:) = (/ 0.d0, 0.d0 /)
   
    e%xx_R = 0.d0

  END SUBROUTINE init_Gradiet_Recovery_QUA_Q1
  !==========================================

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFIC FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%    QUADRANGLE Q2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !============================================================
  SUBROUTINE init_faces_QUA_Q2(e, Nodes, Coords, NU_seg, n_ele)
  !============================================================

    IMPLICIT NONE

    CLASS(quadrangle)                        :: e
    INTEGER,      DIMENSION(:),   INTENT(IN) :: Nodes
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: Coords
    INTEGER,      DIMENSION(:),   INTENT(IN) :: NU_seg
    INTEGER,      DIMENSION(:),   INTENT(IN) :: n_ele
    !-------------------------------------------------

    TYPE(segment), POINTER :: seg

    INTEGER,      DIMENSION(:),   ALLOCATABLE :: loc
    INTEGER,      DIMENSION(:),   ALLOCATABLE :: VV
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: RR

    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: p_D_phi
    
    INTEGER :: id, if, i, j, k, N_ln, istat
    !-------------------------------------------------

    e%N_faces  = 4 ! # Faces

    ALLOCATE( e%faces(e%N_faces) )

    ! # of local nodes
    N_ln = 3

    ! Dofs on the faces
    ALLOCATE( VV(N_ln), loc(N_ln) )

    ! Coordinates on the faces
    ALLOCATE( RR(SIZE(Coords, 1), N_ln) )

    ALLOCATE( p_D_phi(e%N_dim, N_ln, e%N_points) )

    DO if = 1, e%N_faces

       i = MOD(if, MIN(if+1, e%N_faces)) + 1
       j = MOD(i,  MIN(if+2, e%N_faces)) + 1

       k = i + 4

       loc = (/ i, j, k /)

       VV = Nodes(loc)
       RR = Coords(:, loc)

       ALLOCATE(seg, STAT=istat)          
       IF(istat /=0) THEN
          WRITE(*,*) 'ERROR: failed quadrangle allocation'
       ENDIF

       CALL seg%initialize( loc, VV, RR, NU_seg(if), n_ele(if) )

       CALL gradient_trace(e, if, p_D_phi)
       
       CALL seg%face_quadrature(p_D_phi)
          
       e%faces(if)%f => seg
          
    ENDDO

    DEALLOCATE( VV, RR, loc, p_D_phi )

  END SUBROUTINE init_faces_QUA_Q2
  !===============================

  !=====================================================
  FUNCTION basis_function_QUA_Q2(e, i, xi) RESULT(psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(quadrangle)                      :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8) :: psi_i
    !-------------------------------------------

    REAL(KIND=8) :: x, y, t1, t3
    !-------------------------------------------

    x = xi(1); y = xi(2)

    t1 = x*x; t3 = y*y
     
    SELECT CASE (i)
       
    CASE(1)
 
      psi_i = (t1 - x)*(t3 - y)/4.d0

    CASE(2)

      psi_i = (t1 + x)*(t3 - y)/4.d0

    CASE(3)
 
      psi_i = (t1 + x)*(t3 + y)/4.d0

    CASE(4)       

      psi_i = (t1 - x)*(t3 + y)/4.d0

   CASE(5)

      psi_i = (-t1 + 1.d0)*(t3 - y)/2.d0

   CASE(6)

      psi_i = (t1 + x)*(-t3 + 1.d0)/2.d0

   CASE(7)

      psi_i = (-t1 + 1.d0)*(t3 + y)/2.d0

   CASE(8)

      psi_i = (t1 - x)*(-t3 + 1.d0)/2.d0
      
   CASE(9)

      psi_i = (-t1 + 1.d0)*(-t3 + 1.d0)
            
   CASE DEFAULT
      
       WRITE(*,*) 'ERROR: not supported Dof in Quadrangle  baisis function'
       STOP
       
    END SELECT

  END FUNCTION basis_function_QUA_Q2
  !=================================

  !=====================================================
  FUNCTION gradient_ref_QUA_Q2(e, i, xi) RESULT(D_psi_i)
  !=====================================================

    IMPLICIT NONE

    CLASS(quadrangle)                      :: e
    INTEGER,                    INTENT(IN) :: i
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xi

    REAL(KIND=8), DIMENSION(2) :: D_psi_i
    !-------------------------------------------

    REAL(KIND=8) :: x, y, t1, t2
    !-------------------------------------------

    x = xi(1); y = xi(2)

    t1 = x*x; t2 = y*y
    
    SELECT CASE(i)

    CASE(1)
     
       D_psi_i(1) = (x - 1.d0/2.d0)*(t2 - y)/2.d0
       D_psi_i(2) = (t1 - x)*(y - 1.d0/2.d0)/2.d0

      
    CASE(2)

      D_psi_i(1) = (x + 1.d0/2.d0)*(t2-y)/2.d0
      D_psi_i(2) = (t1 + x)*(y - 1.d0/2.d0)/2.d0
             
    CASE(3)

      D_psi_i(1) = (x + 1.d0/2.d0)*(t2 + y)/2.d0
      D_psi_i(2) = (t1 + x)*(y + 1.d0/2.d0)/2.d0

    CASE(4)

      D_psi_i(1) = (x - 1.d0/2.d0)*(t2 + y)/2.d0
      D_psi_i(2) = (t1 - x)*(y + 1.d0/2.d0)/2.d0

   CASE(5)

      D_psi_i(1) = -x*(t2 - y)
      D_psi_i(2) = (-t1 + 1)*(y - 1.d0/2.d0)

   CASE(6)

      D_psi_i(1) = (x + 1.d0/2.d0)*(-t2 + 1)
      D_psi_i(2) = -(t1 + x)*y

   CASE(7)

      D_psi_i(1) = -x*(t2 + y)
      D_psi_i(2) = (-t1 + 1)*(y + 1.d0/2.d0)

   CASE(8)

      D_psi_i(1) = (x - 1.d0/2.d0)*(-t2 + 1)
      D_psi_i(2) = -(t1 - x)*y

   CASE(9)

      D_psi_i(1) = -2*x*(-t2 + 1)
      D_psi_i(2) = -2*(-t1 + 1)*y
      
   CASE DEFAULT
       
       WRITE(*,*) 'ERROR: not supported Dof in Quadrangle gradient'
       STOP

    END SELECT  
    
  END FUNCTION gradient_ref_QUA_Q2
  !===============================

  !========================================
  FUNCTION DOFs_ref_QUA_Q2(e) RESULT(x_dof)
  !========================================

    IMPLICIT NONE

    CLASS(quadrangle) :: e

    REAL(KIND=8), DIMENSION(9,2) :: x_dof
    !------------------------------------

    x_dof(1, :) = (/ -1.d0,  -1.d0  /)
    x_dof(2, :) = (/  1.d0,  -1.d0  /)
    x_dof(3, :) = (/  1.d0,   1.d0  /)
    x_dof(4, :) = (/ -1.d0,   1.d0  /)
    x_dof(5, :) = (/  0.d0,  -1.d0  /)
    x_dof(6, :) = (/  1.d0,   0.d0  /)
    x_dof(7, :) = (/  0.d0,   1.d0  /)
    x_dof(8, :) = (/ -1.d0,   0.d0  /)
    x_dof(9, :) = (/  0.d0,   0.d0  /)    
    
  END FUNCTION DOFs_ref_QUA_Q2
  !===========================

  !===========================================
  FUNCTION rd_normal_QUA_Q2(e, k) RESULT(nn_k)
  !===========================================

    IMPLICIT NONE

    CLASS(quadrangle)   :: e
    INTEGER, INTENT(IN) :: k

    REAL(KIND=8), DIMENSION(e%N_dim) :: nn_k
    !---------------------------------------

    INTEGER :: i, j
    !---------------------------------------

    SELECT CASE(k)
    CASE(1)
       i = 2; j = 4
    CASE(2)
       i = 3; j = 1
    CASE(3)
       i = 4; j = 2
    CASE(4)
       i = 1; j = 3
    CASE(5)
       i = 2; j = 1
    CASE(6)
       i = 3; j = 2
    CASE(7)
       i = 4; j = 3
    CASE(8)
       i = 1; j = 4
    CASE(9)
       i = 1; j = 1        
    CASE DEFAULT
       WRITE(*,*) 'ERROR: wrong nodes for rd_normal'
       STOP       
    END SELECT    

    nn_k(1) = e%Coords(2,i) - e%Coords(2,j)
    nn_k(2) = e%Coords(1,j) - e%Coords(1,i)

  END FUNCTION rd_normal_QUA_Q2
  !============================

  !=============================================
  FUNCTION face_trace_QUA_Q2(e, if) RESULT(xi_f)
  !=============================================

    IMPLICIT NONE

    CLASS(quadrangle)        :: e
    INTEGER,      INTENT(IN) :: if

    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xi_f
    !----------------------------------------------

    ALLOCATE( xi_f(3, 2) )

    SELECT CASE(if)

    CASE(1)

       xi_f(1, :) = (/ 1.d0,  -1.d0 /)
       xi_f(2, :) = (/ 1.d0,   1.d0 /)
       xi_f(3, :) = (/ 1.0d0,  0.d0 /)

    CASE(2)

       xi_f(1, :) = (/  1.d0,  1.d0 /)
       xi_f(2, :) = (/ -1.d0,  1.d0 /)
       xi_f(3, :) = (/  0.d0,  1.d0 /)

    CASE(3)
  
       xi_f(1, :) = (/ -1.d0,  1.d0  /)
       xi_f(2, :) = (/ -1.d0, -1.d0  /)
       xi_f(3, :) = (/ -1.d0,  0.0d0 /)

     CASE(4)

       xi_f(1, :) = (/ -1.d0, -1.d0 /)
       xi_f(2, :) = (/  1.d0, -1.d0 /)
       xi_f(3, :) = (/  0.d0, -1.d0 /)
           
    CASE DEFAULT

       WRITE(*,*) 'ERROR: wrong Quadrangle face in face trace Q2'
       STOP

    END SELECT    
    
  END FUNCTION face_trace_QUA_Q2
  !=============================

  !===================================
  SUBROUTINE init_quadrature_QUA_Q2(e)
  !===================================

    IMPLICIT NONE

    CLASS(quadrangle) :: e
    !----------------------------------------------

    e%N_quad = 9

    ALLOCATE( e%phi_q(e%N_points, e%N_quad) )

    ALLOCATE( e%D_phi_q(e%N_dim, e%N_points, e%N_quad) )

    ALLOCATE( e%w_q(e%N_quad   ) )
    ALLOCATE( e%x_q(e%N_quad, 2) )

    ALLOCATE( e%xx_q(e%N_dim, e%N_quad) )
    
    !-------------------
    ! Quadrature formula
    !--------------------------------------
    e%x_q(1,:) = (/ -DSQRT(0.6d0), -DSQRT(0.6d0) /)
    e%x_q(2,:) = (/        0.d0 ,  -DSQRT(0.6d0) /)
    e%x_q(3,:) = (/  DSQRT(0.6d0), -DSQRT(0.6d0) /)
    e%x_q(4,:) = (/ -DSQRT(0.6d0),         0.d0  /)
    e%x_q(5,:) = (/        0.d0 ,          0.d0  /)
    e%x_q(6,:) = (/  DSQRT(0.6d0),         0.d0  /)
    e%x_q(7,:) = (/ -DSQRT(0.6d0),  DSQRT(0.6d0) /)
    e%x_q(8,:) = (/        0.d0 ,   DSQRT(0.6d0) /)
    e%x_q(9,:) = (/  DSQRT(0.6d0),  DSQRT(0.6d0) /)
    
    e%w_q(1) = 25.d0/81.d0
    e%w_q(2) = 40.d0/81.d0
    e%w_q(3) = 25.d0/81.d0
    e%w_q(4) = 40.d0/81.d0
    e%w_q(5) = 64.d0/81.d0
    e%w_q(6) = 40.d0/81.d0
    e%w_q(7) = 25.d0/81.d0
    e%w_q(8) = 40.d0/81.d0
    e%w_q(9) = 25.d0/81.d0
    !--------------------------------------
        
    e%xx_q = 0.d0
    
  END SUBROUTINE init_quadrature_QUA_Q2
  !====================================

  !=========================================
  SUBROUTINE init_Gradiet_Recovery_QUA_Q2(e)
  !=========================================

    IMPLICIT NONE

    CLASS(quadrangle) :: e
    !----------------------

    e%N_rcv = 4

    ALLOCATE( e%D_phi_R(e%N_dim, e%N_points, e%N_rcv) )
    
    ALLOCATE( e%x_R(e%N_rcv, 2) )

    ALLOCATE( e%xx_R(e%N_dim, e%N_rcv) )

    !----------------
    ! Recovery points
    !---------------------------------------------------------
    e%x_R(1,:) = (/ -1.d0/DSQRT(3.d0), -1.d0/DSQRT(3.d0) /)
    e%x_R(2,:) = (/  1.d0/DSQRT(3.d0), -1.d0/DSQRT(3.d0) /)
    e%x_R(3,:) = (/ -1.d0/DSQRT(3.d0),  1.d0/DSQRT(3.d0) /)
    e%x_R(4,:) = (/  1.d0/DSQRT(3.d0),  1.d0/DSQRT(3.d0) /)

    e%xx_R = 0.d0

  END SUBROUTINE init_Gradiet_Recovery_QUA_Q2
  !==========================================

END MODULE Quadrangle_Class

