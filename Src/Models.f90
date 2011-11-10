MODULE models

  USE geometry,      ONLY: N_dim, N_elements, elements

  IMPLICIT NONE

  PRIVATE

   !===============================================
   INTEGER, PARAMETER :: LINEAR_ADVECTION   = 0, &
                         ROTATION           = 1, &
                         BURGER             = 2, &
                         LIN_VISC_ADVEC     = 3, &
                         SMITH_HUTTON       = 4, &
                         PERAIRE_ADVEC_DIFF = 5, &
                         ADVECTION_SOURCE   = 6
   !===============================================
   
   REAL(KIND=8), PARAMETER :: PI = DACOS(-1.d0)                         
   REAL(KIND=8), PARAMETER :: theta = 100.d0
   !===============================================

   PUBLIC :: advection_speed, flux_function,         & 
             exact_solution, strong_bc, source_term, &
             detect_source
   !===============================================

 CONTAINS
   
   !=====================================================
   FUNCTION advection_speed(type_pb, uu, x, y)  RESULT(a)
   !=====================================================
   !
   ! Advection speed of the bidimensional 
   ! problem written in quasi-linear form
   !
   IMPLICIT NONE
   
      INTEGER,      INTENT(IN) :: type_pb
      REAL(KIND=8), INTENT(IN) :: uu
      REAL(KIND=8), INTENT(IN) :: x, y
      
      REAL(KIND=8), DIMENSION(2) :: a
      !------------------------------------
      
      SELECT CASE (type_pb)
      
         CASE(LINEAR_ADVECTION)
         
            a(1) = 0.d0
            a(2) = 1.d0
      
         CASE(ROTATION)
         
            a(1) = -y
            a(2) =  x
      
         CASE(BURGER)
         
            a(1) = uu 
            a(2) = 1.d0
             
	 CASE(LIN_VISC_ADVEC)
	 
	    a(1) = 0.d0
	    a(2) = 1.d0

         CASE(SMITH_HUTTON)
         
            a(1) =  2.d0*(1 - x**2)*y
            a(2) = -2.d0*(1 - y**2)*x

         CASE(PERAIRE_ADVEC_DIFF)

            a(1) = 1.d0
            a(2) = 2.d0*x

         CASE(ADVECTION_SOURCE)

            a(1) = y
            a(2) = 1.d0 - x
            
         CASE DEFAULT
         
            WRITE(*,*) 'Problem of unknow type.'
            WRITE(*,*) 'STOP!'
            
            STOP
      
      END SELECT
   
   END fUNCTION advection_speed
   !===========================

   !======================================================
   FUNCTION flux_function(type_pb, uu, x, y)  RESULT(flux)
   !======================================================
   !
   ! Flux of the bidimensional problem 
   ! written in quasi-linear form
   !
   IMPLICIT NONE
   
      INTEGER,      INTENT(IN) :: type_pb
      REAL(KIND=8), INTENT(IN) :: uu
      REAL(KIND=8), INTENT(IN) :: x, y
      
      REAL(KIND=8), DIMENSION(2) :: flux
      !------------------------------------
      
      SELECT CASE (type_pb)
      
         CASE(LINEAR_ADVECTION)
         
            flux(1) = 0.d0
            flux(2) = uu
      
         CASE(ROTATION)
         
            flux(1) = -y*uu
            flux(2) =  x*uu
      
         CASE(BURGER)
         
            flux(1) = 0.5d0 * uu**2 
            flux(2) = uu

	 CASE(LIN_VISC_ADVEC)
	 
	    flux(1) = 0.d0
            flux(2) = uu
                  
         CASE(SMITH_HUTTON)
         
            flux(1) =  2.d0*(1 - x**2)*y * uu
            flux(2) = -2.d0*(1 - y**2)*x * uu

         CASE(PERAIRE_ADVEC_DIFF)
            
            flux(1) = 1.d0*uu
            flux(2) = 2.d0*x*uu

         CASE(ADVECTION_SOURCE)

            flux(1) = y*uu
            flux(2) = (1.d0 - x)*uu
      
         CASE DEFAULT
         
            WRITE(*,*) 'Problem of unknow type.'
            WRITE(*,*) 'STOP!'
            
            STOP
      
      END SELECT
   
   END FUNCTION flux_function
   !=========================

   !================================================
   FUNCTION source_term(type_pb, uu, x, y) RESULT(S)
   !================================================
   IMPLICIT NONE

      INTEGER,      INTENT(IN) :: type_pb
      REAL(KIND=8), INTENT(IN) :: uu
      REAL(KIND=8), INTENT(IN) :: x, y
      
      REAL(KIND=8) :: S
      !------------------------------------

      REAL(KIND=8) :: r
      !------------------------------------
      
      SELECT CASE (type_pb)

      CASE(ADVECTION_SOURCE)

         r = SQRT((x - 0.5d0)**2 + (y - 0.3)**2)

         IF ( r <= 0.25 ) THEN
            S = 10.d0
         ELSE
            S = 0.d0
         ENDIF
         
         
      CASE DEFAULT

         WRITE(*,*) 'Problem of unknow type.'         
         WRITE(*,*) 'STOP!'
         
         STOP
         
      END SELECT
      
   END FUNCTION source_term
   !=======================

   !============================================
   FUNCTION detect_source(type_pb) RESULT(logic)
   !============================================

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: type_pb

     LOGICAL :: logic
     !------------------------------

     logic = .FALSE.

     IF (type_pb == ADVECTION_SOURCE) logic = .TRUE.
     
   END FUNCTION detect_source
   !=========================

   
   !===========================================
   SUBROUTINE strong_bc(type_pb, visc, uu, rhs)
   !===========================================
   !
   ! Impose the bondary conditions in the strong way,
   ! according to the chosen type of the problem
   !
   IMPLICIT NONE
   
      INTEGER,                    INTENT(IN)    :: type_pb
      REAL(KIND=8),               INTENT(IN)    :: visc
      REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: uu
      REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: rhs
      !======================================================
      
      INTEGER, DIMENSION(:), ALLOCATABLE :: is, is_loc
      
      REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: u_l, rhs_l
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: coord
      
      INTEGER :: n_loc, Nv
      
      LOGICAL :: b_flag
      
      INTEGER :: jt, i, k
      !=======================================================
      
      DO jt = 1, N_elements

         Nv = elements(jt)%p%N_points

         ALLOCATE( is(Nv), is_loc(Nv), u_l(Nv), rhs_l(Nv) )
         ALLOCATE( coord(N_dim, Nv) )
      
         is = elements(jt)%p%NU
         
         DO i = 1, N_dim
            coord(i, :) = elements(jt)%p%Coords(i, :)
         ENDDO

           u_l(:) =  uu(is)
         rhs_l(:) = rhs(is)

          n_loc = 0
         is_loc = 0
         b_flag = .FALSE.
         
         SELECT CASE(type_pb)
      
         CASE(LINEAR_ADVECTION)
      
            CALL bc_linadv()
               
         CASE(ROTATION)
            
            CALL bc_rotation()
               
         CASE(BURGER)
            
            CALL bc_burger()
               
	 CASE(LIN_VISC_ADVEC)
	    
	    CALL bc_linviscadv()

         CASE(SMITH_HUTTON)
            
            CALL bc_smithutton()

         CASE(PERAIRE_ADVEC_DIFF)

            CALL bc_peraire_advec_diff()

         CASE(ADVECTION_SOURCE)

            CALL bc_advection_source()
                        
         CASE DEFAULT
         
            WRITE(*,*) 'Problem of unknow type.'
            WRITE(*,*) 'STOP!'   
            STOP
            
         END SELECT
         
         IF (b_flag) THEN
               
            DO i = 1, n_loc
                uu( is(is_loc(i)) ) =   u_l( is_loc(i) )            
               rhs( is(is_loc(i)) ) = rhs_l( is_loc(i) )
            ENDDO
                     
         ENDIF

         DEALLOCATE( is, is_loc, coord, u_l, rhs_l )
      
      ENDDO

!=======
CONTAINS
!=======

      !......................
      SUBROUTINE bc_linadv()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( ABS(coord(1, k)) <= 0.d0 .AND. & 
                 ABS(coord(2, k)) >= 0.d0) THEN
               u_l(k) = initial_profile(coord(1, k))
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO
          
         DO k = 1, Nv
            IF ( (ABS(coord(1, k) - 1.d0)) <= 0.d0) THEN
               u_l(k) = initial_profile(coord(1, k))
               rhs_l(k) = 0.d0              
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO
          
         DO k = 1, Nv
            IF ( ABS(coord(2, k)) <= 0.d0 ) THEN
               u_l(k) = initial_profile(coord(1, k))
               rhs_l(k) = 0.d0              
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO        
      
      END SUBROUTINE bc_linadv
      !........................

      !........................
      SUBROUTINE bc_rotation()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( (ABS(coord(1, k) - 1.d0)) <= 0.d0 ) THEN
               u_l(k) = 0.d0
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO
                    
         DO k = 1, Nv
            IF ( ABS(coord(2, k)) <= 0.d0 ) THEN
               IF (coord(1, k) >= 0.25d0 .AND. coord(1, k) <= 0.75d0) THEN
                  u_l(k) = DCOS(2.d0*PI*coord(1, k))**2
               ELSE
                  u_l(k) = 0
               ENDIF
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.           
            ENDIF            
         ENDDO        
      
      END SUBROUTINE bc_rotation
      !........................
      
      !......................
      SUBROUTINE bc_burger()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( ABS(coord(2, k)) <= 0.0 ) THEN
               u_l(k) = 1.5d0 - 2.d0*coord(1, k)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.  
            ENDIF
         ENDDO
                    
         DO k = 1, Nv
            IF ((ABS(coord(1, k)  - 1.d0)) <= 0.0 ) THEN
               u_l(k) = -0.5d0
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO 
                
         DO k = 1, Nv
            IF ( ABS(coord(1, k)) <= 0.d0 .AND. &
                 ABS(coord(2, k)) >= 0.d0) THEN
               u_l(k) = 1.5d0 - 2.d0*coord(1, k)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.         
            ENDIF            
         ENDDO   
             
      END SUBROUTINE bc_burger
      !........................

      !.........................
      SUBROUTINE bc_linviscadv()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( ABS(coord(1, k)) <= 0.d0 .AND. & 
                 ABS(coord(2, k)) >= 0.d0 ) THEN
               u_l(k) = visc_advection(coord(1, k), coord(2, k), visc)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF
         ENDDO
          
         DO k = 1, Nv           
            IF ( (ABS(coord(1, k)  - 1.d0)) <= 0.d0 ) THEN
               u_l(k) = visc_advection(coord(1, k), coord(2, k), visc)
               rhs_l(k) = 0.d0              
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO
          
         DO k = 1, Nv
            IF ( ABS(coord(2, k)) <= 0.d0 ) THEN
               u_l(k) =  visc_advection(coord(1, k), coord(2, k), visc)
               rhs_l(k) = 0.d0              
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO   
             
      END SUBROUTINE bc_linviscadv
      !...........................
      
      !.........................
      SUBROUTINE bc_smithutton()
      !
      ! Rectangle [-1,1]x[0,1]
      !
      IMPLICIT NONE

         DO k = 1, Nv
            IF ( ABS(coord(2, k)) <= 0.d0 .AND. &
                 coord(1, k) <= 0.d0 ) THEN
	       u_l(k) = 1 + TANH(theta*(2.d0*coord(1, k) + 1.d0))
	       rhs_l(k) = 0.d0
	       n_loc = n_loc + 1              
	       is_loc(n_loc) = k
	       b_flag = .TRUE.               
            ENDIF
         ENDDO
                    
         DO k = 1, Nv
            IF ((ABS(coord(2, k) - 1.d0)) <= 0.d0 ) THEN
               u_l(k) = 1.d0 - TANH(theta)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.             
            ENDIF            
         ENDDO 
                
         DO k = 1, Nv
            IF ( (ABS(coord(1, k) - 1.d0)) <= 0.d0 .AND. &
                  ABS(coord(2, k))         >= 0.d0 ) THEN
               u_l(k) = 1.d0 - TANH(theta)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.         
            ENDIF            
         ENDDO   
         
         DO k = 1, Nv
            IF ( (ABS(coord(1, k) + 1.d0)) <= 0.d0 .AND. &
                  ABS(coord(2, k) - 1.d0) >= 0.d0 ) THEN
               u_l(k) = 1.d0 - TANH(theta)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1              
               is_loc(n_loc) = k
               b_flag = .TRUE.         
            ENDIF            
         ENDDO
         
      END SUBROUTINE bc_smithutton
      !............................

      !..................................
      SUBROUTINE  bc_peraire_advec_diff()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE
     
        DO k = 1, Nv
            IF ( ABS(coord(1, k)) <= 0.d0 .AND. & 
                 ABS(coord(2, k)) >= 0.d0) THEN
               u_l(k) = 1.d0 - coord(2, k)
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO

        DO k = 1, Nv
            IF ( ABS(coord(1, k)) >=  0.d0 .AND. & 
                 ABS(coord(2, k)) <= 0.d0) THEN
               u_l(k) = coord(1, k) - 1.d0
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO         
           
      END SUBROUTINE bc_peraire_advec_diff
      !...................................

      !...............................
      SUBROUTINE bc_advection_source()
      !
      ! Square [0,1]x[0,1]
      !
      IMPLICIT NONE

        DO k = 1, Nv
            IF ( ABS(coord(1, k)) <= 0.d0 .AND. & 
                 ABS(coord(2, k)) >= 0.d0) THEN
               u_l(k) = 0.d0
               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO       

        DO k = 1, Nv
            IF ( ABS(coord(1, k)) >= 0.d0 .AND. & 
                 ABS(coord(2, k)) <= 0.d0) THEN

               IF ( coord(1, k) >= 0.3d0 .AND. &
                    coord(1, k) <= 0.8d0 ) THEN
                  u_l(k) = -5.d0
               ELSE                  
                  u_l(k) = -0.d0                  
               ENDIF               

               rhs_l(k) = 0.d0
               n_loc = n_loc + 1
               is_loc(n_loc) = k
               b_flag = .TRUE.
            ENDIF
         ENDDO
      
      END SUBROUTINE bc_advection_source
      !.................................
         
   END SUBROUTINE strong_bc
   !=======================
   
   !========================================================
   FUNCTION exact_solution(type_pb, coord, nu) RESULT(uu_ex)
   !========================================================
   !
   ! Compute the exact solution of the 2d problem, when possible
   !
   IMPLICIT NONE
   
      INTEGER,                    INTENT(IN) :: type_pb
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: coord
      REAL(KIND=8),               INTENT(IN) :: nu
      
      REAL(KIND=8) :: uu_ex
      !-------------------------------------------------
      
      REAL(KIND=8) :: x, y, r, phi    
      !-------------------------------------------------

      x = coord(1)
      y = coord(2)
      
      SELECT CASE (type_pb)
      
      CASE(LINEAR_ADVECTION)
                  
         uu_ex = initial_profile(x)
                                 
      CASE(ROTATION)
                  
         r = SQRT(x**2 + y**2)         
         IF( r >= 0.25 .AND. r <= 0.75 ) THEN
            uu_ex = DCOS(2*PI*r)**2
         ELSE            
            uu_ex = 0.d0            
         ENDIF
         
            
      CASE(BURGER)
                  
         IF (y >= 0.5) THEN                                                                
            IF (-2.d0*(x - 3.d0/4.d0) + (y - 0.5) <= 0) THEN
               uu_ex = -0.5
            ELSE
               uu_ex = 1.5   
            ENDIF               
         ELSE                        
            uu_ex = MAX( -0.5, MIN( 1.5, ((x - 3.d0/4.d0)/(y - 0.5)) ) )                
         ENDIF


      CASE(LIN_VISC_ADVEC)
      
         IF (nu > 0.d0) THEN    
            uu_ex = visc_advection(x, y, nu)	          
         ELSE	    
            WRITE(*,*) 'ERROR: LIN_VISC_ADVEC'            
            WRITE(*,*) 'STOP!'
            STOP       
         ENDIF                 
            
      CASE(SMITH_HUTTON)
                       
         phi = -(1.d0 - x**2)*(1 - y**2)         
         uu_ex = 1.d0 + TANH(theta*(1.d0 - 2.d0*SQRT(1.d0 + phi)))
                        
      CASE DEFAULT
         
         uu_ex = 0.0
      
      END SELECT
   
   END FUNCTION exact_solution
   !==========================   

   !====================================
   FUNCTION initial_profile(x) RESULT(f)
   !====================================
   !
   ! Initial profile of the solution used as 
   ! inlet boundary condition in some problems
   !
   IMPLICIT NONE
   
      REAL(KIND=8), INTENT(IN) :: x
      
      REAL(KIND=8) :: f
      !--------------------------------------
      
      REAL(KIND=8), PARAMETER :: gamma = 100
      REAL(KIND=8) :: a, b      
      !--------------------------------------

      a = 1.d0
      b = a*gamma*(x - 0.5)**2
      
      f = EXP(-b)*SQRT(a)   
   
   END FUNCTION initial_profile
   !===========================
   
   !==========================================
   FUNCTION visc_advection(x, y, nu) RESULT(f)
   !==========================================
   !
   ! Solution of the viscous linear advection problem
   !
   IMPLICIT NONE
   
      REAL(KIND=8),               INTENT(IN) :: x
      REAL(KIND=8),               INTENT(IN) :: y
      REAL(KIND=8),               INTENT(IN) :: nu
      
      REAL(KIND=8) :: f
      !--------------------------------------------
      
      IF (nu > 0.d0) THEN
	 f = -DCOS(2.d0*PI*x)* & 
	     EXP(0.5d0*y*( 1 - SQRT(1.d0 + 16.d0*(PI**2)*nu**2) )/nu)
      ELSE
         WRITE(*,*) 'ERROR: Visc zero'
         WRITE(*,*) 'STOP!'
         STOP
      ENDIF
   
    END FUNCTION visc_advection
    !==========================

END MODULE models