MODULE time_integration

  USE Element_Class
  
  USE geometry,          ONLY: N_dofs, N_elements, elements
  USE init_problem,      ONLY: pb_name, pb_type, visc, time_int, CFL
  USE models,            ONLY: strong_bc
  USE space_integration

  IMPLICIT NONE
  PRIVATE

  REAL(KIND=8) :: res_0
  
  PUBLIC :: time_advance, Nodal_Pe

CONTAINS
  
  !=========================================
  SUBROUTINE time_advance(ite, uu, rhs, res)
  !=========================================

    IMPLICIT NONE

    INTEGER,                    INTENT(IN)    :: ite
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: uu
    REAL(KIND=8), DIMENSION(:), INTENT(OUT)   :: rhs
    REAL(KIND=8),               INTENT(OUT)   :: res
    !-----------------------------------------------

    REAL(KIND=8), DIMENSION(N_dofs) :: Dt_V
    
    INTEGER :: UNIT, ierror
    !-----------------------------------------------

    CALL compute_rhs(uu, rhs, Dt_V)
    
    IF ( time_int == 0 ) THEN

       !---------------
       ! Esplicit Euler
       !----------------------------------
       uu = uu - Dt_V*rhs

    ELSEIF ( time_int == 1 ) THEN

       !---------------
       ! Implicit Euler
       !---------------------------------
       STOP
       
    ENDIF
   
    !-----------------------
    ! Normalized L2 Residual
    !-----------------------------------------------
    res = SQRT(SUM(rhs*rhs)) / REAL(N_dofs)

    IF(ite == 1) res_0 = res
    res = res/res_0

    !------------------------
    ! Save convergence strory
    !------------------------------------------------
    UNIT = 3

    OPEN(UNIT, FILE = 'convergence.'//TRIM(ADJUSTL(pb_name)), &
         ACTION = 'WRITE', POSITION = 'APPEND', IOSTAT = ierror)
            
    IF(ierror /= 0) THEN
       
       WRITE(*,*) 'ERROR: Impossible to open the file converegence'
       WRITE(*,*) 'STOP'
       STOP
       
    ENDIF
    
    WRITE(UNIT, 500) ite, res
    
    CLOSE(UNIT)
     
    !--------------------
    ! Data on the screan
    !------------------------------------------------
    WRITE(*, 600) ite, MINVAL(uu), MAXVAL(uu)
    WRITE(*, 601) res, res_0
    WRITE(*, *)

500 FORMAT(I5, 1x, e24.16)      
600 FORMAT('Iteration # ', I5, '    uu min = ', F10.6, ', uu max =', F10.6)
601 FORMAT('Residual = ', e16.8, '    Residual_0 = ' e16.8)

  END SUBROUTINE time_advance
  !==========================

  !=================================
  FUNCTION Nodal_Pe(uu) RESULT(n_Pe)
  !=================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: uu

    REAL(KIND=8), DIMENSION(SIZE(uu)) :: n_Pe
    !---------------------------------------------

    REAL(KIND=8), DIMENSION(SIZE(uu)) :: Den
    
    TYPE(element) :: loc_ele

    REAL(KIND=8) :: ele_Pe

    INTEGER :: je
    !---------------------------------------------

    n_Pe = 0.d0;  Den = 0.d0

    DO je = 1, N_elements

       loc_ele = elements(je)%p

       ele_Pe = Local_Pe(loc_ele, uu(loc_ele%NU))

       n_Pe(loc_ele%NU) = n_Pe(loc_ele%NU) + &
                          ele_Pe * loc_ele%volume

       Den(loc_ele%NU) = Den(loc_ele%NU) + loc_ele%volume
       
    ENDDO

    n_Pe = n_Pe / Den

  END FUNCTION Nodal_Pe  
  !====================

END MODULE time_integration
