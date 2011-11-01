MODULE init_problem

  !USE Element_Class
  !USE Geometry,       ONLY: N_dofs, element
  !USE Models,         ONLY: strong_bc, detect_source

  IMPLICIT NONE

  PRIVATE

  !=======================================
  CHARACTER(len=64) :: pb_name
  INTEGER           :: order
  INTEGER           :: scheme_type
  INTEGER           :: time_int
  INTEGER           :: pb_type
  INTEGER           :: ite_max      
  REAL(KIND=8)      :: CFL
  REAL(KIND=8)      :: toll_res      
  REAL(KIND=8)      :: visc
  LOGICAL           :: with_source
  !=======================================

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: theta
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: theta_t

  !=======================================

  PUBLIC :: read_param !, initialization
  PUBLIC :: pb_name, order,              &
            scheme_type, time_int,       &
            pb_type, CFL, visc, ite_max, &
            toll_res, with_source,       &
            theta, theta_t
  !=======================================

CONTAINS

  !=======================================
  SUBROUTINE read_param (unit, param_file)
  !=======================================
  !
  ! Read the parameters of the simulation from the input file
  !
    IMPLICIT NONE

    INTEGER,           INTENT(IN) :: unit
    CHARACTER(len=64), INTENT(IN) :: param_file
    !===========================================

    LOGICAL :: is_visc
    INTEGER :: ierror
     
    OPEN(unit, FILE = param_file, ACTION = 'READ', IOSTAT = ierror)
    IF (ierror /= 0 ) THEN
       WRITE(*,*) 'ERROR: Impossible to open param file', TRIM(ADJUSTL(param_file))
       WRITE(*,*) 'STOP!'    
       STOP
    ENDIF
      
    READ(unit, *) pb_name
    READ(unit, *) order
    READ(unit, *) scheme_type
    READ(unit, *) time_int
    READ(unit, *) pb_type
    READ(unit, *) ite_max
    READ(unit, *) toll_res
    READ(unit, *) CFL
    READ(unit, *) is_visc
    IF(is_visc) THEN 
       READ(unit, *) visc
    ELSE
       visc = 0.d0
    ENDIF

    CLOSE (unit)
      
    WRITE(*,*)
    WRITE(*,*) 'Problem name: ', pb_name
    WRITE(*,*) 'Order:', order
    WRITE(*,*) 'Num scheme:', scheme_type
    WRITE(*,*) 'Time integration:', time_int
    WRITE(*,*) 'Problem type:', pb_type
    WRITE(*,*) 'Ite max:', ite_max
    WRITE(*,*) 'Residual:', toll_res
    WRITE(*,*) 'CFL:', CFL
    WRITE(*,*) 'Viscous?', is_visc
    IF(is_visc) WRITE(*,*) 'Viscosity', visc
   
  END SUBROUTINE read_param
  !========================

END MODULE init_problem
