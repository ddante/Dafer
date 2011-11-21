MODULE init_problem

  USE Geometry,       ONLY: N_dofs, elements,  N_elements
  USE Models,         ONLY: strong_bc, detect_source,  exact_solution

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

  PUBLIC :: read_param, initialization
  PUBLIC :: pb_name, order,              &
            scheme_type, time_int,       &
            pb_type, CFL, visc, ite_max, &
            toll_res, with_source,       &
            theta, theta_t
  !=======================================

CONTAINS
  
   !=================================
   SUBROUTINE initialization(uu, rhs)
   !=================================
   !
   ! Initialize the solution at the first time step
   !
   IMPLICIT NONE
   
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: uu
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: rhs
      !-------------------------------------------------------------
      
integer :: jt, i

      INTEGER :: ierror, UNIT
      !-------------------------------------------------------------

      ! Solution and RHS
      ALLOCATE( uu(N_dofs), rhs(N_dofs) )

      ! Shock detector
      ALLOCATE( theta(N_dofs), theta_t(N_dofs) )
      
      uu  = 0.d0;  rhs = 0.d0

      theta = 1.d0;  theta_t = 10.d0

      with_source = detect_source(pb_type)

DO jt = 1, N_elements
DO i = 1, SIZE(elements(jt)%p%NU)
uu(elements(jt)%p%NU(i)) = exact_solution(pb_type, elements(jt)%p%Coords(:, i), visc)
ENDDO
ENDDO
     
      CALL strong_bc(pb_type, visc, uu, rhs)

      ! Delete a possible previous convergence history...
      UNIT = 4
      
      OPEN(UNIT, FILE = 'convergence.'//TRIM(ADJUSTL(pb_name)), & 
          STATUS= 'REPLACE', IOSTAT = ierror)

      IF(ierror /= 0) THEN
         WRITE(*,*) 'ERROR: Impossible to open the file converegence'
         WRITE(*,*) 'STOP'
      
         STOP
      ENDIF     
      CLOSE (UNIT)
      
      ! ... and error file
      UNIT = 9
      
      OPEN(UNIT, FILE = 'error.'//TRIM(ADJUSTL(pb_name)), &
           STATUS= 'REPLACE', IOSTAT = ierror)

      IF(ierror /= 0) THEN
         WRITE(*,*) 'ERROR: Impossible to open the file error'
         WRITE(*,*) 'STOP'
      
         STOP
      ENDIF     
      CLOSE (UNIT)
   
   END SUBROUTINE initialization
   !============================

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
