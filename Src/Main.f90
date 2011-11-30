PROGRAM main

  USE init_problem,      ONLY: read_param, initialization, order, &
                               ite_max, toll_res

  USE geometry,          ONLY: read_gmshMesh, init_elements

  USE time_integration

  USE post_pro

  IMPLICIT NONE

  !------------------------------------------------
  CHARACTER(len=64) :: param_file, mesh_file

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: uu
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhs

  REAL(KIND=8) :: res

  INTEGER :: ite, UNIT
  !------------------------------------------------

  !--------------------
  !Read the imput files
  !------------------------------------------------  
  IF (command_argument_count() < 2) THEN
      WRITE(*,*) 'ERROR: No file param and/or mesh.'
      WRITE(*,*) 'STOP!'
          
       STOP
  ENDIF

  CALL get_command_argument(1, param_file)
  CALL get_command_argument(2, mesh_file)
   
 ! File param   
 UNIT = 1
 CALL read_param(UNIT, param_file)

 ! File mesh
 UNIT = 2
 CALL read_gmshMesh(UNIT, mesh_file)

 !---------------
 ! Pre-processing
 !------------------------------------------------
 CALL Init_Elements(Order)
   
 !---------------
 ! Initialization
 !------------------------------------------------
 CALL Initialization(uu, rhs)

 ite = 1; res = 1.d0

 !----------------
 ! Solution update
 !-----------------------------------------------
 DO 

    CALL time_advance(ite, uu, rhs, res)
    
    IF (ite >= ite_max .OR. res < toll_res) EXIT
 
    ite = ite + 1
 
 ENDDO
 
 !----------------
 ! Post-processing
 !----------------------------------------------
 CALL plot_procedure(uu)

 CALL compute_error(uu)

END PROGRAM main

