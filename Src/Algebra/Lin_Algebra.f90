MODULE Lin_Algebra

  IMPLICIT NONE

CONTAINS

  !====================================
  FUNCTION Cross_Product(AA) RESULT(Vn)
  !====================================

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: AA

    REAL(KIND=8), DIMENSION(SIZE(AA,2)) :: Vn
    !----------------------------------------------

    REAL(KIND=8), DIMENSION(3)   :: a, b, c
    REAL(KIND=8), DIMENSION(3,3) :: a_cross
    !----------------------------------------------

    IF( (SIZE(AA,2) /= SIZE(AA,1)+1) .AND. &
        (SIZE(AA,2) /= 2 .OR. SIZE(AA,2) /= 3) ) THEN
        WRITE(*,*) 'ERROR: inconsistent sizes for the cross product'   
        STOP
    ENDIF

    a = 0.d0
    a(1:SIZE(AA,2)) = AA(1, :)

    b = 0.d0
    IF(SIZE(AA,1) == 1) THEN
       b(3) = 1.d0
    ELSE
       b(:) = AA(2, :)
    ENDIF
        
    a_cross(1, :) = (/  0.d0, -a(3),  a(2) /)
    a_cross(2, :) = (/  a(3),  0.d0, -a(1) /)
    a_cross(3, :) = (/ -a(2),  a(1),  0.d0 /)

    c = MATMUL(a_cross, b)   

    Vn = c(1:SIZE(AA,2))
    
  END FUNCTION Cross_Product
  !=========================

END MODULE Lin_Algebra
