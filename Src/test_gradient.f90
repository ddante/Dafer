module test_gradient

  use Element_Class
  use geometry
  use quadrature_rules
  use init_problem

  implicit none

  real(kind=8), parameter :: Pi = DACOS(-1.d0)

contains

  subroutine accuracy_gradients(uu, D_uu)

    implicit none

    real(kind=8), dimension(:),   intent(in) :: uu
    real(kind=8), dimension(:,:), intent(in) :: D_uu

    real(kind=8), dimension(2) :: D_u_ex
    real(kind=8), dimension(:), allocatable :: err2, d2

    integer :: je, k, i

    real(kind=8) :: err_L2, int_d2

    err_L2 = 0.d0; int_d2 = 0.d0

    do je = 1, N_elements

       allocate( err2(elements(je)%p%N_points), &
                   d2(elements(je)%p%N_points) )

       do k = 1, elements(je)%p%N_points

          D_u_ex = exact_grad(elements(je)%p%coords(1, k), &
                              elements(je)%p%coords(2, k) )
          

          err2(k) = SUM( ( D_uu(:, elements(je)%p%NU(k)) - D_u_ex )**2 )

          d2(k) =  SUM( (  D_u_ex**2 ) )

       enddo

       err_L2 = err_L2 + Int_d(elements(je)%p, err2)
       int_d2 = int_d2 + Int_d(elements(je)%p, d2)

       deallocate(  err2, d2 )

    enddo
    
    err_L2 = sqrt(err_l2)/sqrt(int_d2)

    write(*,22) SIZE(uu), err_l2

22 FORMAT(I6, 1(1x,e24.16))

  end subroutine accuracy_gradients
  


!!$  function exact_grad(x,y) result(gg)
!!$
!!$    implicit none
!!$
!!$    real(kind=8), intent(in) :: x, y
!!$    real(kind=8), dimension(2) :: gg
!!$
!!$    gg(1) = dcos(Pi*x)*dcos(Pi*x)*Pi*dsin(Pi*y)*dcos(Pi*y) - &
!!$            dsin(Pi*x)*dsin(Pi*x)*Pi*dsin(Pi*y)*dcos(Pi*y)
!!$
!!$    gg(2) = dsin(Pi*x)*dcos(Pi*x)*dcos(Pi*y)*dcos(Pi*y)*Pi - &
!!$            dsin(Pi*x)*dcos(Pi*x)*dsin(Pi*y)*dsin(Pi*y)*Pi
!!$
!!$  end function exact_grad


  
 function exact_grad(x,y) result(gg)

    implicit none

    real(kind=8), intent(in) :: x, y
    real(kind=8), dimension(2) :: gg

    gg(1) = 2.d0*dsin(2*Pi*x)*Pi*&
            exp(0.5d0*y* (1.d0 - dsqrt(1.d0+16.d0*(Pi**2)*visc**2) )/visc)

    gg(2) = -dcos(2.d0*Pi*x) * (1.d0 - dsqrt(1.d0 + 16.d0*(Pi**2)*visc**2))/visc * &
             exp(y*(1.d0-dsqrt(1.d0 + 16.d0 * (Pi**2)* visc**2))/visc/2.d0)/2.d0


  end function exact_grad
  

end module test_gradient

