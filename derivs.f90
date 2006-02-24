module derivs
  !Routines which calculate finite-difference derivatives.  These only
  !calculate the numerator of the derivatives; the denominator dx, dz etc. is
  !included in other routines whenever the derivative is used.
  implicit none

  private
  public :: deriv_x, deriv_z, deriv_xx, deriv_zz

  contains

  subroutine deriv_x(f, fx)
    !First x derivative
    use parameters
    implicit none

    real    (r2), intent(in)  :: f(0:nx,0:nz)
    real    (r2), intent(out) :: fx(0:nx,0:nz)
    integer (i1)              :: j, k

    do j = 1, nx1
      fx(j,:) = f(j+1,:) - f(j-1,:)
    end do

    return
  end subroutine deriv_x

  subroutine deriv_xx(f, fxx)
    !Second x derivative
    use parameters
    implicit none

    real    (r2), intent(in)  :: f(0:nx,0:nz)
    real    (r2), intent(out) :: fxx(0:nx,0:nz)
    integer (i1)              :: j, k

    do j = 1, nx1
      fxx(j,:) = f(j+1,:) - 2.0_r2 * f(j,:) + f(j-1,:)
    end do

    return
  end subroutine deriv_xx

  subroutine deriv_z(f, fz)
    !First z derivative
    use parameters
    implicit none

    real    (r2), intent(in)  :: f(0:nx,0:nz)
    real    (r2), intent(out) :: fz(0:nx,0:nz)
    integer (i1)              :: j, k

    do k = 1, nz1
      fz(1:nx1,k) = f(1:nx1,k+1) - f(1:nx1,k-1)
      fz(1:nx1,0) = -f(1:nx1,2) + 4.0_r2 * f(1:nx1,1) - 3.0_r2 * f(1:nx1,0)
      fz(1:nx1,nz) = f(1:nx1,nz-2) - 4.0_r2 * f(1:nx1,nz1) + 3.0_r2 * &
                     f(1:nx1,nz)   !forward/backward difference at boundaries
    end do

    return
  end subroutine deriv_z

  subroutine deriv_zz(f, fzz)
    !Second z derivative
    use parameters
    implicit none

    real    (r2), intent(in)  :: f(0:nx,0:nz)
    real    (r2), intent(out) :: fzz(0:nx,0:nz)
    integer (i1)              :: j, k

    do k = 1, nz1
      fzz(1:nx1,k) = f(1:nx1,k+1) - 2.0_r2 * f(1:nx1,k) + f(1:nx1,k-1)
      fzz(1:nx1,0) = -f(1:nx1,3) + 4.0_r2 * f(1:nx1,2) - &
                     5.0_r2 * f(1:nx1,1) + 2.0_r2 * f(1:nx1,0) !forward/backward
      fzz(1:nx1,nz) = -f(1:nx1,nz-3) + 4.0_r2 * f(1:nx1,nz-2) - & !difference at
                      5.0_r2 * f(1:nx1,nz1) + 2.0_r2 * f(1:nx1,nz)  !boundaries
    end do

    return
  end subroutine deriv_zz

end module derivs
