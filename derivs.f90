MODULE derivs
implicit none

contains

SUBROUTINE deriv_x(f, fx)
use parameters
implicit none
double precision, intent(in) :: f(0:nx,0:nz)
double precision, intent(out) :: fx(0:nx,0:nz)
integer :: j, k

do j = 1, nx1
   fx(j,:) = f(j+1,:) - f(j-1,:)
end do

return
END SUBROUTINE deriv_x

SUBROUTINE deriv_xx(f, fxx)
use parameters
implicit none
double precision, intent(in) :: f(0:nx,0:nz)
double precision, intent(out) :: fxx(0:nx,0:nz)
integer :: j, k

do j = 1, nx1
   fxx(j,:) = f(j+1,:) - 2d0 * f(j,:) + f(j-1,:)
end do

return
END SUBROUTINE deriv_xx

SUBROUTINE deriv_z(f, fz)
use parameters
implicit none
double precision, intent(in) :: f(0:nx,0:nz)
double precision, intent(out) :: fz(0:nx,0:nz)
integer :: j, k

do k = 1, nz1
   fz(1:nx1,k) = f(1:nx1,k+1) - f(1:nx1,k-1)
   fz(1:nx1,0) = -f(1:nx1,2) + 4d0 * f(1:nx1,1) - 3d0 * f(1:nx1,0)
   fz(1:nx1,nz) = f(1:nx1,nz-2) - 4d0 * f(1:nx1,nz1) + 3d0 * &
                   f(1:nx1,nz)
end do

return
END SUBROUTINE deriv_z

SUBROUTINE deriv_zz(f, fzz)
use parameters
implicit none
double precision, intent(in) :: f(0:nx,0:nz)
double precision, intent(out) :: fzz(0:nx,0:nz)
integer :: j, k

do k = 1, nz1
   fzz(1:nx1,k) = f(1:nx1,k+1) - 2d0 * f(1:nx1,k) + f(1:nx1,k-1)
   fzz(1:nx1,0) = -f(1:nx1,3) + 4d0 * f(1:nx1,2) - &
                    5d0 * f(1:nx1,1) + 2d0 * f(1:nx1,0)
   fzz(1:nx1,nz) = -f(1:nx1,nz-3) + 4d0 * f(1:nx1,nz-2) - &
                     5d0 * f(1:nx1,nz1) + 2d0 * f(1:nx1,nz)
end do

return
END SUBROUTINE deriv_zz

END MODULE derivs
