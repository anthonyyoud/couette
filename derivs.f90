MODULE derivs
!use parameters
implicit none

contains

SUBROUTINE deriv_x(f, fx, f0x, f1x)
use parameters
implicit none
double precision, intent(in) :: f(0:nx,0:nz)
double precision, intent(out) :: fx(0:nx,0:nz)
double precision, intent(out), optional :: f0x(0:nx,0:nz), &
                                           f1x(0:nx,0:nz)
integer :: j, k

do j = 1, nx-1
   fx(j,1:nz-1) = f(j+1,1:nz-1) - f(j-1,1:nz-1)
end do

if (present(f0x)) then
   do j = 1, nx-1
      f0x(j,0) = f(j+1,0) - f(j-1,0)
   end do
end if

if (present(f1x)) then
   do j = 1, nx-1
      f1x(j,nz) = f(j+1,nz) - f(j-1,nz)
   end do
end if

return
END SUBROUTINE deriv_x

SUBROUTINE deriv_xx(f, fxx, f0xx, f1xx)
use parameters
implicit none
double precision, intent(in) :: f(0:nx,0:nz)
double precision, intent(out) :: fxx(0:nx,0:nz)
double precision, intent(out), optional :: f0xx(0:nx,0:nz), &
                                           f1xx(0:nx,0:nz)
integer :: j, k

do j = 1, nx-1
   fxx(j,1:nz-1) = f(j+1,1:nz-1) - 2d0 * f(j,1:nz-1) + f(j-1,1:nz-1)
end do

if (present(f0xx)) then
   do j = 1, nx-1
      f0xx(j,0) = f(j+1,0) - 2d0 * f(j,0) + f(j-1,0)
   end do
end if

if (present(f1xx)) then
   do j = 1, nx-1
      f1xx(j,nz) = f(j+1,nz) - 2d0 * f(j,nz) + f(j-1,nz)
   end do
end if


return
END SUBROUTINE deriv_xx

SUBROUTINE deriv_z(f, fz, f0z, f1z)
use parameters
implicit none
double precision, intent(in) :: f(0:nx,0:nz)
double precision, intent(out) :: fz(0:nx,0:nz)
double precision, intent(out), optional :: f0z(0:nx,0:nz), &
                                           f1z(0:nx,0:nz)
integer :: j, k

do k = 1, nz-1
   fz(1:nx-1,k) = f(1:nx-1,k+1) - f(1:nx-1,k-1)
end do

if (present(f0z)) then
   f0z(1:nx-1,0) = -f(1:nx1,2) + 4d0 * f(1:nx1,1) - 3d0 * f(1:nx1,0)
end if

if (present(f1z)) then
   f1z(1:nx-1,nz) = f(1:nx1,nz-2) - 4d0 * f(1:nx1,nz1) + 3d0 * f(1:nx1,nz)
end if

return
END SUBROUTINE deriv_z

SUBROUTINE deriv_zz(f, fzz, f0zz, f1zz)
use parameters
implicit none
double precision, intent(in) :: f(0:nx,0:nz)
double precision, intent(out) :: fzz(0:nx,0:nz)
double precision, intent(out), optional :: f0zz(0:nx,0:nz), &
                                           f1zz(0:nx,0:nz)
integer :: j, k

do k = 1, nz-1
   fzz(1:nx-1,k) = f(1:nx-1,k+1) - 2d0 * f(1:nx-1,k) + f(1:nx-1,k-1)
end do

if (present(f0zz)) then
   f0zz(1:nx-1,0) = -f(1:nx1,3) + 4d0 * f(1:nx1,2) - &
                     5d0 * f(1:nx1,1) + 2d0 * f(1:nx1,0)
end if

if (present(f1zz)) then
   f1zz(1:nx-1,nz) = -f(1:nx1,nz-3) + 4d0 * f(1:nx1,nz-2) - &
                      5d0 * f(1:nx1,nz1) + 2d0 * f(1:nx1,nz)
end if

return
END SUBROUTINE deriv_zz

END MODULE derivs
