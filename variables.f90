MODULE variables
use parameters
implicit none

private
public :: copy_var, vr_vz

type, public :: var
!variables
   double precision :: new(0:nx,0:nz)
   double precision :: old(0:nx,0:nz)
   double precision :: old2(0:nx,0:nz)
   double precision :: inter(0:nx,0:nz)
   double precision :: nlin_new(0:nx,0:nz)
   double precision :: nlin_old(0:nx,0:nz)
end type var

type, public :: deriv
!derivatives
   double precision :: x(0:nx,0:nz)
   double precision :: xx(0:nx,0:nz)
   double precision :: z(0:nx,0:nz)
   double precision :: zz(0:nx,0:nz)
   double precision :: zx(0:nx,0:nz)
   double precision :: zxx(0:nx,0:nz)
   double precision :: zzz(0:nx,0:nz)
end type deriv

type, public :: mat_comp
!matrix components
   double precision :: lo(2:nx1)
   double precision :: di(nx1)
   double precision :: up(nx-2)
end type mat_comp

type, public :: uz_mat_comp
!matrix components for u in the z-direction
   double precision :: lo(nz)
   double precision :: di(0:nz)
   double precision :: up(0:nz1)
end type uz_mat_comp

type, public :: zz_mat_comp
!matrix components for Z in the z-direction
   double precision :: lo(2:nz1)
   double precision :: di(nz1)
   double precision :: up(nz-2)
end type zz_mat_comp

type (var), public, save :: ut, zt, psi, bt, jt
double precision, public, save :: vr(0:nx,0:nz), vz(0:nx,0:nz), &
                          vrold(0:nx,0:nz) = 0d0, vzold(0:nx,0:nz) = 0d0

contains

SUBROUTINE copy_var(var_out, var_in)
!Copy a variable
use parameters
implicit none

double precision, intent(in) :: var_in(0:nx,0:nz)
double precision, intent(out) :: var_out(0:nx,0:nz)

var_out = var_in

return
END SUBROUTINE copy_var

SUBROUTINE vr_vz(p, vr, vz)
!Calculate radial and axial velocity components from the stream-function
use parameters
use ic_bc, only : s
implicit none
double precision, intent(in) :: p(0:nx,0:nz)
double precision, intent(out) :: vr(0:nx,0:nz), vz(0:nx,0:nz)
integer :: j, k

do k = 1, nz1
   vr(:,k) = (-1d0 / (2d0 * s(:) * delz)) * (p(:,k+1) - p(:,k-1))
end do

if (tau /= 1) then
   vr(:,0) = (-1d0 / (2d0 * s(:) * delz)) * &
             (-3d0 * p(:,0) + 4d0 * p(:,1) - p(:,2))
   vr(:,nz) = (-1d0 / (2d0 * s(:) * delz)) * &
              (3d0 * p(:,nz) - 4d0 * p(:,nz1) + p(:,nz-2))
else
   vr(:,0) = 0d0
   vr(:,nz) = 0d0
end if

do j = 1, nx1
   vz(j,:) = (1d0 / (2d0 * s(j) * delx)) * (p(j+1,:) - p(j-1,:))
end do

vz(0,:) = 0d0 !(1d0 / (2d0 * s(0) * delx)) * &
              !(-3d0 * p(0,k) + 4d0 * p(1,k) - p(2,k))
vz(nx,:) = 0d0 !(1d0 / (2d0 * s(nx) * delx)) * &
              !(3d0 * p(nx,k) - 4d0 * p(nx1,k) + p(nx-2,k))

return
END SUBROUTINE vr_vz

END MODULE variables
