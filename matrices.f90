MODULE matrices
!use parameters
implicit none

contains

SUBROUTINE matrix_setup(ux, uz, zx, zz, s)
use parameters
use io
implicit none
double precision, intent(in) :: s(0:nx)
type (mat_comp), intent(out) :: ux, zx
type (uz_mat_comp), intent(out) :: uz
type (zz_mat_comp), intent(out) :: zz

ux%di(:) = 1d0 + rxx + (((1d0 - eta)**2) * dt) / (2d0 * s(1:nx-1)**2)
ux%lo(:) = -0.5d0 * rxx + ((1d0 - eta) * rx) / (4d0 * s(2:nx-1))
ux%up(:) = -0.5d0 * rxx - ((1d0 - eta) * rx) / (4d0 * s(1:nx-2))

uz%di(:) = 1d0 + rzz
uz%lo(:) = -0.5d0 * rzz
uz%up(:) = -0.5d0 * rzz

if (tau /= 1) then
   uz%lo(nz) = -rzz
   uz%up(0) = -rzz
end if

zx%di(:) = 1d0 + rxx
zx%lo(:) = -0.5d0 * rxx + (3d0 * (1d0 - eta) * rx) / (4d0 * s(2:nx-1))
zx%up(:) = -0.5d0 * rxx - (3d0 * (1d0 - eta) * rx) / (4d0 * s(1:nx-2))

zz%di(:) = 1d0 + rzz
zz%lo(:) = -0.5d0 * rzz
zz%up(:) = -0.5d0 * rzz

return
END SUBROUTINE matrix_setup

SUBROUTINE ABC_mat_setup(AB, IPIV, s)
use parameters
implicit none

double precision, intent(in) :: s(0:nx)
double precision :: alp(0:nx), gam(0:nx)
double precision :: beta, delta
double precision, intent(out) :: AB(2*nx1+nx1+1,nx1*nz1)
integer :: j, k, info
integer, intent(out) :: IPIV(nx1*nz1)

alp(:) = dz2 + 0.5d0 * delx * dz2 * (1d0 - eta) / s(:)
gam(:) = dz2 - 0.5d0 * delx * dz2 * (1d0 - eta) / s(:)
beta = -2d0 * (dz2 + dx2)
delta = dx2

do j = 1, nx1*nz1
   AB(2*nx1+1,j) = beta
end do

do j = 1, nx1*nz1-1
   AB(2*nx1,j+1) = gam(mod(j, nx1))
end do

do j = nx1, nx1*nz1-nx1, nx1
   AB(2*nx1,j+1) = 0d0
end do

do j = 2, nx1*nz1
   AB(2*nx1+2,j-1) = alp(mod(j-1, nx1) + 1)
end do

do j = nx, nx1*nz1-nx1+1, nx1
   AB(2*nx1+2,j-1) = 0d0
end do

do j = 1, nx1*nz1-nx1
   AB(2*nx1+1-nx1,j+nx1) = delta
end do

do j = nx, nx1*nz1
   AB(2*nx1+1+nx1,j-nx1) = delta
end do

call DGBTRF(nx1*nz1, nx1*nz1, nx1, nx1, AB, 2*nx1+nx1+1, IPIV, info)

!open (60, file = 'AB_mat.dat')
!write(60, '(741e17.9)') ((AB(j,k), k = 1, nx1*nz1), j = 1, 2*nx1+nx1+1)
!close (60)

return
END SUBROUTINE ABC_mat_setup

END MODULE matrices
