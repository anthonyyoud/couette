MODULE matrices
implicit none

contains

SUBROUTINE matrix_setup(ux, uz, zx, zz)
use parameters
use io
use ic_bc
implicit none
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
   uz%di(0) = 1d0 + rzz + (rz * tau / (1d0 - tau))
   uz%di(nz) = 1d0 + rzz + (rz * tau / (1d0 - tau))
   uz%lo(nz) = -rzz
   uz%up(0) = -rzz
end if

zx%di(:) = 1d0 + rxx + (((1d0 - eta)**2) * dt) / (2d0 * s(1:nx-1)**2)
zx%lo(:) = -0.5d0 * rxx + ((1d0 - eta) * rx) / (4d0 * s(2:nx-1))
zx%up(:) = -0.5d0 * rxx - ((1d0 - eta) * rx) / (4d0 * s(1:nx-2))

zz%di(:) = 1d0 + rzz
zz%lo(:) = -0.5d0 * rzz
zz%up(:) = -0.5d0 * rzz

return
END SUBROUTINE matrix_setup

SUBROUTINE psi_mat_setup(p_mat, IPIV)
use parameters
use ic_bc
implicit none

double precision :: alp(0:nx), gam(0:nx)
double precision :: beta, delta
double precision, intent(out) :: p_mat(2*nx1+nx1+1,nx1*nz1)
integer :: j, k, info
integer, intent(out) :: IPIV(nx1*nz1)

alp(:) = dz2 + 0.5d0 * delx * dz2 * (1d0 - eta) / s(:)
gam(:) = dz2 - 0.5d0 * delx * dz2 * (1d0 - eta) / s(:)
beta = -2d0 * (dz2 + dx2)
delta = dx2

do j = 1, nx1*nz1
   p_mat(2*nx1+1,j) = beta
end do

do j = 1, nx1*nz1-1
   p_mat(2*nx1,j+1) = gam(mod(j, nx1))
end do

do j = nx1, nx1*nz1-nx1, nx1
   p_mat(2*nx1,j+1) = 0d0
end do

do j = 2, nx1*nz1
   p_mat(2*nx1+2,j-1) = alp(mod(j-1, nx1) + 1)
end do

do j = nx, nx1*nz1-nx1+1, nx1
   p_mat(2*nx1+2,j-1) = 0d0
end do

do j = 1, nx1*nz1-nx1
   p_mat(2*nx1+1-nx1,j+nx1) = delta
end do

do j = nx, nx1*nz1
   p_mat(2*nx1+1+nx1,j-nx1) = delta
end do

call DGBTRF(nx1*nz1, nx1*nz1, nx1, nx1, p_mat, 2*nx1+nx1+1, IPIV, info)

!open (60, file = 'p_mat_mat.dat')
!write(60, '(741e17.9)') ((p_mat(j,k), k = 1, nx1*nz1), j = 1, 2*nx1+nx1+1)
!close (60)

return
END SUBROUTINE psi_mat_setup

SUBROUTINE b_mat_setup(b_mat, IPIV)
use parameters
use ic_bc
implicit none

double precision :: alp(0:nx), beta(0:nx), gam(0:nx), delta
double precision, intent(out) :: b_mat(2*nxp1+nxp1+1,0:nxp1*nz1-1)
integer :: j, k, info
integer, intent(out) :: IPIV(nxp1*nz1)

alp(:) = dz2 - 0.5d0 * delx * dz2 * (1d0 - eta) / s(:)
beta(:) = -2d0 * (dz2 + dx2) - dx2 * dz2 * (1d0 - eta)**2 / s(:)**2
gam(:) = dz2 + 0.5d0 * delx * dz2 * (1d0 - eta) / s(:)
delta = dx2

do j = 0, nxp1*nz1-1
   b_mat(2*nx+3,j) = beta(mod(j, nxp1))
end do

do j = 0, nxp1*nz1-2
   b_mat(2*nx+2,j+1) = gam(mod(j, nxp1))
end do

do j = 1, nxp1*nz1-1
   b_mat(2*nx+4,j-1) = alp(mod(j, nxp1))
end do

do j = 0, nxp1*nz1-nxp1-1
   b_mat(nx+2,j+nxp1) = delta
end do

do j = nxp1, nxp1*nz1-1
   b_mat(3*nx+4,j-nxp1) = delta
end do

do j = nx, nxp1*nz1-nxp1-1, nxp1
   b_mat(2*nx+2,j+1) = 0d0
end do

do j = nxp1, nxp1*nz1-nxp1, nxp1
   b_mat(2*nx+4,j-1) = 0d0
end do

do j = 0, nxp1*nz1-nxp1, nxp1
   b_mat(2*nx+3,j) = (2d0 * alp(0) * delx * (1d0 - eta) / s(0)) - beta(0)
end do

do j = nx, nxp1*nz1-1, nxp1
   b_mat(2*nx+3,j) = (-2d0 * gam(nx) * delx * (1d0 - eta) / s(nx)) - beta(nx)
end do

do j = 0, nxp1*nz1-nxp1, nxp1
   b_mat(2*nx+2,j+1) = alp(0) + gam(0)
end do

do j = nx, nxp1*nz1-1, nxp1
   b_mat(2*nx+4,j-1) = alp(nx) + gam(nx)
end do

call DGBTRF(nxp1*nz1, nxp1*nz1, nxp1, nxp1, b_mat, 2*nxp1+nxp1+1, IPIV, info)

return
END SUBROUTINE b_mat_setup

SUBROUTINE j_mat_setup(j_mat, IPIV)
use parameters
use ic_bc
implicit none

double precision :: alp(0:nx), beta(0:nx), gam(0:nx), delta
double precision, intent(out) :: j_mat(2*nx1+nx1+1,nx1*nzp1)
integer :: j, k, info
integer, intent(out) :: IPIV(nx1*nzp1)

alp(:) = dz2 + 0.5d0 * delx * dz2 * (1d0 - eta) / s(:)
gam(:) = dz2 - 0.5d0 * delx * dz2 * (1d0 - eta) / s(:)
beta = -2d0 * (dz2 + dx2)
delta = dx2

do j = 1, nx1*nzp1
   j_mat(2*nx1+1,j) = beta(mod(j-1, nx-1)+1)
end do

do j = 1, nx1*nzp1-1
   j_mat(2*nx1,j+1) = gam(mod(j, nx1))
end do

do j = nx1, nx1*nzp1-nx1, nx1
   j_mat(2*nx1,j+1) = 0d0
end do

do j = 2, nx1*nzp1
   j_mat(2*nx1+2,j-1) = alp(mod(j-1, nx1) + 1)
end do

do j = nx, nx1*nzp1-nx1+1, nx1
   j_mat(2*nx1+2,j-1) = 0d0
end do

do j = 1, nx1*nzp1-nx1
   j_mat(2*nx1+1-nx1,j+nx1) = delta
end do

do j = nx, nx1*nzp1
   j_mat(2*nx1+1+nx1,j-nx1) = delta
end do

do j = 1, nx1
   j_mat(2*nx1+1-nx1,j+nx1) = 2d0 * delta
end do

do j = nx1*nzp1-nx1+1, nx1*nzp1
   j_mat(2*nx1+1+nx1,j-nx1) = 2d0 * delta
end do

call DGBTRF(nx1*nzp1, nx1*nzp1, nx1, nx1, j_mat, 2*nx1+nx1+1, IPIV, info)

return
END SUBROUTINE j_mat_setup

END MODULE matrices
