MODULE matrices
implicit none

private
public :: matrix_setup, psi_mat_setup, b_mat_setup, fin_b_mat_setup, &
          j_mat_setup, fin_j_mat_setup

contains

SUBROUTINE matrix_setup(ux, uz, zx, zz)
!Set up the upper, lower and diagonal parts of the tridiagonal matrices
!for the solution of the azimuthal velocity and vorticity equations
use parameters
use variables
use ic_bc, only : s
implicit none

type (mat_comp),    intent(out) :: ux, zx
type (uz_mat_comp), intent(out) :: uz
type (zz_mat_comp), intent(out) :: zz

ux%di(:) = 1.0_r2 + rxx + (((1.0_r2 - eta)**2) * dt) / (2.0_r2 * s(1:nx-1)**2)
ux%lo(:) = -0.5_r2 * rxx + ((1.0_r2 - eta) * rx) / (4.0_r2 * s(2:nx-1))
ux%up(:) = -0.5_r2 * rxx - ((1.0_r2 - eta) * rx) / (4.0_r2 * s(1:nx-2))

uz%di(:) = 1.0_r2 + rzz
uz%lo(:) = -0.5_r2 * rzz
uz%up(:) = -0.5_r2 * rzz

if (tau /= 1) then
   uz%di(0) = 1.0_r2 + rzz + (rz * tau / (1.0_r2 - tau))   !extra entries due
   uz%di(nz) = 1.0_r2 + rzz + (rz * tau / (1.0_r2 - tau))  !to Neumann BCS at ends
   uz%lo(nz) = -rzz
   uz%up(0) = -rzz
end if

zx%di(:) = 1.0_r2 + rxx + (((1.0_r2 - eta)**2) * dt) / (2.0_r2 * s(1:nx-1)**2)
zx%lo(:) = -0.5_r2 * rxx + ((1.0_r2 - eta) * rx) / (4.0_r2 * s(2:nx-1))
zx%up(:) = -0.5_r2 * rxx - ((1.0_r2 - eta) * rx) / (4.0_r2 * s(1:nx-2))

zz%di(:) = 1.0_r2 + rzz
zz%lo(:) = -0.5_r2 * rzz
zz%up(:) = -0.5_r2 * rzz

return
END SUBROUTINE matrix_setup

SUBROUTINE psi_mat_setup(p_mat, desc_p, af)
!Setup of LHS matrix in solution of stream-function Poisson equation
use parameters
use ic_bc, only : s
implicit none

integer (i1), intent(out) :: desc_p(7)
real (r2),    intent(out) :: p_mat(p_M,p_N), af(laf)
real (r2)                 :: alp(0:nx), gam(0:nx), beta, delta, work(lwork_fac)
integer (i1)              :: i, j, k, info, cpcol

alp(:) = dz2 + 0.5_r2 * delx * dz2 * (1.0_r2 - eta) / s(:)   !coefficients
gam(:) = dz2 - 0.5_r2 * delx * dz2 * (1.0_r2 - eta) / s(:)   !in matrix
beta = -2.0_r2 * (dz2 + dx2)
delta = dx2

desc_p(1) = 501
desc_p(2) = ictxt
desc_p(3) = nx1*nz1   !ScaLAPACK descriptor vector for LHS matrix
desc_p(4) = nb
desc_p(5) = 0
desc_p(6) = 2*nx1+1

cpcol = 0

!distribute the matrix over the process grid using the block-column
!distribution for banded matrices
do j = 1, nx1*nz1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nx1*nz1-j+1)
         i = k + j - 1
         p_mat(nx,k) = beta   !diagonal
         if (i /= 1) then
            p_mat(nx1,k) = gam(modulo(i-1, nx1))   !upper-diagonal
            if(modulo(i-1, nx1) == 0) then
               p_mat(nx1,k) = 0.0_r2   !upper-diagonal, BCS
            end if
         end if
         if (i /= nx1*nz1) then
            p_mat(nxp1,k) = alp(modulo(i, nx1) + 1) !lower-diagonal
            if (modulo(i, nx1) == 0) then
               p_mat(nxp1,k) = 0.0_r2   !lower-diagonal, BCS
            end if
         end if
         if (i > nx1) then
            p_mat(1,k) = delta  !upper band
         end if
         if (i <= nx1*nz1-nx1) then
            p_mat(2*nx1+1,k) = delta   !lower band
        end if
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

!LU factorisation of the matrix for use in PDDBTRS to solve
call PDDBTRF(nx1*nz1, nx1, nx1, p_mat, 1, desc_p, af, laf, &
             work, lwork_fac, info)
if (info /= 0) print*, 'psi_PDDBTRF ', info

return
END SUBROUTINE psi_mat_setup

SUBROUTINE b_mat_setup(b_mat, desc_b, af)
!Setup of LHS matrix in solution of magnetic Poisson equation.
!Algorithm as for stream-function above.
use parameters
use ic_bc, only : s
implicit none

integer (i1), intent(out) :: desc_b(7)
real (r2),    intent(out) :: b_mat(b_M,b_N), af(b_laf)
real (r2)                 :: alp(0:nx), beta(0:nx), gam(0:nx), delta, &
                             work(lwork_b_fac)
integer (i1)              :: i, j, k, info, cpcol

alp(:) = dz2 - 0.5_r2 * delx * dz2 * (1.0_r2 - eta) / s(:)
beta(:) = -2.0_r2 * (dz2 + dx2) - dx2 * dz2 * (1.0_r2 - eta)**2 / s(:)**2
gam(:) = dz2 + 0.5_r2 * delx * dz2 * (1.0_r2 - eta) / s(:)
delta = dx2

desc_b(1) = 501
desc_b(2) = ictxt
desc_b(3) = nxp1*nz1
desc_b(4) = nb
desc_b(5) = 0
desc_b(6) = 2*nxp1+1

cpcol = 0

do j = 0, nxp1*nz1-1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nxp1*nz1-j)
         i = k + j - 1
         b_mat(nx+2,k) = beta(modulo(i, nxp1))   !diagonal
         if (modulo(i, nxp1) == 0) then
            b_mat(nx+2,k) = (2.0_r2 * alp(0) * delx * (1.0_r2 - eta) / &
                          s(0)) + beta(0)   !diagonal, BCS
         end if
         if (modulo(i+1, nxp1) == 0) then
            b_mat(nx+2,k) = (-2.0_r2 * gam(nx) * delx * (1.0_r2 - eta) / &
                          s(nx)) + beta(nx) !diagonal, BCS
         end if
         if (i /= 0) then
            b_mat(nxp1,k) = gam(modulo(i-1, nxp1))   !upper-diagonal
            if(modulo(i, nxp1) == 0) then
               b_mat(nxp1,k) = 0.0_r2   !upper-diagonal, BCS
            end if
            if(modulo(i-1, nxp1) == 0) then
               b_mat(nxp1,k) = alp(0) + gam(0)   !upper-diagonal, BCS
            end if
         end if
         if (i /= nxp1*nz1-1) then
            if (modulo(i+1, nxp1) /= 0) then
               b_mat(nx+3,k) = alp(modulo(i+1,nxp1))   !lower-diagonal
            end if
            if (modulo(i+1, nxp1) == 0) then
               b_mat(nx+3,k) = 0.0_r2   !lower-diagonal, BCS
            end if
            if (modulo(i+2, nxp1) == 0) then
               b_mat(nx+3,k) = alp(nx) + gam(nx)   !lower-diagonal, BCS
            end if
         end if
         if (i > nx) then
            b_mat(1,k) = delta   !upper band
         end if
         if (i < nxp1*nz1-nxp1) then
            b_mat(2*nxp1+1,k) = delta   !lower band
         end if
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call PDDBTRF(nxp1*nz1, nxp1, nxp1, b_mat, 1, desc_b, af, b_laf, &
             work, lwork_b_fac, info)
if (info /= 0) print*, 'b_infinite_PDDBTRF ', info

return
END SUBROUTINE b_mat_setup

SUBROUTINE fin_b_mat_setup(b_mat, desc_b, af)
use parameters
use ic_bc, only : s
implicit none

integer (i1), intent(out) :: desc_b(7)
real (r2),    intent(out) :: b_mat(b_M,b_N), af(b_laf)
real (r2)                 :: alp(0:nx), beta(0:nx), gam(0:nx), delta, &
                             work(lwork_b_fac)
integer (i1)              :: i, j, k, info, cpcol

alp(:) = dz2 - 0.5_r2 * delx * dz2 * (1.0_r2 - eta) / s(:)
beta(:) = -2.0_r2 * (dz2 + dx2) - dx2 * dz2 * (1.0_r2 - eta)**2 / s(:)**2
gam(:) = dz2 + 0.5_r2 * delx * dz2 * (1.0_r2 - eta) / s(:)
delta = dx2

desc_b(1) = 501
desc_b(2) = ictxt
desc_b(3) = nxp1*nzp1
desc_b(4) = nb
desc_b(5) = 0
desc_b(6) = 2*nxp1+1

cpcol = 0

do j = 0, nxp1*nzp1-1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nxp1*nzp1-j)
         i = k + j - 1
         !diagonal
         b_mat(nx+2,k) = beta(modulo(i, nxp1))
         !diagonal, j=0
         if (modulo(i, nxp1) == 0) then
            b_mat(nx+2,k) = (2.0_r2 * alp(0) * delx * (1.0_r2 - eta) / &
                          s(0)) + beta(0)
         end if
         !diagonal, k=0
         if (i < nxp1) then
            b_mat(nx+2,k) = beta(modulo(i, nxp1)) - &
                         2.0_r2 * delta * delz * (1.0_r2 - tau) / tau
         end if
         !diagonal, j=nx
         if (modulo(i+1, nxp1) == 0) then
            b_mat(nx+2,k) = (-2.0_r2 * gam(nx) * delx * (1.0_r2 - eta) / &
                          s(nx)) + beta(nx)
         end if
         !diagonal, k=nz
         if (i > nxp1*nzp1-nxp1) then
           b_mat(nx+2,k) = beta(modulo(i, nxp1)) - &
                         2.0_r2 * delta * delz * (1.0_r2 - tau) / tau
         end if
         !diagonal, j=0, k=0
         if (i == 0) then
            b_mat(nx+2,k) = beta(0) - &
                         2.0_r2 * delta * delz * (1.0_r2 - tau) / tau + &
                        (2.0_r2 * alp(0) * delx * (1.0_r2 - eta) / s(0))
         end if
         !diagonal, j=0, k=nz
         if (i == nxp1*nzp1-nxp1) then
            b_mat(nx+2,k) = beta(0) - &
                         2.0_r2 * delta * delz * (1.0_r2 - tau) / tau + &
                        (2.0_r2 * alp(0) * delx * (1.0_r2 - eta) / s(0))
         end if
         !diagonal, j=nx, k=0
         if (i == nx) then
            b_mat(nx+2,k) = beta(nx) - &
                         2.0_r2 * delta * delz * (1.0_r2 - tau) / tau - &
                        (2.0_r2 * gam(nx) * delx * (1.0_r2 - eta) / s(nx))
         end if
         !diagonal, j=nx, k=nz
         if (i == nxp1*nzp1-1) then
            b_mat(nx+2,k) = beta(nx) - &
                         2.0_r2 * delta * delz * (1.0_r2 - tau) / tau - &
                        (2.0_r2 * gam(nx) * delx * (1.0_r2 - eta) / s(nx))
         end if
         if (i /= 0) then
            !super-diagonal
            b_mat(nxp1,k) = gam(modulo(i-1, nxp1))
            !super-diagonal, j=nx, so j=0 not present
            if(modulo(i, nxp1) == 0) then
               b_mat(nxp1,k) = 0.0_r2
            end if
            !super-diagonal, j=0
            if(modulo(i-1, nxp1) == 0) then
               b_mat(nxp1,k) = alp(0) + gam(0)
            end if
         end if
         if (i /= nxp1*nzp1-1) then
            !sub-diagonal
            if (modulo(i+1, nxp1) /= 0) then
               b_mat(nx+3,k) = alp(modulo(i+1, nxp1))
            end if
            !sub-diagonal, j=0, so j=nx not present
            if (modulo(i+1, nxp1) == 0) then
               b_mat(nx+3,k) = 0.0_r2
            end if
            !sub-diagonal, j=nx
            if (modulo(i+2, nxp1) == 0) then
               b_mat(nx+3,k) = alp(nx) + gam(nx)
            end if
         end if
         !upper branch
         if (i > nx) then
            b_mat(1,k) = delta
            if (i < 2*nxp1) then
               b_mat(1,k) = 2.0_r2 * delta
            end if
         end if
         !lower branch
         if (i < nxp1*nzp1-nxp1) then
            b_mat(2*nxp1+1,k) = delta
            if (i >= nxp1*nzp1-2*nxp1) then
               b_mat(2*nxp1+1,k) = 2.0_r2 * delta
            end if
         end if
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call PDDBTRF(nxp1*nzp1, nxp1, nxp1, b_mat, 1, desc_b, af, b_laf, &
             work, lwork_b_fac, info)
if (info /= 0) print*, 'b_finite_PDDBTRF ', info

return
END SUBROUTINE fin_b_mat_setup

SUBROUTINE j_mat_setup(j_mat, desc_j, af)
use parameters
use ic_bc, only : s
implicit none

integer (i1), intent(out) :: desc_j(7)
real (r2),    intent(out) :: j_mat(j_M,j_N), af(laf)
real (r2)                 :: alp(0:nx), beta(0:nx), gam(0:nx), delta, &
                             work(lwork_fac)
integer (i1)              :: i, j, k, info, cpcol

alp(:) = dz2 - 0.5_r2 * delx * dz2 * (1.0_r2 - eta) / s(:)
beta(:) = -2.0_r2 * (dz2 + dx2) - dx2 * dz2 * (1.0_r2 - eta)**2 / s(:)**2
gam(:) = dz2 + 0.5_r2 * delx * dz2 * (1.0_r2 - eta) / s(:)
delta = dx2

desc_j(1) = 501
desc_j(2) = ictxt
desc_j(3) = nx1*nzp1
desc_j(4) = nb
desc_j(5) = 0
desc_j(6) = 2*nx1+1

cpcol = 0

do j = 1, nx1*nzp1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nx1*nzp1-j+1)
         i = k + j - 1
         j_mat(nx,k) = beta(modulo(i-1, nx1) + 1)   !diagonal
         if (i <= nx1) then
            j_mat(nx,k) = beta(modulo(i-1, nx1) + 1) - &      !extra diagonal
                       2.0_r2 * delta * delz * tau / (1.0_r2 - tau) !due to tau
         end if
         if (i >= nx1*nzp1-nx1+1) then
            j_mat(nx,k) = beta(modulo(i-1, nx1) + 1) - &      !extra diagonal
                       2.0_r2 * delta * delz * tau / (1.0_r2 - tau) !due to tau
         end if
         if (i /= 1) then
            j_mat(nx1,k) = gam(modulo(i-1, nx1))   !upper-diagonal
            if(modulo(i-1, nx1) == 0) then
               j_mat(nx1,k) = 0.0_r2  !upper-diagonal, BCS
            end if
         end if
         if (i /= nx1*nzp1) then
            j_mat(nxp1,k) = alp(modulo(i,nx1) + 1)   !lower-diagonal
            if (modulo(i, nx1) == 0) then
               j_mat(nxp1,k) = 0.0_r2   !lower-diagonal, BCS
            end if
         end if
         if (i > nx1) then
            j_mat(1,k) = delta   !upper band
            if (i <= 2*nx1) then
               j_mat(1,k) = 2.0_r2 * delta   !upper band, BCS
            end if
         end if
         if (i <= nx1*nzp1-nx1) then
            j_mat(2*nx1+1,k) = delta   !lower band
            if (i > nx1*nzp1-2*nx1) then
               j_mat(2*nx1+1,k) = 2.0_r2 * delta   !lower band, BCS
            end if
         end if
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call PDDBTRF(nx1*nzp1, nx1, nx1, j_mat, 1, desc_j, af, laf, &
             work, lwork_fac, info)
if (info /= 0) print*, 'j_infinite_PDDBTRF ', info

return
END SUBROUTINE j_mat_setup

SUBROUTINE fin_j_mat_setup(j_mat, desc_j, af)
use parameters
use ic_bc, only : s
implicit none

integer (i1), intent(out) :: desc_j(7)
real (r2),    intent(out) :: j_mat(j_M,j_N), af(laf)
real (r2)                 :: alp(0:nx), beta(0:nx), gam(0:nx), delta, &
                             work(lwork_fac)
integer (i1)              :: i, j, k, info, cpcol

alp(:) = dz2 - 0.5_r2 * delx * dz2 * (1.0_r2 - eta) / s(:)
beta(:) = -2.0_r2 * (dz2 + dx2) - dx2 * dz2 * (1.0_r2 - eta)**2 / s(:)**2
gam(:) = dz2 + 0.5_r2 * delx * dz2 * (1.0_r2 - eta) / s(:)
delta = dx2

desc_j(1) = 501
desc_j(2) = ictxt
desc_j(3) = nx1*nz1
desc_j(4) = nb
desc_j(5) = 0
desc_j(6) = 2*nx1+1

cpcol = 0

do j = 1, nx1*nz1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nx1*nz1-j+1)
         i = k + j - 1
         j_mat(nx,k) = beta(modulo(i-1, nx1) + 1)   !diagonal
         if (i /= 1) then
            j_mat(nx1,k) = gam(modulo(i-1, nx1))   !upper-diagonal
            if(modulo(i-1, nx1) == 0) then
               j_mat(nx1,k) = 0.0_r2   !upper-diagonal, BCS
            end if
         end if
         if (i /= nx1*nz1) then
            j_mat(nxp1,k) = alp(modulo(i, nx1) + 1)   !lower-diagonal
            if (modulo(i, nx1) == 0) then
               j_mat(nxp1,k) = 0.0_r2   !lower-diagonal, BCS
            end if
         end if
         if (i > nx1) then
            j_mat(1,k) = delta   !upper band
         end if
         if (i <= nx1*nz1-nx1) then
            j_mat(2*nx1+1,k) = delta   !lower band
         end if
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call PDDBTRF(nx1*nz1, nx1, nx1, j_mat, 1, desc_j, af, laf, &
             work, lwork_fac, info)
if (info /= 0) print*, 'j_infinite_PDDBTRF ', info

return
END SUBROUTINE fin_j_mat_setup

END MODULE matrices
