MODULE current 
implicit none

contains

SUBROUTINE j_poisson(p_mat, jn, j_mat, desc_j, af)
use parameters
use variables
use ic_bc
use derivs
implicit none

integer :: desc_rp(7)
integer, intent(in) :: desc_j(7)
double precision, intent(in) :: af(laf)
double precision, intent(in) :: p_mat(0:nx,0:nz)
double precision, intent(in) :: j_mat(j_M,j_N)
double precision, intent(out) :: jn(0:nx,0:nz)
double precision :: p_vec(nb)
integer :: h, i, j, k, l, info, cpcol
double precision :: work(lwork_sol)
type (deriv) :: dp

desc_rp(1) = 502
desc_rp(2) = ictxt
desc_rp(3) = nx1*nzp1
desc_rp(4) = nb
desc_rp(5) = 0
desc_rp(6) = nb

if (mycol == 0) then
   call deriv_z(p_mat, dp%z)
   call deriv_x(dp%z, dp%zx)
   call deriv_xx(dp%z, dp%zxx)
   call deriv_zz(dp%z, dp%zzz)
end if

call SLTIMER(9)
call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zx, nxp1, 0, 0)
call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zxx, nxp1, 0, 0)
call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zzz, nxp1, 0, 0)
call SLTIMER(9)

cpcol = 0

do j = 1, nx1*nzp1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nx1*nzp1-j+1)
         i = k + j - 1
         h = modulo(i-1, nx1) + 1
         l = (i-1)/nx1
         p_vec(k) = dx2 * dz2 * (0.5d0 * dp%zzz(h,l) / (s(h) * delz**3) + &
                    0.5d0 * dp%zxx(h,l) / (s(h) * dx2 * delz) - &
                    0.25d0 * (1d0 - eta) * dp%zx(h,l) / &
                    (s(h)**2 * delx * delz))
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call PDDBTRS('N', nx1*nzp1, nx1, nx1, 1, j_mat, 1, desc_j, p_vec, 1, &
              desc_rp, af, laf, work, lwork_sol, info)
if (info /= 0) print*, 'j_infinite_PDDBTRS ', info

cpcol = 0

jn = 0d0
do j = 1, nx1*nzp1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nx1*nzp1-j+1)
         i = k + j - 1
         h = modulo(i-1, nx1) + 1
         l = (i-1)/nx1
         jn(h,l) = p_vec(k)
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call SLTIMER(7)
call DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, jn, nxp1, 0, 0)
call SLTIMER(7)

if (mycol == 0) then
   call j_BCS(jn)
end if

return
END SUBROUTINE j_poisson

SUBROUTINE fin_j_poisson(p_mat, jn, j_mat, desc_j, af)
use parameters
use variables
use ic_bc
use derivs
implicit none

integer :: desc_rp(7)
integer, intent(in) :: desc_j(7)
double precision, intent(in) :: af(laf)
double precision, intent(in) :: p_mat(0:nx,0:nz)
double precision, intent(in) :: j_mat(j_M,j_N)
double precision, intent(out) :: jn(0:nx,0:nz)
double precision :: p_vec(nb)
integer :: h, i, j, k, l, info, cpcol
double precision :: work(lwork_sol)
type (deriv) :: dp

desc_rp(1) = 502
desc_rp(2) = ictxt
desc_rp(3) = nx1*nz1
desc_rp(4) = nb
desc_rp(5) = 0
desc_rp(6) = nb

if (mycol == 0) then
   call deriv_z(p_mat, dp%z)
   call deriv_x(dp%z, dp%zx)
   call deriv_xx(dp%z, dp%zxx)
   call deriv_zz(dp%z, dp%zzz)
end if

call SLTIMER(9)
call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zx, nxp1, 0, 0)
call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zxx, nxp1, 0, 0)
call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zzz, nxp1, 0, 0)
call SLTIMER(9)

cpcol = 0

do j = 1, nx1*nz1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nx1*nz1-j+1)
         i = k + j - 1
         h = modulo(i-1, nx1) + 1
         l = (i-1)/nx1 + 1
         p_vec(k) = dx2 * dz2 * (0.5d0 * dp%zzz(h,l) / &
                    (s(h) * delz**3) + &
                    0.5d0 * dp%zxx(h,l) / (s(h) * dx2 * delz) - &
                    0.25d0 * (1d0 - eta) * dp%zx(h,l) / &
                    (s(h)**2 * delx * delz))
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call PDDBTRS('N', nx1*nz1, nx1, nx1, 1, j_mat, 1, desc_j, p_vec, 1, &
              desc_rp, af, laf, work, lwork_sol, info)
if (info /= 0) print*, 'j_finite_PDDBTRS ', info

cpcol = 0

jn = 0d0
do j = 1, nx1*nz1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nx1*nz1-j+1)
         i = k + j - 1
         h = modulo(i-1, nx1) + 1
         l = (i-1)/nx1 + 1
         jn(h,l) = p_vec(k)
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call SLTIMER(7)
call DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, jn, nxp1, 0, 0)
call SLTIMER(7)

if (mycol == 0) then
   call j_BCS(jn)
end if

return
END SUBROUTINE fin_j_poisson

END MODULE current
