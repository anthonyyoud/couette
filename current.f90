MODULE current 
implicit none

private
public :: j_poisson, fin_j_poisson

contains

SUBROUTINE j_poisson(p_mat, jn, j_mat, desc_j, af)
!Solve Poisson equation for the azimuthal current when tau/=1
use parameters
use variables
use ic_bc, only : j_BCS, s
use derivs
implicit none

integer (i1), intent(in)  :: desc_j(7)
real (r2),    intent(in)  :: af(laf), p_mat(0:nx,0:nz), j_mat(j_M,j_N)
real (r2),    intent(out) :: jn(0:nx,0:nz)
real (r2)                 :: p_vec(nb), work(lwork_sol)
integer (i1)              :: h, i, j, k, l, info, cpcol, desc_rp(7)
type (deriv)              :: dp

desc_rp(1) = 502
desc_rp(2) = ictxt
desc_rp(3) = nx1*nzp1   !ScaLAPACK descriptor vector for RHS
desc_rp(4) = nb
desc_rp(5) = 0
desc_rp(6) = nb

if (mycol == 0) then
   call deriv_z(p_mat, dp%z)
   call deriv_x(dp%z, dp%zx)   !get derivatives for RHS
   call deriv_xx(dp%z, dp%zxx)
   call deriv_zz(dp%z, dp%zzz)
end if

call SLTIMER(9)
if (npcol > 1) then
   call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zx, nxp1, 0, 0)   !broadcast RHS
   call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zxx, nxp1, 0, 0)  !to all
   call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zzz, nxp1, 0, 0)  !processes
end if
call SLTIMER(9)

cpcol = 0   !initialise current process column

!distribute the RHS as a vector over the process grid, in block column format
do j = 1, nx1*nzp1, nb
   if (mycol == cpcol) then  !make sure each process gets right pieces
      do k = 1, min(nb, nx1*nzp1-j+1)
         i = k + j - 1
         h = modulo(i-1, nx1) + 1
         l = (i-1)/nx1
         p_vec(k) = dx2 * dz2 * (0.5d0 * dp%zzz(h,l) / (s(h) * delz**3) + &
                    0.5d0 * dp%zxx(h,l) / (s(h) * dx2 * delz) - &
                    0.25d0 * (1d0 - eta) * dp%zx(h,l) / &
                    (s(h)**2 * delx * delz))  !transform RHS matrix into vector
      end do
   end if
   if (cpcol == npcol) exit  !if last process then exit
   cpcol = cpcol + 1  !next process column
end do

!Solve Poisson equation using factorised matrix from PDDBTRF
call PDDBTRS('N', nx1*nzp1, nx1, nx1, 1, j_mat, 1, desc_j, p_vec, 1, &
              desc_rp, af, laf, work, lwork_sol, info)
if (info /= 0) print*, 'j_infinite_PDDBTRS ', info

cpcol = 0   !reset current process column

jn = 0d0   !set matrix to zero on all processes
do j = 1, nx1*nzp1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nx1*nzp1-j+1)   !transform distributed RHS vector
         i = k + j - 1                  !into distributed matrix
         h = modulo(i-1, nx1) + 1
         l = (i-1)/nx1
         jn(h,l) = p_vec(k)
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call SLTIMER(7)                                            !collect distributed
if (npcol > 1) then
   call DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, jn, nxp1, 0, 0)  !matrix onto
end if
call SLTIMER(7)                                            !master process

if (mycol == 0) then
   call j_BCS(jn)   !update boundary conditions
end if

return
END SUBROUTINE j_poisson

SUBROUTINE fin_j_poisson(p_mat, jn, j_mat, desc_j, af)
!Solve Poisson equation for the azimuthal current when tau=1.
!Algorithm as above but indices change to reflect different dimensions.
use parameters
use variables
use ic_bc, only : j_BCS, s
use derivs
implicit none

integer (i1), intent(in)  :: desc_j(7)
real (r2),    intent(in)  :: af(laf), p_mat(0:nx,0:nz), j_mat(j_M,j_N)
real (r2),    intent(out) :: jn(0:nx,0:nz)
real (r2)                 :: p_vec(nb)
integer (i1)              :: h, i, j, k, l, info, cpcol, desc_rp(7)
real (r2)                 :: work(lwork_sol)
type (deriv)              :: dp

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
if (npcol > 1) then
   call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zx, nxp1, 0, 0)
   call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zxx, nxp1, 0, 0)
   call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zzz, nxp1, 0, 0)
end if
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
if (npcol > 1) then
   call DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, jn, nxp1, 0, 0)
end if
call SLTIMER(7)

if (mycol == 0) then
   call j_BCS(jn)
end if

return
END SUBROUTINE fin_j_poisson

END MODULE current
