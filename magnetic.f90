MODULE magnetic 
!Algorithms for the solution of the magnetic Poisson equation are as
!for the current in current.f90
implicit none

private
public :: b_poisson, fin_b_poisson

contains

SUBROUTINE b_poisson(u_mat, bn, b_mat, desc_b, af)
!Solve Poisson equation for azimuthal magnetic field when tau=0
use parameters
use ic_bc, only : b_BCS
implicit none

integer :: desc_u(7)
integer, intent(in) :: desc_b(7)
double precision, intent(in) :: af(b_laf)
double precision, intent(in) :: u_mat(0:nx,0:nz)
double precision, intent(in) :: b_mat(b_M,b_N)
double precision, intent(out) :: bn(0:nx,0:nz)
double precision :: u_vec(nb)
integer :: i, j, k, info, cpcol
double precision :: work(lwork_b_sol)

desc_u(1) = 502
desc_u(2) = ictxt
desc_u(3) = nxp1*nz1
desc_u(4) = nb
desc_u(5) = 0
desc_u(6) = nb

cpcol = 0

do j = 0, nxp1*nz1-1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nxp1*nz1-j)
         i = k + j - 1
         u_vec(k) = 0.5d0 * dx2 * delz * &
                              (u_mat(modulo(i, nxp1), i/nxp1) - &
                               u_mat(modulo(i, nxp1), i/nxp1 + 2))
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call PDDBTRS('N', nxp1*nz1, nxp1, nxp1, 1, b_mat, 1, desc_b, u_vec, 1, &
              desc_u, af, b_laf, work, lwork_b_sol, info)
if (info /= 0) print*, 'b_infinite_PDDBTRS ', info

cpcol = 0

bn = 0d0
do j = 0, nxp1*nz1-1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nxp1*nz1-j)
         i = k + j - 1
         bn(modulo(i, nxp1), i/nxp1 + 1) = u_vec(k)
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call SLTIMER(8)
if (npcol > 1) then
   call DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, bn, nxp1, 0, 0)
end if
call SLTIMER(8)

if (mycol == 0) then
   call b_BCS(bn)
end if

return
END SUBROUTINE b_poisson

SUBROUTINE fin_b_poisson(u_mat, bn, b_mat, desc_b, af)
!Solve Poisson equation for azimuthal magnetic field when tau/=0
use parameters
use derivs, only : deriv_z
implicit none

integer :: desc_u(7)
integer, intent(in) :: desc_b(7)
double precision, intent(in) :: af(b_laf)
double precision, intent(in) :: u_mat(0:nx,0:nz)
double precision, intent(in) :: b_mat(b_M,b_N)
double precision, intent(out) :: bn(0:nx,0:nz)
double precision :: u_vec(nb), u_mat_z(0:nx,0:nz)
integer :: i, j, k, info, cpcol
double precision :: work(lwork_b_sol)

if (mycol == 0) then
   call deriv_z(u_mat, u_mat_z)
end if

if (npcol > 1) then
   call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, u_mat_z, nxp1, 0, 0)
end if

desc_u(1) = 502
desc_u(2) = ictxt
desc_u(3) = nxp1*nzp1
desc_u(4) = nb
desc_u(5) = 0
desc_u(6) = nb

cpcol = 0

do j = 0, nxp1*nzp1-1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nxp1*nzp1-j)
         i = k + j - 1
         u_vec(k) = -0.5d0 * dx2 * delz * &
                     u_mat_z(modulo(i, nxp1), i/nxp1)
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call PDDBTRS('N', nxp1*nzp1, nxp1, nxp1, 1, b_mat, 1, desc_b, u_vec, 1, &
              desc_u, af, b_laf, work, lwork_b_sol, info)
if (info /= 0) print*, 'b_finite_PDDBTRS ', info

cpcol = 0

bn = 0d0
do j = 0, nxp1*nzp1-1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nxp1*nzp1-j)
         i = k + j - 1
         bn(modulo(i, nxp1), i/nxp1) = u_vec(k)
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call SLTIMER(8)
if (npcol > 1) then
   call DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, bn, nxp1, 0, 0)
end if
call SLTIMER(8)

return
END SUBROUTINE fin_b_poisson

END MODULE magnetic
