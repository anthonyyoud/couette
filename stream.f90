MODULE stream
!Algorithms for the solution of the stream-function Poisson equation are as
!for the current in current.f90
implicit none

contains

SUBROUTINE p_poisson(Z_mat, psi, p_mat, desc_p, af)
!Solve Poisson equation for the stream-function, psi for all tau
use parameters
use ic_bc
implicit none

integer :: desc_z(7)
integer, intent(in) :: desc_p(7)
double precision, intent(in) :: af(laf)
double precision, intent(in) :: Z_mat(0:nx,0:nz)
double precision, intent(in) :: p_mat(p_M,p_N)
double precision, intent(out) :: psi(0:nx,0:nz)
double precision :: zvec(nb)
integer :: i, j, k, info, cpcol
double precision :: work(lwork_sol)

desc_z(1) = 502
desc_z(2) = ictxt
desc_z(3) = nx1*nz1
desc_z(4) = nb
desc_z(5) = 0
desc_z(6) = nb

cpcol = 0

do j = 1, nx1*nz1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nx1*nz1-j+1)
         i = k + j - 1
         zvec(k) = -s(modulo(i-1, nx1) + 1) * dx2 * dz2 * &
                    Z_mat(modulo(i-1, nx1) + 1, (i-1)/nx1 + 1)
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call PDDBTRS('N', nx1*nz1, nx1, nx1, 1, p_mat, 1, desc_p, zvec, 1, &
              desc_z, af, laf, work, lwork_sol, info)
if (info /= 0) print*, 'psi_PDDBTRS ', info

cpcol = 0

psi = 0d0
do j = 1, nx1*nz1, nb
   if (mycol == cpcol) then
      do k = 1, min(nb, nx1*nz1-j+1)
         i = k + j - 1
         psi(modulo(i-1, nx1) + 1, (i-1)/nx1 + 1) = zvec(k)
      end do
   end if
   if (cpcol == npcol) exit
   cpcol = cpcol + 1
end do

call SLTIMER(6)
if (npcol > 1) then
   call DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, psi, nxp1, 0, 0)
end if
call SLTIMER(6)

if (mycol == 0) then
   call p_BCS(psi)
end if

return
END SUBROUTINE p_poisson

END MODULE stream
