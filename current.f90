MODULE current 
implicit none

contains

SUBROUTINE j_poisson(p_mat, jn, j_mat, IPIV)
use parameters
use io
use ic_bc
use derivs
implicit none

double precision, intent(in) :: p_mat(0:nx,0:nz)
double precision, intent(in) :: j_mat(2*nx1+nx1+1,nx1*nzp1)
double precision, intent(out) :: jn(0:nx,0:nz)
double precision :: p_vec(nx1*nzp1)
integer, intent(in) :: IPIV(nx1*nzp1)
integer :: j, k, info
type (deriv) :: dp

call deriv_z(p_mat, dp%z)
call deriv_x(dp%z, dp%zx)
call deriv_xx(dp%z, dp%zxx)
call deriv_zz(dp%z, dp%zzz)

do k = 0, nz
   do j = 1, nx1
      p_vec(nx1*k+j) = dx2 * dz2 * 0.5d0 * dp%zzz(j,k) / (s(j) * delz**3) + &
                       0.5d0 * dp%zxx(j,k) / (s(j) * dx2 * delz) - &
                       0.25d0 * (1d0 - eta) * dp%zx(j,k) / &
                       (s(j)**2 * delx * delz)
   end do
end do

call DGBTRS('N', nx1*nzp1, nx1, nx1, 1, j_mat, 2*nx1+nx1+1, &
             IPIV, p_vec, nx1*nzp1, info)

do k = 0, nz
   do j = 1, nx1
      jn(j,k) = p_vec(nx1*k+j)
   end do
end do

call j_BCS(jn)

END SUBROUTINE j_poisson

END MODULE current
