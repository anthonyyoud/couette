MODULE magnetic 
implicit none

contains

SUBROUTINE b_poisson(u_mat, bn, b_mat, IPIV)
use parameters
use ic_bc
implicit none

double precision, intent(in) :: u_mat(0:nx,0:nz)
double precision, intent(in) :: b_mat(2*nxp1+nxp1+1,0:nxp1*nz1-1)
double precision, intent(out) :: bn(0:nx,0:nz)
double precision :: u_vec(0:nxp1*nz1-1)
integer, intent(in) :: IPIV(nxp1*nz1)
integer :: j, k, info

do k = 1, nz1
   do j = 0, nx
      u_vec(nxp1*(k-1)+j) = 0.5d0 * dx2 * delz * &
                            (u_mat(j,k-1) - u_mat(j,k+1))
   end do
end do

call DGBTRS('N', nxp1*nz1, nxp1, nxp1, 1, b_mat, 2*nxp1+nxp1+1, &
             IPIV, u_vec, nxp1*nz1, info)

do k = 1, nz1
   do j = 0, nx
      bn(j,k) = u_vec(nxp1*(k-1)+j)
   end do
end do

call b_BCS(bn)

return
END SUBROUTINE b_poisson

SUBROUTINE fin_b_poisson(u_mat, bn, b_mat, IPIV)
use parameters
use ic_bc
use derivs
implicit none

double precision, intent(in) :: u_mat(0:nx,0:nz)
double precision, intent(in) :: b_mat(2*nxp1+nxp1+1,0:nxp1*nzp1-1)
double precision, intent(out) :: bn(0:nx,0:nz)
double precision :: u_vec(0:nxp1*nzp1-1), u_mat_z(0:nx,0:nz)
integer, intent(in) :: IPIV(nxp1*nzp1)
integer :: j, k, info

call deriv_z(u_mat, u_mat_z)

do k = 0, nz
   do j = 0, nx
      u_vec(nxp1*k+j) = -0.5d0 * dx2 * delz * u_mat_z(j,k)
   end do
end do

call DGBTRS('N', nxp1*nzp1, nxp1, nxp1, 1, b_mat, 2*nxp1+nxp1+1, &
             IPIV, u_vec, nxp1*nzp1, info)

do k = 0, nz
   do j = 0, nx
      bn(j,k) = u_vec(nxp1*k+j)
   end do
end do

return
END SUBROUTINE fin_b_poisson

END MODULE magnetic
