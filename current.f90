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
      p_vec(nx1*k+j) = dx2 * dz2 * (0.5d0 * dp%zzz(j,k) / (s(j) * delz**3) + &
                       0.5d0 * dp%zxx(j,k) / (s(j) * dx2 * delz) - &
                       0.25d0 * (1d0 - eta) * dp%zx(j,k) / &
                       (s(j)**2 * delx * delz))
   end do
end do

!do k = 0, nz
!   do j = 1, nx1
!      p_vec(nx1*k+j) = dx2 * dz2 * ( &
!                       -(pi**2+alpha**2+((1d0-eta)**2)/(s(j)**2))*&
!                       sin(pi*x(j))*cos(alpha*z(k))+&
!                       (1d0-eta)*pi*cos(pi*x(j))*cos(alpha*z(k))/s(j))
!   end do
!end do
!write(6,'(e17.9)') p_vec
!write(6,*)
!write(6,'(e17.9)') j_mat(2*nx-1,1), j_mat(2*nx-2,2), j_mat(nx,1+nx1)
!write(6,*)
call DGBTRS('N', nx1*nzp1, nx1, nx1, 1, j_mat, 2*nx1+nx1+1, &
             IPIV, p_vec, nx1*nzp1, info)
!write(6,'(e17.9)') p_vec
!write(6,*)

do k = 0, nz
   do j = 1, nx1
      jn(j,k) = p_vec(nx1*k+j)
   end do
end do

!open(73,file='j_test.dat')
!do j = 0, nx
!   write(73, '(741e19.7)') (x(j), z(k), jn(j,k), k = 0, nz)
!   write(73, *)
!end do
!close(73)

call j_BCS(jn)

return
END SUBROUTINE j_poisson

SUBROUTINE fin_j_poisson(p_mat, jn, j_mat, IPIV)
use parameters
use io
use ic_bc
use derivs
implicit none

double precision, intent(in) :: p_mat(0:nx,0:nz)
double precision, intent(in) :: j_mat(2*nx1+nx1+1,nx1*nz1)
double precision, intent(out) :: jn(0:nx,0:nz)
double precision :: p_vec(nx1*nz1)
integer, intent(in) :: IPIV(nx1*nz1)
integer :: j, k, info
type (deriv) :: dp

call deriv_z(p_mat, dp%z)
call deriv_x(dp%z, dp%zx)
call deriv_xx(dp%z, dp%zxx)
call deriv_zz(dp%z, dp%zzz)

do k = 1, nz1
   do j = 1, nx1
      p_vec(nx1*(k-1)+j) = dx2 * dz2 * (0.5d0 * dp%zzz(j,k) / &
                       (s(j) * delz**3) + &
                       0.5d0 * dp%zxx(j,k) / (s(j) * dx2 * delz) - &
                       0.25d0 * (1d0 - eta) * dp%zx(j,k) / &
                       (s(j)**2 * delx * delz))
   end do
end do

!do k = 0, nz
!   do j = 1, nx1
!      p_vec(nx1*k+j) = dx2 * dz2 * ( &
!                       -(pi**2+alpha**2+((1d0-eta)**2)/(s(j)**2))*&
!                       sin(pi*x(j))*cos(alpha*z(k))+&
!                       (1d0-eta)*pi*cos(pi*x(j))*cos(alpha*z(k))/s(j))
!   end do
!end do
!write(6,'(e17.9)') p_vec
!write(6,*)
!write(6,'(e17.9)') j_mat(2*nx-1,1), j_mat(2*nx-2,2), j_mat(nx,1+nx1)
!write(6,*)
call DGBTRS('N', nx1*nz1, nx1, nx1, 1, j_mat, 2*nx1+nx1+1, &
             IPIV, p_vec, nx1*nz1, info)
!write(6,'(e17.9)') p_vec
!write(6,*)

do k = 1, nz1
   do j = 1, nx1
      jn(j,k) = p_vec(nx1*(k-1)+j)
   end do
end do

!open(73,file='j_test.dat')
!do j = 0, nx
!   write(73, '(741e19.7)') (x(j), z(k), jn(j,k), k = 0, nz)
!   write(73, *)
!end do
!close(73)

call j_BCS(jn)

return
END SUBROUTINE fin_j_poisson

END MODULE current
