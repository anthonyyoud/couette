MODULE pressure
implicit none

contains

SUBROUTINE poisson(Z_mat, psi, AB, IPIV)
use parameters
use ic_bc
implicit none

logical, parameter :: write_file = .false.
double precision, intent(in) :: Z_mat(0:nx,0:nz)
double precision, intent(in) :: AB(2*nx1+nx1+1,nx1*nz1)
double precision, intent(out) :: psi(0:nx,0:nz)
double precision :: zvec(nx1*nz1)
integer, intent(in) :: IPIV(nx1*nz1)
integer :: j, k, info

if (write_file) then
   open (61, file = 'zvecin.dat')
   open (62, file = 'zvecout.dat')
end if

do k = 1, nz1
   do j = 1, nx1
      zvec(nx1*(k-1)+j) = -s(j) * dx2 * dz2 * Z_mat(j,k)
!      zvec(nx1*(k-1)+j) = -2d0 * pi**2 * dx2 * dz2 * &
!                          sin(pi*x(j)) * sin(pi*z(k))
   end do
end do

if (write_file) then
   write(61, '(e17.9)') (zvec(j), j = 1, nx1*nz1)
end if

call DGBTRS('N', nx1*nz1, nx1, nx1, 1, AB, 2*nx1+nx1+1, &
             IPIV, zvec, nx1*nz1, info)

if (write_file) then
   write(62, '(e17.9)') (zvec(j), j = 1, nx1*nz1)
end if

do k = 1, nz1
   do j = 1, nx1
      psi(j,k) = zvec(nx1*(k-1)+j)
   end do
end do

call p_BCS(psi)

if (write_file) then
   close (61)
   close (62)
end if

END SUBROUTINE poisson

END MODULE pressure
