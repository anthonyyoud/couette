MODULE pressure
implicit none

contains

SUBROUTINE poisson(Z_mat, psi, AB, IPIV, s)
use parameters
!use matrices
implicit none

logical, parameter :: write_file = .false.
double precision, intent(in) :: Z_mat(0:nx,0:nz), s(0:nx)
double precision, intent(in) :: AB(2*nx1+nx1+1,nx1*nz1)
double precision, intent(out) :: psi(0:nx,0:nz)
double precision :: zvec(nx1*nz1) !, save_AB(2*nx1+nx1+1,nx1*nz1)
integer, intent(in) :: IPIV(nx1*nz1)
integer :: j, k, info
double precision :: x(0:nx), z(0:nz)

do j = 0, nx
   x(j) = j * delx
end do
do k = 0, nz
   z(k) = k * delz
end do


if (write_file) then
   open (61, file = 'zvecin.dat')
   open (62, file = 'zvecout.dat')
end if

do k = 1, nz1
   do j = 1, nx1
      zvec(nx1*(k-1)+j) = -s(j)**2 * dx2 * dz2 * Z_mat(j,k)
!      zvec(nx1*(k-1)+j) = -2d0 * pi**2 * dx2 * dz2 * &
!                          sin(pi*x(j)) * sin(pi*z(k))
   end do
end do

if (write_file) then
   write(61, '(e17.9)') (zvec(j), j = 1, nx1*nz1)
end if

!save_AB = AB

!call DGBSV(nx1*nz1, nx1, nx1, 1, AB, 2*nx1+nx1+1, &
!           IPIV, zvec, nx1*nz1, info)

call DGBTRS('N', nx1*nz1, nx1, nx1, 1, AB, 2*nx1+nx1+1, &
             IPIV, zvec, nx1*nz1, info)

!AB = save_AB

if (write_file) then
   write(62, '(e17.9)') (zvec(j), j = 1, nx1*nz1)
end if

do k = 1, nz1
   do j = 1, nx1
      psi(j,k) = zvec(nx1*(k-1)+j)
   end do
end do

psi(0,:) = 0d0
psi(nx,:) = 0d0
psi(:,0) = 0d0
psi(:,nz) = 0d0

if (write_file) then
   close (61)
   close (62)
end if

END SUBROUTINE poisson

SUBROUTINE p_BCS(p)
use parameters
implicit none
double precision, intent(out) :: p(0:nx,0:nz)

p(0,:) = 0d0
p(nx,:) = 0d0

p(:,0) = 0d0
p(:,nz) = 0d0

return
END SUBROUTINE p_BCS

END MODULE pressure
