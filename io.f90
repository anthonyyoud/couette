MODULE io
use parameters
implicit none

type grid
   double precision :: x(0:nx)
   double precision :: z(0:nz)
end type grid

type var
   double precision :: new(0:nx,0:nz)
   double precision :: old(0:nx,0:nz)
   double precision :: old2(0:nx,0:nz)
   double precision :: int(0:nx,0:nz)
end type var

type mat_comp
   double precision :: lo(2:nx1)
   double precision :: di(nx1)
   double precision :: up(nx-2)
end type mat_comp

type uz_mat_comp
   double precision :: lo(nz)
   double precision :: di(0:nz)
   double precision :: up(0:nz1)
end type uz_mat_comp   

type zz_mat_comp
   double precision :: lo(2:nz1)
   double precision :: di(nz1)
   double precision :: up(nz-2)
end type zz_mat_comp

contains

!SUBROUTINE grid(x, z, s)
!use parameters
!implicit none
!integer :: j, k
!double precision, save :: x(0:nx), z(0:nz), s(0:nx)

!do j = 0, nx
!   x(j) = dble(j) * delx
!   s(j) = eta + ((1d0 - eta) * x(j))
!   do k = 0, nz
!      z(k) = dble(k) * delz
!   end do
!end do

!END SUBROUTINE grid

FUNCTION itos(n)
implicit none
character(7) :: itos
integer, intent(in) :: n
integer   :: i, n_, d(7)
character :: c(0:9) =  &
   (/'0','1','2','3','4','5','6','7','8','9'/)

n_ = n
do i = 7, 1, -1
   d(i) = mod(n_,10)
   n_ = n_ / 10
end do

itos = c(d(1))//c(d(2))//c(d(3))//c(d(4))//c(d(5))//c(d(6))//c(d(7))

return
END FUNCTION itos

SUBROUTINE open_files()
use parameters
implicit none
open (20, status = 'unknown', file = 'u_growth.dat')
open (33, status = 'unknown', file = 'torque.dat')
open (51, file = 'time_tau.dat')
if (diag) then
   open (90, status = 'unknown', file = 'lhs_u.dat')
   open (91, status = 'unknown', file = 'lhs_Z.dat')
   open (92, status = 'unknown', file = 'lhs_psi.dat')
end if
open (99, file = 'RUNNING')
close (99)

return
END SUBROUTINE open_files

SUBROUTINE close_files()
use parameters
implicit none
close (20)
close (33)
close (51)
if (diag) then
   close (90)
   close (91)
   close (92)
end if

return
END SUBROUTINE close_files

SUBROUTINE save_growth(t, ur, ur_prev, uz, pn, v, z, growth)
use parameters
implicit none
double precision, intent(in) :: t, ur(0:nx,0:nz), uz(0:nx,0:nz), &
                                pn(0:nx,0:nz), v(0:nx,0:nz), &
                                z(0:nx,0:nz), ur_prev(0:nx,0:nz)
double precision, intent(out) :: growth
integer :: zpos

if ((Re1_mod == 0d0) .and. (Re2_mod == 0d0) .and. &
    (om1 == 0d0) .and. (om2 == 0d0)) then
   growth = log(abs(ur(nx/2,nz/2)/ur_prev(nx/2,nz/2))) / (dt * save_rate)
end if

zpos = nz - (nz / (2 * gamma)) !(nz * (gamma - 1)) / (2 * gamma)

write(20, '(8e19.7)') t, ur(nx/2,nz/2), growth, uz(nx/10,zpos), &
                      pn(nx/2,3*nz/4), v(nx/2,nz/2), &
                      z(nx/2,nz/4), Re1 + Re1_mod * dcos(om1 * t)

return
END SUBROUTINE save_growth

SUBROUTINE save_torque(t, v)
use parameters
implicit none
double precision, intent(in) :: t, v(0:nx,0:nz)
integer :: k
double precision :: xi, C1, C2, G1(0:nz), G2(0:nz), G1_, G2_

xi = Re1_mod * cos(om1 * t) - eta * Re2_mod * cos(om2 * t)
xi = xi / (Re1 - eta * Re2)

C1 = (-2d0 * (1d0 + xi)) / (eta * (1d0 + eta))

C2 = 1d0 / (Re1 - eta * Re2)

   G1(:) = C1 + C2 * (0.5d0 * (4d0 * v(1,:) - v(2,:))) / delx
   G2(:) = C1 + (C2 / eta**2) * (0.5d0 * (v(nx-2,:) - 4d0 * &
                                    v(nx-1,:))) / delx

G1_ = sum(G1(:))
G2_ = sum(G2(:))
!print*,C1,C2
write (33, '(3e17.9)') t, G1_, G2_ !G1(nz/4), G2(nz/4)

return
END SUBROUTINE save_torque

SUBROUTINE save_xsect(ur, uz, x, z, p)
use parameters
implicit none

integer, intent(in) :: p
double precision, intent(in) :: ur(0:nx,0:nz), uz(0:nx,0:nz), &
                                x(0:nx), z(0:nz)
integer :: j, k

open (32, status = 'unknown', file = 'xsect'//itos(p)//'.dat')

write (32, '(2i5)') nx, nz
write (32, '(e17.9)') ((ur(j,k), j = 0, nx), k = 0, nz)
write (32, '(e17.9)') ((uz(j,k), j = 0, nx), k = 0, nz)
write (32, '(e17.9)') (x(j), j = 0, nx)
write (32, '(e17.9)') (z(k), k = 0, nz)

close (32)

return
END SUBROUTINE

SUBROUTINE save_surface(pn, v, zn, ur, uz, x, z, p, t)
use parameters
implicit none
integer, intent(in) :: p
double precision, intent(in) :: t, pn(0:nx,0:nz), v(0:nx,0:nz), &
                                zn(0:nx,0:nz), ur(0:nx,0:nz), &
                                uz(0:nx,0:nz), x(0:nx), z(0:nz)
integer :: j, k

open (19, status = 'unknown', file = 'p'//itos(p)//'.dat')
open (21, status = 'unknown', file = 'z'//itos(p)//'.dat')
open (23, status = 'unknown', file = 'u'//itos(p)//'.dat')
open (30, status = 'unknown', file = 'vr'//itos(p)//'.dat')
open (31, status = 'unknown', file = 'vz'//itos(p)//'.dat')

write(30, '(2A, i10, e19.7)') '#', 'p=', p, t
write(31, '(2A, i10, e19.7)') '#', 'p=', p, t
write(23, '(2A, i10, e19.7)') '#', 'p=', p, t
write(19, '(2A, i10, e19.7)') '#', 'p=', p, t
write(21, '(2A, i10, e19.7)') '#', 'p=', p, t

do j = 0, nx
   write(30, '(3e19.7)') (x(j), z(k), ur(j,k), k = 0, nz)
   write(30, *)
   write(31, '(3e19.7)') (x(j), z(k), uz(j,k), k = 0, nz)
   write(31, *)
   write(23, '(3e19.7)') (x(j), z(k), v(j,k), k = 0, nz)
   write(23, *)
   write(19, '(3e19.7)') (x(j), z(k), pn(j,k), k = 0, nz)
   write(19, *)
   write(21, '(3e19.7)') (x(j), z(k), zn(j,k), k = 0, nz)
   write(21, *)
end do

close (19)
close (21)
close (23)
close (30)
close (31)

return
END SUBROUTINE save_surface

SUBROUTINE end_state(u, zn, pn, p)
use parameters
implicit none
double precision, intent(in) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                pn(0:nx,0:nz)
integer, intent(in) :: p
integer :: j, k

open (50, file = 'end_state.dat')

write(50, '(i7)') p
write(50, '(e19.7)') ((u(j,k), k = 0, nz), j = 0, nx)
write(50, '(e19.7)') ((zn(j,k), k = 0, nz), j = 0, nx)
write(50, '(e19.7)') ((pn(j,k), k = 0, nz), j = 0, nx)

close (50)
open (99, file = 'RUNNING')
close (99, status = 'delete')

return
END SUBROUTINE end_state

SUBROUTINE state_restart(u, zn, pn, p)
use parameters
implicit none
double precision, intent(out) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                 pn(0:nx,0:nz)
integer, intent(out) :: p
integer :: j, k

open (50, file = 'end_state.dat')

read(50, *) p
read(50, *) ((u(j,k), k = 0, nz), j = 0, nx)
read(50, *) ((zn(j,k), k = 0, nz), j = 0, nx)
read(50, *) ((pn(j,k), k = 0, nz), j = 0, nx)

close (50)

return
END SUBROUTINE state_restart

SUBROUTINE save_time_tau (tau, t)
implicit none
double precision, intent(in) :: tau, t

write(51, '(2e17.9)') t, tau

END SUBROUTINE save_time_tau

SUBROUTINE thomas (lb, m, up, di, lo, r)
implicit none
integer :: j
integer, intent(in) :: m, lb
double precision, intent(in) :: up(lb:m-1), di(lb:m), lo(lb+1:m)
double precision, intent(inout) :: r(lb:m)
double precision :: dnew(lb:m), aa = 0d0

dnew = di
do j = lb+1, m
   aa = -lo(j) / dnew(j-1)
   dnew(j) = dnew(j) + aa * up(j-1)
   r(j) = r(j) + aa * r(j-1)
end do

r(m) = r(m) / dnew(m)

do j = m-1, lb, -1
   r(j) = (r(j) - up(j) * r(j+1)) / dnew(j)
end do

return
END SUBROUTINE thomas

SUBROUTINE get_timestep()
use parameters
implicit none

double precision :: timestep

if (Re1 * Re2 >= 0d0) then
   timestep = max(1d0, dabs(Re1) + dabs(Re1_mod), &
                  dabs(Re2) + dabs(Re2_mod))
else
   timestep = max(1d0, dabs(Re1) + dabs(Re1_mod) + &
                    dabs(Re2) + dabs(Re2_mod))
end if

timestep = 1d-4 * eta * 2d0 * pi / ((1d0 - eta) * timestep)

print*, timestep

END SUBROUTINE get_timestep

END MODULE io
