MODULE io
use parameters
implicit none

type var
   double precision :: new(0:nx, 0:nz)
   double precision :: old(0:nx, 0:nz)
   double precision :: old2(0:nx, 0:nz)
   double precision :: int(0:nx, 0:nz)
   double precision :: nlin_new(0:nx, 0:nz)
   double precision :: nlin_old(0:nx, 0:nz)
end type var

type deriv
   double precision :: x(0:nx, 0:nz)
   double precision :: xx(0:nx, 0:nz)
   double precision :: z(0:nx, 0:nz)
   double precision :: zz(0:nx, 0:nz)
   double precision :: zx(0:nx, 0:nz)
   double precision :: zxx(0:nx, 0:nz)
   double precision :: zzz(0:nx, 0:nz)
end type deriv

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
open (22, status = 'unknown', file = 'max_psi.dat')
open (24, status = 'unknown', file = 'particle.dat')
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
close (22)
close (24)
close (33)
close (51)
if (diag) then
   close (90)
   close (91)
   close (92)
end if

return
END SUBROUTINE close_files

SUBROUTINE save_growth(t, ur, ur_prev, uz, pn, v, zn, bn, jn, growth)
use parameters
implicit none
double precision, intent(in) :: t, ur(0:nx,0:nz), uz(0:nx,0:nz), &
                                pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                jn(0:nx,0:nz), v(0:nx,0:nz), &
                                zn(0:nx,0:nz), ur_prev(0:nx,0:nz)
double precision, intent(out) :: growth
double precision, save :: min_p, max_p, min_ur, max_ur, min_uz, max_uz
integer :: zpos, xpos

growth = log(abs(ur(nx/2,nz/2)/ur_prev(nx/2,nz/2))) / (dt * save_rate)

xpos = nx/2
zpos = nz/2

write(20, '(11e17.9)') t, ur(nx/2,nz/2), ur(nx/2,0), growth, &
                      uz(xpos,zpos), &
                      pn(nx/4,3*nz/4), v(nx/2,nz/2), &
                      zn(nx/2,nz/4), bn(nx/2,nz/2), jn(nx/2,nz/2), &
                      Re1 + Re1_mod * dcos(om1 * t)

if (maxval(ur) > max_ur) then
   max_ur = maxval(ur)
end if
if (minval(ur) < min_ur) then
   min_ur = minval(ur)
end if

if (maxval(uz) > max_uz) then
   max_uz = maxval(uz)
end if
if (minval(uz) < min_uz) then
   min_uz = minval(uz)
end if

if (maxval(pn) > max_p) then
   max_p = maxval(pn)
end if
if (minval(pn) < min_p) then
   min_p = minval(pn)
end if

write (22, '(7e17.9)') t, max_p, min_p, max_ur, min_ur, max_uz, min_uz

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

SUBROUTINE save_xsect(ur, uz, pn, t, p)
use parameters
use ic_bc
implicit none

integer, intent(in) :: p
double precision, intent(in) :: ur(0:nx,0:nz), uz(0:nx,0:nz), &
                                pn(0:nx,0:nz), t
integer :: j, k

open (32, status = 'unknown', file = 'xsect'//itos(p)//'.dat')

write (32, '(e17.9)') t
write (32, '(2i5)') nx, nz
write (32, '(e17.9)') ((ur(j,k), j = 0, nx), k = 0, nz)
write (32, '(e17.9)') ((uz(j,k), j = 0, nx), k = 0, nz)
write (32, '(e17.9)') ((pn(j,k), j = 0, nx), k = 0, nz)
write (32, '(e17.9)') (x(j), j = 0, nx)
write (32, '(e17.9)') (z(k), k = 0, nz)

close (32)

return
END SUBROUTINE

SUBROUTINE save_surface(pn, v, zn, ur, uz, bn, jn, p, t)
use parameters
use ic_bc
implicit none
integer, intent(in) :: p
double precision, intent(in) :: t, pn(0:nx,0:nz), v(0:nx,0:nz), &
                                zn(0:nx,0:nz), bn(0:nx,0:nz), &
                                jn(0:nx,0:nz), ur(0:nx,0:nz), &
                                uz(0:nx,0:nz)
integer :: j, k

open (19, status = 'unknown', file = 'p'//itos(p)//'.dat')
!open (70, status = 'unknown', file = 'sin.dat')
open (21, status = 'unknown', file = 'z'//itos(p)//'.dat')
open (23, status = 'unknown', file = 'u'//itos(p)//'.dat')
open (30, status = 'unknown', file = 'vr'//itos(p)//'.dat')
open (31, status = 'unknown', file = 'vz'//itos(p)//'.dat')
open (34, status = 'unknown', file = 'b'//itos(p)//'.dat')
open (35, status = 'unknown', file = 'j'//itos(p)//'.dat')

write(30, '(2A, i10, e19.7)') '#', 'p=', p, t
write(31, '(2A, i10, e19.7)') '#', 'p=', p, t
write(34, '(2A, i10, e19.7)') '#', 'p=', p, t
write(35, '(2A, i10, e19.7)') '#', 'p=', p, t
write(23, '(2A, i10, e19.7)') '#', 'p=', p, t
write(19, '(2A, i10, e19.7)') '#', 'p=', p, t
write(21, '(2A, i10, e19.7)') '#', 'p=', p, t
!write(70, '(2A, i10, e19.7)') '#', 'p=', p, t

do j = 0, nx
   write(30, '(3e19.7)') (x(j), z(k), ur(j,k), k = 0, nz)
   write(30, *)
   write(31, '(3e19.7)') (x(j), z(k), uz(j,k), k = 0, nz)
   write(31, *)
   write(34, '(3e19.7)') (x(j), z(k), bn(j,k), k = 0, nz)
   write(34, *)
   write(35, '(3e19.7)') (x(j), z(k), jn(j,k), k = 0, nz)
   write(35, *)
   write(23, '(3e19.7)') (x(j), z(k), v(j,k), k = 0, nz)
   write(23, *)
   write(19, '(3e19.7)') (x(j), z(k), pn(j,k), k = 0, nz)
   write(19, *)
   write(21, '(3e19.7)') (x(j), z(k), zn(j,k), k = 0, nz)
   write(21, *)
!   write(70, '(3e19.7)') (x(j), z(k), dsin(pi*x(j))*&
!                          dsin(alpha*z(k)), k = 0, nz)
!   write(70, *)
end do

close (19)
close (21)
close (23)
close (30)
close (31)
close (34)
close (35)
!close (70)

return
END SUBROUTINE save_surface

SUBROUTINE end_state(u, zn, pn, bn, jn, p)
use parameters
implicit none
double precision, intent(in) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                jn(0:nx,0:nz)
integer, intent(in) :: p
integer :: j, k

open (50, file = 'end_state.dat')

write(50, '(i7)') p
write(50, '(e19.7)') dt
write(50, '(e19.7)') ((u(j,k), k = 0, nz), j = 0, nx)
write(50, '(e19.7)') ((zn(j,k), k = 0, nz), j = 0, nx)
write(50, '(e19.7)') ((pn(j,k), k = 0, nz), j = 0, nx)
write(50, '(e19.7)') ((bn(j,k), k = 0, nz), j = 0, nx)
write(50, '(e19.7)') ((jn(j,k), k = 0, nz), j = 0, nx)

close (50)
open (99, file = 'RUNNING')
close (99, status = 'delete')

return
END SUBROUTINE end_state

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

SUBROUTINE particle (vr, vrold, vz, vzold, xold, zold)
use parameters
implicit none

double precision, intent(in) :: vr(0:nx, 0:nz), vz(0:nx, 0:nz), &
                                vrold(0:nx, 0:nz), vzold(0:nx, 0:nz)
double precision, intent(inout) :: xold, zold
integer :: xmin, xplu, zmin, zplu, j
double precision :: c1, c2, rvel, zvel, xnew, znew, del_t

del_t = dt / 1d0

do j = 1, 1
xmin = int(xold)
xplu = int(xold + 1d0)
zmin = int(zold)
zplu = int(zold + 1d0)

c1 = (xold - xmin) / (xplu - xmin)
c2 = (zold - zmin) / (zplu - zmin)

rvel = (1d0 - c1) * (1d0 - c2) * vrold(xmin, zmin) + &
        c1 * (1d0 - c2) * vrold(xplu, zmin) + &
        c1 * c2 * vrold(xplu, zplu) + &
        (1d0 - c1) * c2 * vrold(xmin, zplu)

zvel = (1d0 - c1) * (1d0 - c2) * vzold(xmin, zmin) + &
        c1 * (1d0 - c2) * vzold(xplu, zmin) + &
        c1 * c2 * vzold(xplu, zplu) + &
        (1d0 - c1) * c2 * vzold(xmin, zplu)

xnew = xold + del_t * rvel
znew = zold + del_t * zvel
xold = xnew
zold = znew
end do

write (24, '(2e17.9)') xnew / nx, znew * gamma / nz

return
END SUBROUTINE particle

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
