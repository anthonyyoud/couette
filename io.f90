MODULE io
implicit none

contains

FUNCTION itos(n)
!Function to convert an integer into a string of length 7
implicit none

character(7) :: itos
integer, intent(in) :: n
integer   :: i, n_, d(7)
character :: c(0:9) = (/'0','1','2','3','4','5','6','7','8','9'/)

n_ = n
do i = 7, 1, -1
   d(i) = mod(n_,10)
   n_ = n_ / 10
end do

itos = c(d(1))//c(d(2))//c(d(3))//c(d(4))//c(d(5))//c(d(6))//c(d(7))

return
END FUNCTION itos

SUBROUTINE open_files()
!Open runtime files
use parameters
implicit none

open (20, status = 'unknown', file = 'u_growth.dat')
open (22, status = 'unknown', file = 'max_psi.dat')
open (24, status = 'unknown', file = 'particle.dat')
open (33, status = 'unknown', file = 'torque.dat')
if (save3d) open (35, status = 'unknown', file = 'isosurface.dat')
open (51, file = 'time_tau.dat')
open (99, file = 'RUNNING')
close (99)

return
END SUBROUTINE open_files

SUBROUTINE close_files()
!Close runtime files
use parameters
implicit none

close (20)
close (22)
close (24)
close (33)
if (save3d) close (35)
close (51)

return
END SUBROUTINE close_files

SUBROUTINE save_growth(t, ur, ur_prev, uz, pn, v, zn, bn, jn, growth)
!Save fields at particular points in (x,z)-plane
use parameters
implicit none

double precision, intent(in) :: t, ur(0:nx,0:nz), uz(0:nx,0:nz), &
                                pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                jn(0:nx,0:nz), v(0:nx,0:nz), &
                                zn(0:nx,0:nz), ur_prev(0:nx,0:nz)
double precision, intent(out) :: growth
double precision, save :: min_p, max_p, min_ur, max_ur, min_uz, max_uz
integer :: zpos, xpos

!growth rate of vortices
growth = log(abs(ur(nx/2,nz/2)/ur_prev(nx/2,nz/2))) / (dt * save_rate)

xpos = nx/2
zpos = nz/2   !position at which to save fields

write(20, '(11e17.9)') t, ur(nx/2,nz/2), ur(nx/2,0), growth, &
                      uz(xpos,zpos), &
                      pn(nx/4,3*nz/4), v(nx/2,nz/2), &
                      zn(nx/2,nz/4), bn(nx/2,nz/4), jn(nx/2,nz/2), &
                      Re1 + Re1_mod * dcos(om1 * t)

if (maxval(ur) > max_ur) then
   max_ur = maxval(ur)
end if
if (minval(ur) < min_ur) then
   min_ur = minval(ur)
end if

if (maxval(uz) > max_uz) then        !calculate maximum/minimum
   max_uz = maxval(uz)               !values of fields
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
!Save torques
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
!Save cross-sections of fields for use in IDL
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

SUBROUTINE save_3d(u_r, u_t, u_z, pn, p)
use parameters
use ic_bc
implicit none

integer, intent(in) :: p
double precision, intent(in) :: u_r(0:nx,0:nz), u_t(0:nx,0:nz), &
                                u_z(0:nx,0:nz), pn(0:nx,0:nz)
double precision :: hel(0:nx,0:nz)
integer :: j, k, l

!open (35, status = 'unknown', file = 'p3d'//itos(p)//'.dat')

if (iso_hel) call helicity(u_r, u_t, u_z, hel)

write(35,*) 'nx, nt, nz = ', nx, nt, nz, p

if (iso_hel) then
   do j = 0, nx
      do l = 0, nt
         do k = 0, nz
            write(35,'(4e11.3)') x_(j) * dcos(th(l)), x_(j) * dsin(th(l)), &
	                         z(k), hel(j,k)
         end do
      end do
   end do
else
   do j = 0, nx
      do l = 0, nt
         do k = 0, nz
            write(35,'(4e11.3)') x_(j) * dcos(th(l)), x_(j) * dsin(th(l)), &
	                         z(k), pn(j,k)
         end do
      end do
   end do
end if

write(35,*)

!close(35)

return
END SUBROUTINE save_3d

SUBROUTINE helicity(u_r, u_t, u_z, hel)
use parameters
use ic_bc
use derivs
implicit none

double precision, intent(in) :: u_r(0:nx,0:nz), u_t(0:nx,0:nz), u_z(0:nx,0:nz)
double precision, intent(out) :: hel(0:nx, 0:nz)
double precision :: u_r_z(0:nx,0:nz), u_t_z(0:nx,0:nz), &
                    u_t_x(0:nx,0:nz), u_z_x(0:nx,0:nz)
integer :: k

call deriv_z(u_r, u_r_z)
call deriv_z(u_t, u_t_z)
call deriv_x(u_t, u_t_x)
call deriv_x(u_z, u_z_x)

do k = 0, nz
   hel(:,k) = 0.5d0 * delx * (-u_r(:,k) * u_t_z(:,k) + &
                               u_t(:,k) * u_r_z(:,k) - &
                               u_t(:,k) * u_z_x(:,k) + &
                               u_z(:,k) * u_t_x(:,k)) + &
                               (s(:) * u_t(:,k) * u_z(:,k)) / (1d0 - eta)
end do

return
END SUBROUTINE helicity

SUBROUTINE save_surface(pn, v, zn, ur, uz, bn, jn, p, t)
!Save surfaces of fields for use in gnuplot
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

SUBROUTINE write_data(p, p_start, t)
!Calculations of various quantities; save fields
use parameters
use variables
implicit none

integer, intent(in) :: p, p_start
double precision, intent(in) :: t
double precision :: growth_rate

call vr_vz(psi%old, vr, vz)   !get radial, axial velocities
if (save_part) call particle(vr, vrold, vz, vzold, x_pos, z_pos) !save particle
!if ((Re1 /= 0d0) .or. (Re2 /= 0d0)) then                        !path
!   call save_torque(t, unew)
!end if
if ((p /= p_start) .and. ((p - p_start) > save_rate)) then
   call save_growth(t, vr, vrold, vz, psi%old, ut%new, zt%new, &
                    bt%old, jt%old, growth_rate)
   if ((om1 == 0d0) .and. (om2 == 0d0)) then
      if ((dabs(growth_rate) < 1d-8) .and. &  !if vr saturated
          (dabs(vr(nx/2, nz/2)) > 1d-3)) then
         if ((.not. auto_tau) .or. (tau == tau_end)) then  !if tau not auto
            call save_time_tau(tau, t)                     !or tau at end
            call end_state(ut%old, zt%old, psi%old, bt%old, jt%old, p) !finish
         else if (tau < 1d0) then
            call save_time_tau(tau, t)
            tau = tau + tau_step   !increment tau
            print*, 'tau = ', tau
            call save_xsect(vr, vz, psi%old, t, p)
            call save_surface(psi%old, ut%new, zt%new, vr, vz, &
                              bt%old, jt%old, p, t)
         end if
      end if
   end if
end if

return
END SUBROUTINE write_data

SUBROUTINE terminate(p, t)
!Terminate run if file 'RUNNING' does not exist in run directory
use parameters
use variables
implicit none

integer, intent(in) :: p
double precision, intent(in) :: t
logical :: run_exist

if (mycol == 0) then
   inquire(file='RUNNING', exist=run_exist)  !does 'RUNNING' exist?
   if (.not. run_exist) then  !if not then finish
      print*, 'Stop requested.  Ending process ', mycol
      print*, 'Saving end state'
      call end_state(ut%old, zt%old, psi%old, bt%old, jt%old, p)
      call save_xsect(vr, vz, psi%old, t, p)
      call save_surface(psi%old, ut%old, zt%old, &
                        vr, vz, bt%old, jt%old, p, t)
      end_proc = 1  !flag to send to all processes
   end if
end if

call SLTIMER(3)                                          !tell all other
call IGEBR2D(ictxt, 'A', ' ', 1, 1, end_proc, 1, 0, 0)   !processes to finish
call SLTIMER(3)                                          !if end_proc=1

return
END SUBROUTINE terminate

SUBROUTINE save_run(p, t)
!Save cross-sections, surfaces if file 'SAVE' exists in run directory
use parameters
use variables
implicit none

integer, intent(in) :: p
double precision, intent(in) :: t
logical :: save_exist

inquire(file='SAVE', exist=save_exist)   !does 'SAVE' exist?
if (save_exist) then   !if so then save fields
   call save_xsect(vr, vz, psi%old, t, p)
   call save_surface(psi%old, ut%old, zt%old, &
                     vr, vz, bt%old, jt%old, p, t)
   open (98, file = 'SAVE')
   close (98, status = 'delete')
end if

return
END SUBROUTINE save_run

SUBROUTINE end_state(u, zn, pn, bn, jn, p)
!Save variables for use in a restarted run
use parameters
implicit none

double precision, intent(in) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                jn(0:nx,0:nz)
integer, intent(in) :: p
integer :: j, k

open (50, file = 'end_state.dat')

write(50, '(i8)') p
write(50, '(e19.7)') dt
write(50, '(e19.7)') ((u(j,k), k = 0, nz), j = 0, nx)
write(50, '(e19.7)') ((zn(j,k), k = 0, nz), j = 0, nx)
write(50, '(e19.7)') ((pn(j,k), k = 0, nz), j = 0, nx)
write(50, '(e19.7)') ((bn(j,k), k = 0, nz), j = 0, nx)
write(50, '(e19.7)') ((jn(j,k), k = 0, nz), j = 0, nx)

close (50)
open (99, file = 'RUNNING')  !delete 'RUNNING' to finish run
close (99, status = 'delete')

return
END SUBROUTINE end_state

SUBROUTINE save_time_tau (tau, t)
!Save the times at which tau is incremented
implicit none

double precision, intent(in) :: tau, t

write(51, '(2e17.9)') t, tau

END SUBROUTINE save_time_tau

SUBROUTINE particle (vr, vrold, vz, vzold, xold, zold)
!Calculate and save a particle path
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
