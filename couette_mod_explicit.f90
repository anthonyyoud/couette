PROGRAM couette_mod
use parameters
use pressure
use matrices
use io
use ccf
implicit none
logical :: file_exist, file_exist2
integer, parameter :: tfile = 11
type (mat_comp) :: Ux, Zx
type (uz_mat_comp) :: Uz
type (zz_mat_comp) :: Zz
!type (grid_setup) :: grid
double precision :: x(0:nx), z(0:nz), s(0:nx), growth_rate, &
!double precision :: s(0:nx), &
                    t = 0d0, A = 0d0, A_ = 0d0, B = 0d0, B_ = 0d0, &
u_nlin_new(0:nx,0:nz), u_nlin_old(0:nx,0:nz), &
z_nlin_new(0:nx,0:nz), z_nlin_old(0:nx,0:nz), &
unew(0:nx,0:nz), uold(0:nx,0:nz), uold2(0:nx,0:nz), u_int(0:nx,0:nz), &
znew(0:nx,0:nz), zold(0:nx,0:nz), zold2(0:nx,0:nz), z_int(0:nx,0:nz), &
pnew(0:nx,0:nz), pold(0:nx,0:nz), pold2(0:nx,0:nz), &
vr(0:nx,0:nz), vz(0:nx,0:nz), &
vc(0:nx), vc_(0:nx), vr2(0:nx,0:nz) = 0d0, &
AB(2*nx1+nx1+1,nx1*nz1), F(0:nx) !, A_mat(nx1*nz1,nx1*nz1)
integer :: pivot(nx1*nz1)

logical, parameter :: write_ofile = .true.
integer :: j, k, p = 0, p_start = 0

call open_files()
open (93, file = 'psi_deriv_test.dat')

print*, 'Setting up ICS...'
call ICS(unew, znew, pnew, x, z, s, p_start)

if (.not. restart) then
   print*, 'Setting up BCS...'
   call u_BCS(unew, 0d0)
   call p_BCS(pnew)
   call z_BCS(znew, pnew, s, 0d0)
end if

print*, 'Setting up matrices...'
!call matrix_setup(Ux, Uz, Zx, Zz, s)
!call A_mat_setup(A_mat, s)
call ABC_mat_setup(AB, pivot, s)

uold = unew
uold2 = unew
u_int = unew
zold = znew
zold2 = znew
z_int = znew
pold = pnew
pold2 = pnew

!call save_surface(pold, uold, zold, vr, vz, x, z, p, t)

print*, 'Entering time loop'

call get_timestep()

do p = p_start, Ntot
inquire(file='RUNNING', exist=file_exist)
file_exist2 = .not. file_exist
if (file_exist2) then
   print*, 'Stop requested'
   print*, 'Saving end state'
   call end_state(uold, zold, pold, p)
   call save_xsect(vr, vz, pold, x, z, p)
   call save_surface(pold, uold, zold, vr, vz, x, z, p, t)
   exit
end if

   t = p * dt

if (mod(p, save_rate) == 0) then
   call r_vel(pnew, s, vr, vz)
   !call save_torque(t, unew)
   if (p /= 0) then
      call save_growth(t, vr, vr2, vz, pnew, unew, znew, growth_rate)
   end if
end if

vr2 = vr

if (xsect_save) then
   if (mod(p, save_rate_2) == 0) then
      call save_xsect(vr, vz, pnew, x, z, p)
   end if
end if

call get_A(t, A, A_)
call get_B(t, B, B_)
call get_vc(A, A_, B_, B_, vc, vc_, s)
!call get_F(t, s, F)

!A = -eta*Re1/(1-eta**2)
!B =  eta*Re1/(1-eta**2)
!vc(:) = A*s(:)+B/s(:)
F=0d0

uold = unew
zold = znew
pold = pnew

call p_BCS(pnew)

call get_new_u(unew, uold, pold, A, F, s, t)
call get_new_Z(znew, zold, pold, uold, vc, s, t)

!do j = 0, nx
!   do k = 0, nz
!      znew(j,k) = ((2d0 * pi**2d0) / (s(j)**2d0)) * &
!                  sin(pi*x(j))*sin(pi*z(k)) + &
!                  ((pi * (1d0 - eta)) / (s(j)**3d0)) * &
!                  cos(pi*x(j))*sin(pi*z(k))
!      !pold(j,k) = 0d0 !sin(pi*x(j))*sin(pi*z(k))   
!   end do
!end do
!call psi_deriv_test(pnew, s)

call poisson(znew, pnew, AB, pivot, s)


if (diag) then
   call calc_rhs_u(unew, uold, pnew, s, t, A, p)
   call calc_rhs_Z(znew, zold, uold, pold, s, t, A, vc, p)
   call calc_rhs_psi(pnew, znew, s, t, p)
end if

if (p == Ntot) then
   call save_surface(pold, uold, zold, vr, vz, x, z, p, t)
   call end_state(uold, zold, pold, p)
end if

!if ((log(vr(nx/2,nz/2)/vr2(nx/2,nz/2)) / (dt * save_rate)) < 1d-6) then
!   print*, 'SATURATED'
!   call end_state(uold, zold, pold, p)
!   call save_xsect(vr, vz, x, z, p)
!   call save_surface(pold, uold, zold, vr, vz, x, z, p, t)
!   exit
!end if
   
end do
!*** END TIME LOOP***

call close_files()
close (93)

END PROGRAM couette_mod

SUBROUTINE ICS(u, zn, pn, x, z, s, p)
use parameters
use io
implicit none
double precision, intent(inout) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                   pn(0:nx,0:nz)
double precision, intent(out) :: x(0:nx), z(0:nz), s(0:nx)
real :: rand
integer, intent(out) :: p
integer :: j, k

do j = 0, nx
   x(j) = dble(j) * delx
   do k = 0, nz
      z(k) = dble(k) * delz
   end do
end do

s = eta + ((1d0 - eta) * x)

if (restart) then
   print*, 'Getting restart conditions'
   call state_restart(u, zn, pn, p)

else

   do k = 0, nz
      u(:,k) = seed * dsin(pi*x(:)) * dcos(alpha*z(k))
      pn(:,k) = seed * dsin(pi*x(:)) * dsin(alpha*z(k))
!      u(:,k) = seed * x(:)**2 * (1d0 - x(:))**2 * dcos(alpha*z(k))
!      pn(:,k) = seed * x(:)**2 * (1d0 - x(:))**2 * dsin(alpha*z(k))
   end do

   do k = 1, nz1
      do j = 1, nx1
         zn(j,k) = -(pn(j+1,k) - 2d0 * pn(j,k) + pn(j-1,k)) / &
                    (s(j)**2 * dx2) + &
                    0.5d0 * (1d0 - eta) * (pn(j+1,k) - pn(j-1,k)) / &
                    (s(j)**3 * delx) - &
                    (pn(j,k+1) - 2d0 * pn(j,k) + pn(j,k-1)) / &
                    (s(j)**2 * dz2)
      end do
   end do
end if

return
END SUBROUTINE ICS

SUBROUTINE u_BCS(u, t)
use parameters
implicit none
double precision, intent(out) :: u(0:nx,0:nz)
double precision, intent(in) :: t
integer :: k

!u(0,:) = 0d0
!u(nx,:) = 0d0

u(0,:) = Re1 + Re1_mod * dcos(om1 * t)
u(nx,:) = Re2 + Re2_mod * dcos(om2 * t)

return
END SUBROUTINE u_BCS

SUBROUTINE z_BCS(zn, pn, s, t)
use parameters
implicit none
double precision, intent(out) :: zn(0:nx,0:nz)
double precision, intent(in) :: t, s(0:nx), pn(0:nx,0:nz)
integer :: k

do k = 0, nz
zn(0,k) = (-1d0 / (dx2 * s(0)**2)) * (-pn(3,k) + 4d0 * pn(2,k) - &
          5d0 * pn(1,k) + 2d0 * pn(0,k)) !+ &
!          ((1d0 - eta) / (2d0 * delx * s(0)**3)) * (-pn(2,k) + 4d0 * &
!          pn(1,k) - 3d0 * pn(0,k)) - 0*(1d0 / (dz2 * s(0)**2)) * &
!          (pn(0,k+1) - 2d0 * pn(0,k) + pn(0,k-1))

zn(nx,k) = (-1d0 / (dx2 * s(nx)**2)) * (-pn(nx-3,k) + 4d0 * pn(nx-2,k) - &
          5d0 * pn(nx1,k) + 2d0 * pn(nx,k)) !+ &
!          ((1d0 - eta) / (2d0 * delx * s(nx)**3)) * (pn(nx-2,k) - 4d0 * 
!&
!          pn(nx1,k) + 3d0 * pn(nx,k)) - 0*(1d0 / (dz2 * s(nx)**2)) * &
!          (pn(nx,k+1) - 2d0 * pn(nx,k) + pn(nx,k-1))
end do


zn(0,:) = -(8d0 * pn(1,:) - pn(2,:)) / (2d0 * (eta**2) * dx2)
zn(nx,:) = -(8d0 * pn(nx1,:) - pn(nx-2,:)) / (2d0 * dx2)

!zn(0,:) = -(2d0 * pn(1,:)) / ((eta**2) * dx2)
!zn(nx,:) = -(2d0 * pn(nx-1,:)) / dx2

!zn(0,:) = 0d0
!zn(nx,:) = 0d0

zn(:,0) = 0d0 
zn(:,nz) = 0d0 

return
END SUBROUTINE z_BCS

SUBROUTINE r_vel(p, s, vr, vz)
use parameters
implicit none
double precision, intent(in) :: p(0:nx,0:nz), s(0:nx)
double precision, intent(out) :: vr(0:nx,0:nz), vz(0:nx,0:nz)
integer :: j, k

do k = 1, nz-1
   do j = 0, nx
      vr(j,k) = (-1d0 / (2d0 * s(j) * delz)) * (p(j,k+1) - p(j,k-1))
   end do
end do

do j = 0, nx
   vr(j,0) = (-1d0 / (2d0 * s(j) * delz)) * &
              (-3d0 * p(j,0) + 4d0 * p(j,1) - p(j,2))
   vr(j,nz) = (-1d0 / (2d0 * s(j) * delz)) * &
              (3d0 * p(j,nz) - 4d0 * p(j,nz-1) + p(j,nz-2))
end do

do k = 0, nz
   do j = 1, nx-1
      vz(j,k) = (1d0 / (2d0 * s(j) * delx)) * (p(j+1,k) - p(j-1,k))
   end do
end do

do k = 0, nz
   vz(0,k) = (1d0 / (2d0 * s(0) * delx)) * &
              (-3d0 * p(0,k) + 4d0 * p(1,k) - p(2,k))
   vz(nx,k) = (1d0 / (2d0 * s(nx) * delx)) * &
              (3d0 * p(nx,k) - 4d0 * p(nx-1,k) + p(nx-2,k))
end do

return
END SUBROUTINE r_vel

SUBROUTINE calc_rhs_u(up1, u, psi, s, t, A, p)
use parameters
use derivs
use ccf
implicit none
integer, intent(in) :: p
double precision, intent(in) :: t, up1(0:nx,0:nz), u(0:nx,0:nz), &
                                s(0:nx), A, &
                                psi(0:nx,0:nz)
double precision :: u_x(0:nx,0:nz), u_xx(0:nx,0:nz), &
                    u_z(0:nx,0:nz), u_zz(0:nx,0:nz), &
                    psi_x(0:nx,0:nz), psi_z(0:nx,0:nz)
double precision :: lhs_u(0:nx,0:nz)
integer :: j, k

call deriv_x(u, u_x)
call deriv_xx(u, u_xx)
call deriv_x(psi, psi_x)
call deriv_z(u, u_z)
call deriv_zz(u, u_z)
call deriv_z(psi, psi_z)

do j = 1, nx-1
   do k = 1, nz-1
!      lhs_u(j,k) = (u_xx(j,k) / dx2) + &
!                 ((0.5d0 * (1d0 - eta)) / (delx * s(j))) * &
!                 u_x(j,k) - &
!                 (((1d0 - eta)**2) / (s(j)**2)) * &
!                 u(j,k) + &
!                 (u_zz(j,k) / dz2) - &
!                 (0.25d0 / (s(j) * delx * delz)) * &
!                 ((psi_x(j,k) * u_z(j,k)) - &
!                 (psi_z(j,k) * u_x(j,k))) + &
!                 ((0.5d0 * (1d0 - eta)) / (delz * s(j))) * &
!                 ((u(j,k) / s(j)) + (2d0 * A) * psi_z(j,k))

      lhs_u(j,k) = (rxx * u_xx(j,k)) + &
                 (((1d0 - eta) * rx) / (2d0 * s(j))) * &
                 u_x(j,k) - &
                 ((dt * (1d0 - eta)**2) / (s(j)**2)) * &
                 u(j,k) + &
                 (rzz * u_zz(j,k)) - &
                 (rx / (4d0 * s(j) * delz)) * &
                 ((psi_x(j,k) * u_z(j,k)) - &
                 (psi_z(j,k) * u_x(j,k))) + &
                 (((1d0 - eta) * rz) / (2d0 * s(j))) * &
                 ((u(j,k) / s(j)) + (2d0 * A) * psi_z(j,k)) !- &
!                 up1(j,k) + u(j,k)
   end do
end do

if (mod(p, save_rate) == 0) then
   write (90, '(2e17.9)') t, lhs_u(nx/2,nz/2)
end if

return
END SUBROUTINE calc_rhs_u

SUBROUTINE calc_rhs_Z(zp1, z, u, psi, s, t, A, vc, p)
use parameters
use derivs
use ccf
implicit none

integer, intent(in) :: p
double precision, intent(in) :: t, A, zp1(0:nx,0:nz), z(0:nx,0:nz), &
                                u(0:nx,0:nz), &
                                psi(0:nx,0:nz), &
                                s(0:nx), vc(0:nx)

double precision :: u_z(0:nx,0:nz), &
                    psi_x(0:nx,0:nz), psi_z(0:nx,0:nz), &
                    z_x(0:nx,0:nz), &
                    z_xx(0:nx,0:nz), z_z(0:nx,0:nz), &
                    z_zz(0:nx,0:nz)
double precision :: lhs_z(0:nx,0:nz)
integer :: j, k

call deriv_x(z, z_x)
call deriv_xx(z, z_xx)
call deriv_x(psi, psi_x)
call deriv_z(z, z_z)
call deriv_zz(z, z_zz)
call deriv_z(u, u_z)
call deriv_z(psi, psi_z)

do j = 1, nx-1
   do k = 1, nz-1
      lhs_z(j,k) = rxx * (z_xx(j,k)) + &
                   ((1.5d0 * rx * (1d0 - eta)) / s(j)) * &
                   (z_x(j,k)) + &
                   rzz * (z_zz(j,k)) + &
                   (((1d0 - eta) * rz) / (s(j)**2)) * &
                   ((vc(j) + u(j,k)) * u_z(j,k)) - &
                   (rx / (4d0 * s(j) * delz)) * &
                   ((psi_x(j,k) * z_z(j,k) - &
                   psi_z(j,k) * z_x(j,k))) + &
                   z(j,k) - zp1(j,k)
   end do
end do

if (mod(p, save_rate) == 0) then
   write (91, '(2e17.9)') t, lhs_z(nx/2,nz/4)
end if

return
END SUBROUTINE calc_rhs_Z

SUBROUTINE calc_rhs_psi(psip1, zp1, s, t, p)
use parameters
use derivs
use ccf
implicit none

integer, intent(in) :: p
double precision, intent(in) :: t, psip1(0:nx,0:nz), zp1(0:nx,0:nz), &
                                s(0:nx)
double precision :: lhs_psi(0:nx,0:nz), psip1_x(0:nx,0:nz), &
                    psip1_xx(0:nx,0:nz), psip1_zz(0:nx,0:nz)
integer :: j, k

call deriv_x(psip1, psip1_x)
call deriv_xx(psip1, psip1_xx)
call deriv_zz(psip1, psip1_zz)

do j = 1, nx-1
   do k = 1, nz-1
      lhs_psi(j,k) = (-1d0 / dx2) * psip1_xx(j,k) + &
                     (0.5d0 * (1d0 - eta) / (s(j) * delx)) * &
                     psip1_x(j,k) - (1d0 / dz2) * &
                     psip1_zz(j,k) - (s(j)**2) * zp1(j,k)
   end do
end do

if (mod(p, save_rate) == 0) then
   write (92, '(2e17.9)') t, lhs_psi(nx/2,3*nz/4)
end if

return
END SUBROUTINE calc_rhs_psi

SUBROUTINE get_new_u(u, uo, po, A, F, s, t)
use parameters
use derivs
implicit none

integer :: j, k
double precision, intent(in) :: uo(0:nx,0:nz), po(0:nx,0:nz), &
                                s(0:nx), A, F(0:nx), t
double precision, intent(out) :: u(0:nx,0:nz)
double precision :: uo_x(0:nx,0:nz), uo_z(0:nx,0:nz), &
                    po_x(0:nx,0:nz), po_z(0:nx,0:nz), &
                    uo_0x(0:nx,0:nz), uo_1x(0:nx,0:nz), &
                    po_0x(0:nx,0:nz), po_1x(0:nx,0:nz), &
                    uo_0z(0:nx,0:nz), uo_1z(0:nx,0:nz), &
                    po_0z(0:nx,0:nz), po_1z(0:nx,0:nz), &
                    uo_xx(0:nx,0:nz), uo_zz(0:nx,0:nz), &
                    uo_0xx(0:nx,0:nz), uo_1xx(0:nx,0:nz), &
                    uo_0zz(0:nx,0:nz), uo_1zz(0:nx,0:nz)

call deriv_x(uo, uo_x, uo_0x, uo_1x)
call deriv_x(po, po_x, po_0x, po_1x)
call deriv_z(uo, uo_z, uo_0z, uo_1z)
call deriv_z(po, po_z, po_0z, po_1z)
call deriv_xx(uo, uo_xx, uo_0xx, uo_1xx)
call deriv_zz(uo, uo_zz, uo_0zz, uo_1zz)

!po_0z(1:nx-1,0) = 4d0 * po(1:nx-1,1) - po(1:nx-1,2)
!po_1z(1:nx-1,nz) = po(1:nx-1,nz-2) - 4d0 * po(1:nx-1,nz-1)

call u_BCS(u, t)

!do j = 1, nx1
!   u(j,:) = rxx * uo_xx(j,:) + &
!            (0.5d0 * (1d0 - eta) * rx * uo_x(j,:) / s(j)) - &
!            (1d0 - eta)**2 * dt * uo(j,:) / s(j)**2 + &
!            rzz * uo_zz(j,:) + uo(j,:) - &
!            0.25d0 * rx * (po_x(j,:) * uo_z(j,:) - po_z(j,:) * &
!            uo_x(j,:)) / (s(j) * delz) + (0.5d0 * (1d0 - eta) * rz / &
!            s(j)) * &
!            (uo(j,:) / s(j) + 2d0 * A) * po_z(j,:) - F(j)
!end do

!do j = 1, nx1
!   u(j,0) = rxx * uo_0xx(j,0) + &
!            (0.5d0 * (1d0 - eta) * rx * uo_0x(j,0) / s(j)) - &
!            (1d0 - eta)**2 * dt * uo(j,0) / s(j)**2 + &
!            rzz * uo_0zz(j,0) + uo(j,0) - &
!            0.25d0 * rx * (po_0x(j,0) * uo_0z(j,0) - po_0z(j,0) * &
!            uo_0x(j,0)) / (s(j) * delz) + (0.5d0 * (1d0 - eta) * rz / &
!            s(j)) * &
!            (uo(j,0) / s(j) + 2d0 * A) * po_0z(j,0) - F(j)
!end do

!do j = 1, nx1
!   u(j,nz) = rxx * uo_1xx(j,nz) + &
!            (0.5d0 * (1d0 - eta) * rx * uo_1x(j,nz) / s(j)) - &
!            (1d0 - eta)**2 * dt * uo(j,nz) / s(j)**2 + &
!            rzz * uo_1zz(j,nz) + uo(j,nz) - &
!            0.25d0 * rx * (po_1x(j,nz) * uo_1z(j,nz) - po_1z(j,nz) * &
!            uo_1x(j,nz)) / (s(j) * delz) + (0.5d0 * (1d0 - eta) * rz / &
!            s(j)) * &
!            (uo(j,nz) / s(j) + 2d0 * A) * po_1z(j,nz) - F(j)
!end do

do k = 1, nz1
do j = 1, nx1
   u(j,k) = rxx * uo_xx(j,k) + &
            (0.5d0 * rx * (1d0 - eta) / s(j)) * uo_x(j,k) - &
            ((1d0 - eta)**2 * dt / s(j)**2) * uo(j,k) + &
            rzz * uo_zz(j,k) + uo(j,k) - &
            ((1d0 * 0.25d0 * rx) / (s(j) * delz)) * &
            (po_x(j,k) * uo_z(j,k) - po_z(j,k) * uo_x(j,k)) + &
            (0.5d0 * (1d0 - eta) * rz / s(j)) * &
            ((1d0 * uo(j,k) / s(j)) + 0d0 * A) * &
            po_z(j,k) - F(j)
end do
end do

do j = 1, nx-1
  u(j,0) = rxx * uo_0xx(j,0) + &
            (0.5d0 * rx * (1d0 - eta) / s(j)) * uo_0x(j,0) - &
            ((1d0 - eta)**2 * dt / s(j)**2) * uo(j,0) + &
            rzz * uo_0zz(j,0) + uo(j,0) - &
            ((1d0 * 0.25d0 * rx) / (s(j) * delz)) * &
            (po_0x(j,0) * uo_0z(j,0) - po_0z(j,0) * uo_0x(j,0)) + &
            (0.5d0 * (1d0 - eta) * rz / s(j)) * &
            ((1d0 * uo(j,0) / s(j)) + 0d0 * A) * &
            po_0z(j,0) - F(j)
!u(j,0)=0d0
end do

do j = 1, nx-1
   u(j,nz) = rxx * uo_1xx(j,nz) + &
            (0.5d0 * rx * (1d0 - eta) / s(j)) * uo_1x(j,nz) - &
            ((1d0 - eta)**2 * dt / s(j)**2) * uo(j,nz) + &
            rzz * uo_1zz(j,nz) + uo(j,nz) - &
            ((1d0 * 0.25d0 * rx) / (s(j) * delz)) * &
            (po_1x(j,nz) * uo_1z(j,nz) - po_1z(j,nz) * uo_1x(j,nz)) + &
            (0.5d0 * (1d0 - eta) * rz / s(j)) * &
            ((1d0 * uo(j,nz) / s(j)) + 0d0 * A) * &
            po_1z(j,nz) - F(j)
!u(j,nz)=0d0
end do

return
END SUBROUTINE get_new_u

SUBROUTINE get_new_Z(zn, zo, po, uo, vc, s, t)
use parameters
use derivs
implicit none

integer :: j
double precision, intent(in) :: zo(0:nx,0:nz), uo(0:nx,0:nz), &
                                po(0:nx,0:nz), vc(0:nx), s(0:nx), t
double precision, intent(out) :: zn(0:nx,0:nz)
double precision :: zo_xx(0:nx,0:nz), zo_x(0:nx,0:nz), &
                    zo_zz(0:nx,0:nz), uo_z(0:nx,0:nz), &
                    po_x(0:nx,0:nz), zo_z(0:nx,0:nz), &
                    po_z(0:nx,0:nz)

call deriv_xx(zo, zo_xx)                    
call deriv_x(zo, zo_x)
call deriv_zz(zo, zo_zz)
call deriv_z(uo, uo_z)
call deriv_x(po, po_x)
call deriv_z(zo, zo_z)
call deriv_z(po, po_z)

call z_BCS(zn, po, s, t)

do j = 1, nx-1
   zn(j,:) = rxx * zo_xx(j,:) + &
        (1.5d0 * rx * (1d0 - eta) / s(j)) * zo_x(j,:) + &
        rzz * zo_zz(j,:) + zo(j,:) + &
        ((1d0 - eta) * rz / s(j)**2) * (0d0*vc(j) + 1d0 * uo(j,:)) * &
        uo_z(j,:) - &
        ((1d0 * 0.25d0 * rx) / (s(j) * delz)) * &
        (po_x(j,:) * zo_z(j,:) - po_z(j,:) * zo_x(j,:))
end do
call z_BCS(zn, po, s, t)

return
END SUBROUTINE get_new_Z

SUBROUTINE psi_deriv_test(pn, s)
use parameters
implicit none

double precision, intent(in) :: pn(0:nx,0:nz), s(0:nx)
double precision :: xR1_deriv(0:nz), xR2_deriv(0:nz), zx0(0:nz), &
                    zx1(0:nz), zx0_alt(0:nz), zx1_alt(0:nz)
integer, parameter :: zpos = nz/2-1

xR1_deriv(:) = (-pn(2,:) + 4d0 * pn(1,:) - 3d0 * pn(0,:)) / (2d0 * delx)
xR2_deriv(:) = (pn(nx-2,:) - 4d0 * pn(nx1,:) + 3d0 * pn(nx,:)) / (2d0 * delx)

!write (93, '(2e17.9)') xR1_deriv(zpos), xR2_deriv(zpos)

zx0(:) = (-1d0 / (dx2 * s(0)**2)) * (-pn(3,:) + 4d0 * pn(2,:) - &
          5d0 * pn(1,:) + 2d0 * pn(0,:))+ &
          ((1d0 - eta) / (2d0 * delx * s(0)**3)) * (-pn(2,:) + 4d0 * &
          pn(1,:) - 3d0 * pn(0,:))

zx1(:) = (-1d0 / (dx2 * s(nx)**2)) * (-pn(nx-3,:) + 4d0 * pn(nx-2,:) - &
          5d0 * pn(nx1,:) + 2d0 * pn(nx,:))+ &
          ((1d0 - eta) / (2d0 * delx * s(nx)**3)) * (pn(nx-2,:) - 4d0 * &
          pn(nx1,:) + 3d0 * pn(nx,:))

zx0_alt(:) = -(8d0 * pn(1,:) - pn(2,:)) / (2d0 * (eta**2) * dx2)
zx1_alt(:) = -(8d0 * pn(nx1,:) - pn(nx-2,:)) / (2d0 * dx2)

write (93, '(4e17.9)') zx0(zpos), zx0_alt(zpos), zx1(zpos), zx1_alt(zpos)


return
END SUBROUTINE psi_deriv_test
