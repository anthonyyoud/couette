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
double precision :: x(0:nx), z(0:nz), s(0:nx), &
!double precision :: s(0:nx), &
                    t = 0d0, A = 0d0, A_ = 0d0, B = 0d0, B_ = 0d0, &
u_nlin_new(0:nx,0:nz), u_nlin_old(0:nx,0:nz), &
z_nlin_new(0:nx,0:nz), z_nlin_old(0:nx,0:nz), &
unew(0:nx,0:nz), uold(0:nx,0:nz), uold2(0:nx,0:nz), u_int(0:nx,0:nz), &
znew(0:nx,0:nz), zold(0:nx,0:nz), zold2(0:nx,0:nz), z_int(0:nx,0:nz), &
pnew(0:nx,0:nz), pold(0:nx,0:nz), pold2(0:nx,0:nz), &
vr(0:nx,0:nz), vz(0:nx,0:nz), &
vc(0:nx), vc_(0:nx), vr2(0:nx,0:nz) = 0d0, &
AB(2*nx1+nx1+1,nx1*nz1), F(0:nx)
integer :: pivot(nx1*nz1)

logical, parameter :: write_ofile = .true.
integer :: j, k, p = 0, p_start = 1

print*
if (tau == 0) then
   write(6, '(A7, i2, A21)') 'tau = ', tau, '- Infinite cylinder'
else if (tau == 1) then
   write(6, '(A7, i2, A22)') 'tau = ', tau, '- Finite aspect ratio'
end if   

call open_files()

print*, 'Setting up ICS...'
call ICS(unew, znew, pnew, x, z, s)

if (.not. restart) then
   print*, 'Setting up BCS...'
   call u_BCS(unew, vc, 0d0)
   call p_BCS(pnew)
   call z_BCS(znew, pnew, s, 0d0)
end if

print*, 'Setting up matrices...'
call matrix_setup(Ux, Uz, Zx, Zz, s)
call ABC_mat_setup(AB, pivot, s)

uold = unew
uold2 = unew
u_int = unew
zold = znew
zold2 = znew
z_int = znew
pold = pnew
pold2 = pnew

print*, 'Entering time loop'

if (restart) then
   call get_p(p_start)
end if

!call get_timestep()

do p = p_start, Ntot
   inquire(file='RUNNING', exist=file_exist)
   file_exist2 = .not. file_exist
   if (file_exist2) then
      print*, 'Stop requested'
      print*, 'Saving end state'
      call end_state(uold, zold, pold, p)
      call save_xsect(vr, vz, x, z, p)
      call save_surface(pold, uold, zold, vr, vz, x, z, p, t)
      exit
   end if

   t = p * dt

   call get_A(t, A, A_)
   call get_B(t, B, B_)
   call get_vc(A, A_, B_, B_, vc, vc_, s)
   call get_F(t, s, F)

   uold = unew
   zold = znew

!*** Get RHS for v in x-direction ************

   call get_rhs_ux(uold, unew, s)

!*********************************************

!*** Get non-lin terms for v in x-direction **

   call get_nlin_ux(s, A, A_, F, &
                    uold, uold2, pold, pold2, u_nlin_new)

!*********************************************

!*** Get RHS for Z in x-direction ************

   call get_rhs_Zx(zold, znew, s)

!*********************************************

!*** Get non-lin terms for Z in x-direction **

   call get_nlin_Zx(s, t, A, A_, vc, vc_, &
                    uold, uold2, pold, pold2, zold, zold2, z_nlin_new)

!*********************************************

!*** Solve for v in x-direction **************

   call solve_ux(uold, unew, u_nlin_new, s, vc, t, Ux)

!*********************************************

!*** Solve for Z in x-direction **************

   call solve_Zx(zold, znew, z_nlin_new, pold, s, t, Zx)

!*********************************************

   uold2 = u_int
   zold2 = z_int

!*** Solve for v in z-direction **************

   call solve_uz(uold, unew, vc, t, Uz)

!*********************************************

!*** Solve for Z in z-direction **************

   call solve_Zz(zold, pold, znew, s, t, Zz)

!*********************************************

   u_int = unew
   z_int = znew

   pold2 = pold
   call poisson(znew, pnew, AB, pivot, s)

   if (diag) then
      call calc_rhs_u(unew, uold2, pold, s, t, A, p)
      call calc_rhs_Z(znew, zold2, uold2, pold, s, t, A, vc, p)
      call calc_rhs_psi(pnew, znew, s, t, p)
   end if

   pold = pnew

   if (mod(p, save_rate) == 0) then
      call r_vel(pold, s, vr, vz)
      call save_torque(t, unew)
      if (p /= save_rate) then
         call save_growth(t, vr, vr2, vz, pold, unew, znew)
      end if
   end if

   vr2 = vr

   if (xsect_save) then
      if (mod(p, save_rate_2) == 0) then
         call save_xsect(vr, vz, x, z, p)
         call save_surface(pold, unew, znew, vr, vz, x, z, p, t)
      end if
   end if

   if (p == Ntot) then
      call end_state(unew, znew, pold, p)
      call save_xsect(vr, vz, x, z, p)
      call save_surface(pold, unew, znew, vr, vz, x, z, p, t)
   end if

end do
!*** END TIME LOOP***

call close_files()

END PROGRAM couette_mod

SUBROUTINE ICS(u, zn, pn, x, z, s)
use parameters
use io
implicit none
double precision, intent(inout) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                   pn(0:nx,0:nz)
double precision, intent(out) :: x(0:nx), z(0:nz), s(0:nx)
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
   call state_restart(u, zn, pn)

else
   if (tau == 1) then
      do k = 0, nz
         u(:,k) = seed * dsin(pi*z(k)/gamma) * dsin(pi*x(:))
      end do
      do j = 0, nx
         pn(j,:) = seed * dsin(2d0*pi*z(:)/gamma) * dsin(pi*x(j))
      end do
      zn(:,:) = seed
   else
      u(:,:nz/2) = seed
      u(:,nz/2+1:nz) = seed
      pn(:,:nz/2) = -seed 
      pn(:,nz/2+1:nz) = seed
   end if
end if

return
END SUBROUTINE ICS

SUBROUTINE u_BCS(u, vc, t)
use parameters
implicit none
double precision, intent(out) :: u(0:nx,0:nz)
double precision, intent(in) :: vc(0:nx), t
integer :: k

u(0,:) = 0d0
u(nx,:) = 0d0

if (tau == 1) then
   u(:,0) = -vc(:)
   u(:,nz) = -vc(:)
end if

return
END SUBROUTINE u_BCS

SUBROUTINE z_BCS(zn, pn, s, t)
use parameters
implicit none
double precision, intent(out) :: zn(0:nx,0:nz)
double precision, intent(in) :: t, pn(0:nx,0:nz), s(0:nx)

zn(0,:) = -(8d0 * pn(1,:) - pn(2,:)) / (2d0 * (eta**2) * dx2)
zn(nx,:) = -(8d0 * pn(nx1,:) - pn(nx-2,:)) / (2d0 * dx2)

if (tau == 1) then
   zn(:,0) = -(8d0 * pn(:,1) - pn(:,2)) / &
              (2d0 * (s(:)**2) * dz2)  
   zn(:,nz) = -(8d0 * pn(:,nz1) - pn(:,nz-2)) / &
               (2d0 * (s(:)**2) * dz2)
else
   zn(:,0) = 0d0 
   zn(:,nz) = 0d0 
end if

return
END SUBROUTINE z_BCS

SUBROUTINE get_rhs_ux(uo, u, s)
use parameters
use derivs
implicit none

double precision, intent(in) :: s(0:nx), uo(0:nx,0:nz)
double precision, intent(out) :: u(0:nx,0:nz)
double precision :: uo_x(0:nx,0:nz), uo_0x(0:nx,0:nz), &
                    uo_1x(0:nx,0:nz),  uo_xx(0:nx,0:nz), &
                    uo_0xx(0:nx,0:nz),  uo_1xx(0:nx,0:nz), &
                    uo_zz(0:nx,0:nz), uo_0zz(0:nx,0:nz), &
                    uo_1zz(0:nx,0:nz)
integer :: j, k

call deriv_x(uo, uo_x, uo_0x, uo_1x)
call deriv_xx(uo, uo_xx, uo_0xx, uo_1xx)
call deriv_zz(uo, uo_zz, uo_0zz, uo_1zz)

do j = 1, nx1
   u(j,1:nz1) = uo(j,1:nz1) + (0.5d0 * rxx * uo_xx(j,1:nz1)) + &
            (((1d0 - eta) * rx) / (4d0 * s(j))) * uo_x(j,1:nz1) - &
            (((1d0 - eta)**2 * dt) / (2d0 * s(j)**2)) * &
            uo(j,1:nz1) + 0.5d0 * rzz * uo_zz(j,1:nz1)
end do

if (tau /= 1) then
   do j = 1, nx1
      u(j,0) = uo(j,0) + (0.5d0 * rxx * uo_0xx(j,0)) + &
               (((1d0 - eta) * rx) / (4d0 * s(j))) * uo_0x(j,0) - &
               (((1d0 - eta)**2 * dt) / (2d0 * s(j)**2)) * &
               uo(j,0) + 0.5d0 * rzz * uo_0zz(j,0)

      u(j,nz) = uo(j,nz) + (0.5d0 * rxx * uo_1xx(j,nz)) + &
               (((1d0 - eta) * rx) / (4d0 * s(j))) * uo_1x(j,nz) - &
               (((1d0 - eta)**2 * dt) / (2d0 * s(j)**2)) * &
               uo(j,nz) + 0.5d0 * rzz * uo_1zz(j,nz)
   end do
end if

return
END SUBROUTINE get_rhs_ux

SUBROUTINE get_nlin_ux(s, A, A_, F, &
                       uo, uo2, po, po2, u_nl_n)
use parameters
use derivs
use ccf
implicit none
double precision, intent(in) :: A, A_, F(0:nx), &
                                uo(0:nx,0:nz), uo2(0:nx,0:nz), &
                                po(0:nx,0:nz), po2(0:nx,0:nz), &
                                s(0:nx)
double precision, intent(out) :: u_nl_n(0:nx,0:nz)
double precision :: uo_x(0:nx,0:nz), uo2_x(0:nx,0:nz), &
                    uo_z(0:nx,0:nz), uo2_z(0:nx,0:nz), &
                    po_x(0:nx,0:nz), po2_x(0:nx,0:nz), &
                    po_z(0:nx,0:nz), po2_z(0:nx,0:nz), &
                    uo_0x(0:nx,0:nz), uo2_0x(0:nx,0:nz), &
                    uo_1x(0:nx,0:nz), uo2_1x(0:nx,0:nz), &
                    po_0z(0:nx,0:nz), po2_0z(0:nx,0:nz), &
                    po_1z(0:nx,0:nz), po2_1z(0:nx,0:nz)
                   
integer :: j, k

call deriv_x(uo, uo_x, uo_0x, uo_1x)
call deriv_x(uo2, uo2_x, uo2_0x, uo2_1x)
call deriv_z(uo, uo_z)
call deriv_z(uo2, uo2_z)
call deriv_x(po, po_x)
call deriv_x(po2, po2_x)
call deriv_z(po, po_z, po_0z, po_1z)
call deriv_z(po2, po2_z, po2_0z, po2_1z)

do j = 1, nx1
   u_nl_n(j,1:nz1) = 1d0*( (-rx / (8d0 * s(j) * delz)) * &
                 (3d0 * (po_x(j,1:nz1) * uo_z(j,1:nz1) - &
                 po_z(j,1:nz1) * uo_x(j,1:nz1)) - &
                 (po2_x(j,1:nz1) * uo2_z(j,1:nz1) - &
                 po2_z(j,1:nz1) * uo2_x(j,1:nz1))) + &
                 ((1d0 - eta) * rz / (4d0 * s(j)**2)) * &
                 (3d0 * uo(j,1:nz1) * po_z(j,1:nz1) - &
                 uo2(j,1:nz1) * po2_z(j,1:nz1)) )+ &
                 ((1d0 - eta) * rz / (2d0 * s(j))) * &
                 (3d0 * A * po_z(j,1:nz1) - &
                 A_ * po2_z(j,1:nz1)) - F(j)
end do

if (tau /= 1) then
   do j = 1, nx1
      u_nl_n(j,0) = 1d0*( (-rx / (8d0 * s(j) * delz)) * &
                 (-3d0 * po_0z(j,0) * uo_0x(j,0) + &
                 po2_0z(j,0) * uo2_0x(j,0)) + &
                 ((1d0 - eta) * rz / (4d0 * s(j)**2)) * &
                 (3d0 * uo(j,0) * po_0z(j,0) - &
                 uo2(j,0) * po2_0z(j,0)) )+ &
                 ((1d0 - eta) * rz / (2d0 * s(j))) * &
                 (3d0 * A * po_0z(j,0) - &
                 A_ * po2_0z(j,0)) - F(j)

   u_nl_n(j,nz) = 1d0*( (-rx / (8d0 * s(j) * delz)) * &
                  (-3d0 * po_1z(j,nz) * uo_1x(j,nz) + &
                  po2_1z(j,nz) * uo2_1x(j,nz)) + &
                  ((1d0 - eta) * rz / (4d0 * s(j)**2)) * &
                  (3d0 * uo(j,nz) * po_1z(j,nz) - &
                  uo2(j,nz) * po2_1z(j,nz)) )+ &
                  ((1d0 - eta) * rz / (2d0 * s(j))) * &
                  (3d0 * A * po_1z(j,nz) - &
                  A_ * po2_1z(j,nz)) - F(j)
   end do
end if

return
END SUBROUTINE get_nlin_ux

SUBROUTINE get_rhs_Zx(zo, zn, s)
use parameters
use derivs
implicit none
double precision, intent(in) :: s(0:nx), zo(0:nx,0:nz)
double precision, intent(out) :: zn(0:nx,0:nz)
double precision :: zo_x(0:nx,0:nz), zo_xx(0:nx,0:nz), zo_zz(0:nx,0:nz)
integer :: j, k

call deriv_x(zo, zo_x)
call deriv_xx(zo, zo_xx)
call deriv_zz(zo, zo_zz)

do j = 1, nx1
   zn(j,1:nz1) = zo(j,1:nz1) + (0.5d0 * rxx * zo_xx(j,1:nz1)) + &
                  ((3d0 * (1d0 - eta) * rx) / (4d0 * s(j))) * &
                  zo_x(j,1:nz1) + 0.5d0 * rzz * zo_zz(j,1:nz1)
end do

return
END SUBROUTINE get_rhs_Zx

SUBROUTINE get_nlin_Zx(s, t, A, A_, vc, vc_, &
                       uo, uo2, po, po2, zo, zo2, z_nl_n)
use parameters
use derivs
use ccf
implicit none
double precision, intent(in) :: t, A, A_,&
                                uo(0:nx,0:nz), uo2(0:nx,0:nz), &
                                po(0:nx,0:nz), po2(0:nx,0:nz), &
                                zo(0:nx,0:nz), zo2(0:nx,0:nz), &
                                s(0:nx), vc(0:nx), vc_(0:nx)
double precision, intent(out) :: z_nl_n(0:nx,0:nz)
double precision :: uo_z(0:nx,0:nz), uo2_z(0:nx,0:nz), &
                    po_x(0:nx,0:nz), po2_x(0:nx,0:nz), &
                    po_z(0:nx,0:nz), po2_z(0:nx,0:nz), &
                    zo_x(0:nx,0:nz), zo2_x(0:nx,0:nz), &
                    zo_z(0:nx,0:nz), zo2_z(0:nx,0:nz)

integer :: j, k

call deriv_z(uo, uo_z)
call deriv_z(uo2, uo2_z)
call deriv_x(po, po_x)
call deriv_x(po2, po2_x)
call deriv_z(po, po_z)
call deriv_z(po2, po2_z)
call deriv_x(zo, zo_x)
call deriv_x(zo2, zo2_x)
call deriv_z(zo, zo_z)
call deriv_z(zo2, zo2_z)

do j = 1, nx1
   z_nl_n(j,1:nz1) = (((1d0 - eta) * rz) / (2d0 * s(j)**2)) * &
                 (3d0 * (vc(j) + 1d0*uo(j,1:nz1)) * uo_z(j,1:nz1) - &
                 (vc_(j) + 1d0*uo2(j,1:nz1)) * uo2_z(j,1:nz1)) - &
                 (rx / (8d0 * s(j) * delz)) * &
                 1d0*((3d0 * (po_x(j,1:nz1) * zo_z(j,1:nz1) - &
                 po_z(j,1:nz1) * zo_x(j,1:nz1))) - &
                 (po2_x(j,1:nz1) * zo2_z(j,1:nz1) - &
                 po2_z(j,1:nz1) * zo2_x(j,1:nz1)))
end do

return
END SUBROUTINE get_nlin_Zx

SUBROUTINE solve_ux(uo, u, u_nl, s, vc, t, ux)
use parameters
use io
implicit none
double precision, intent(in) :: t
double precision, intent(in) :: u(0:nx,0:nz), u_nl(0:nx,0:nz), &
                                s(0:nx), vc(0:nx)
type (mat_comp), intent(in) :: ux
double precision, intent(inout) :: uo(0:nx,0:nz)
double precision :: ux_rhs(nx1)
integer :: j, k

call u_BCS(uo, vc, t)

do k = 1, nz1
   ux_rhs(:) = u(1:nx1,k) + u_nl(1:nx1,k)
   
   ux_rhs(1) = ux_rhs(1) + (0.5d0 * rxx * uo(0,k)) - &
               (((1d0 - eta) * rx) / (4d0 * s(1))) * uo(0,k)
   ux_rhs(nx1) = ux_rhs(nx1) + (0.5d0 * rxx * uo(nx,k)) + &
               (((1d0 - eta) * rx) / (4d0 * s(nx1))) * uo(nx,k)

   call thomas(xlb, nx1, ux%up, ux%di, ux%lo, ux_rhs)

   uo(1:nx1,k) = ux_rhs(:)
end do

if (tau /= 1) then
   do k = 0, nz, nz
      ux_rhs(:) = u(1:nx1,k) + u_nl(1:nx1,k)
                                
      ux_rhs(1) = ux_rhs(1) + (0.5d0 * rxx * uo(0,k)) - &
                  (((1d0 - eta) * rx) / (4d0 * s(1))) * uo(0,k)
      ux_rhs(nx1) = ux_rhs(nx1) + (0.5d0 * rxx * uo(nx,k)) + &
                  (((1d0 - eta) * rx) / (4d0 * s(nx1))) * uo(nx,k)

      call thomas(xlb, nx1, ux%up, ux%di, ux%lo, ux_rhs)

      uo(1:nx1,k) = ux_rhs(:)
   end do
end if

return
END SUBROUTINE solve_ux

SUBROUTINE solve_Zx(zo, zn, z_nl, po, s, t, zx)
use parameters
use io
implicit none
double precision, intent(in) :: t
double precision, intent(in) :: zn(0:nx,0:nz), po(0:nx,0:nz), &
                                z_nl(0:nx,0:nz), s(0:nx)
type (mat_comp), intent(in) :: zx
double precision, intent(inout) :: zo(0:nx,0:nz)
double precision :: zx_rhs(nx1)
integer :: j, k

call z_BCS(zo, po, s, t)

do k = 1, nz1
   zx_rhs(:) = zn(1:nx1,k) + z_nl(1:nx1,k)
   
   zx_rhs(1) = zx_rhs(1) + (0.5d0 * rxx * zo(0,k)) - &
               ((3d0 * (1d0 - eta) * rx) / (4d0 * s(1))) * zo(0,k)
   zx_rhs(nx1) = zx_rhs(nx1) + (0.5d0 * rxx * zo(nx,k)) + &
               ((3d0 * (1d0 - eta) * rx) / (4d0 * s(nx1))) * zo(nx,k)

   call thomas(xlb, nx1, zx%up, zx%di, zx%lo, zx_rhs)

   zo(1:nx1,k) = zx_rhs(:)
end do

return
END SUBROUTINE solve_Zx

SUBROUTINE solve_uz(uo, u, vc, t, uz)
use parameters
use io
implicit none
double precision, intent(in) :: t
double precision, intent(in) :: uo(0:nx,0:nz), vc(0:nx)
type (uz_mat_comp), intent(in) :: uz
double precision, intent(inout) :: u(0:nx,0:nz)
double precision :: uz_rhs(0:nz), uz_rhs_t1(nz1), &
                    up(nz-2), di(nz1), lo(2:nz1)
integer :: j, k

if (tau == 1) then
   call u_BCS(u, vc, t)
   up(:) = uz%up(1:nz-2)
   di(:) = uz%di(1:nz1)
   lo(:) = uz%lo(2:nz1)

   do j = 1, nx1
      uz_rhs_t1(:) = uo(j,1:nz1)

      uz_rhs_t1(1) = uz_rhs_t1(1) + 0.5d0 * rzz * u(j,0)
      uz_rhs_t1(nz1) = uz_rhs_t1(nz1) + 0.5d0 * rzz * u(j,nz)

      call thomas(zlb+1, nz1, up, di, lo, uz_rhs_t1)

      u(j,1:nz1) = uz_rhs_t1(:)
   end do
else
   do j = 1, nx1
      uz_rhs(:) = uo(j,:)
   
      call thomas(zlb, nz, uz%up, uz%di, uz%lo, uz_rhs)

      u(j,:) = uz_rhs(:)
   end do
end if

return
END SUBROUTINE solve_uz

SUBROUTINE solve_Zz(zo, po, zn, s, t, zz)
use parameters
use io
implicit none
double precision, intent(in) :: t
double precision, intent(in) :: zo(0:nx,0:nz), po(0:nx,0:nz), s(0:nx)
type (zz_mat_comp), intent(in) :: zz
double precision, intent(inout) :: zn(0:nx,0:nz)
double precision :: Zz_rhs(nz1)
integer :: j, k

call z_BCS(zn, po, s, t)

do j = 1, nx1
   Zz_rhs(:) = zo(j,1:nz1)
   
   Zz_rhs(1) = Zz_rhs(1) + 0.5d0 * rzz * zn(j,0)
   Zz_rhs(nz1) = Zz_rhs(nz1) + 0.5d0 * rzz * zn(j,nz)

   call thomas(zlb+1, nz1, zz%up, zz%di, zz%lo, Zz_rhs)

   zn(j,1:nz1) = Zz_rhs(:)
end do

return
END SUBROUTINE solve_Zz

SUBROUTINE r_vel(p, s, vr, vz)
use parameters
implicit none
double precision, intent(in) :: p(0:nx,0:nz), s(0:nx)
double precision, intent(out) :: vr(0:nx,0:nz), vz(0:nx,0:nz)
integer :: j, k

do k = 1, nz1
   do j = 0, nx
      vr(j,k) = (-1d0 / (2d0 * s(j) * delz)) * (p(j,k+1) - p(j,k-1))
   end do
end do

do j = 0, nx
   vr(j,0) = (-1d0 / (2d0 * s(j) * delz)) * &
              (-3d0 * p(j,0) + 4d0 * p(j,1) - p(j,2))
   vr(j,nz) = (-1d0 / (2d0 * s(j) * delz)) * &
              (3d0 * p(j,nz) - 4d0 * p(j,nz1) + p(j,nz-2))
end do

do k = 0, nz
   do j = 1, nx1
      vz(j,k) = (1d0 / (2d0 * s(j) * delx)) * (p(j+1,k) - p(j-1,k))
   end do
end do

do k = 0, nz
   vz(0,k) = 0d0 !(1d0 / (2d0 * s(0) * delx)) * &
              !(-3d0 * p(0,k) + 4d0 * p(1,k) - p(2,k))
   vz(nx,k) = 0d0 !(1d0 / (2d0 * s(nx) * delx)) * &
              !(3d0 * p(nx,k) - 4d0 * p(nx1,k) + p(nx-2,k))
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

do j = 1, nx1
   do k = 1, nz1
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

do j = 1, nx1
   do k = 1, nz1
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

do j = 1, nx1
   do k = 1, nz1
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
