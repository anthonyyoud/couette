MODULE ic_bc
use parameters
implicit none

double precision :: x(0:nx), x_(0:nx), th(0:nt), z(0:nz), s(0:nx) 
                                                !finite-difference mesh
                                                !s=eta+(1-eta)*x
contains

SUBROUTINE get_xzs()
!Finite-difference mesh
use parameters
implicit none

integer :: j, k, l

do k = 0, nz
   z(k) = dble(k) * delz
end do
do l = 0, nt
   th(l) = dble(l) * delt
end do
do j = 0, nx
   x(j) = dble(j) * delx
end do

x_ = x + 1d0
s = eta + ((1d0 - eta) * x)

return
END SUBROUTINE get_xzs

SUBROUTINE ICS(u, zn, pn, bn, jn, p)
!Initial conditions
use parameters
implicit none
double precision, intent(inout) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                   pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                   jn(0:nx,0:nz)
integer, intent(out) :: p
real :: rand
integer :: j, k

if (restart) then
   print*, 'Getting restart conditions'
   call state_restart(u, zn, pn, bn, jn, p)  !get saved data if restart
else                   !put in estimate of eigen-function shape
   if (tau == 1) then  !which satisfies BCS multiplied by small seed
      do k = 0, nz
         u(:,k) = seed * dsin(2d0*pi*z(k)/gamma) * dsin(pi*x(:))
         pn(:,k) = seed * dsin(2d0*pi*z(k)/gamma) * dsin(pi*x(:))
         bn(:,k) = seed * dcos(2d0*pi*z(k)/gamma) / s(:)
         jn(:,k) = seed * dsin(pi*x(:)) * dsin(2d0*pi*z(k)/gamma)
      end do
   else
      do k = 0, nz
         u(:,k) = seed * dsin(pi*x(:)) * dcos(2d0*pi*z(k)/gamma)
         pn(:,k) = seed * dsin(pi*x(:)) * dsin(2d0*pi*z(k)/gamma)
         bn(:,k) = seed * dsin(2d0*pi*z(k)/gamma) / s(:)
         jn(:,k) = seed * dsin(pi*x(:)) * dcos(2d0*pi*z(k)/gamma)
      end do

      do k = 1, nz1
        do j = 1, nx1
           zn(j,k) = -(pn(j+1,k) - 2d0 * pn(j,k) + pn(j-1,k)) / &
                      (s(j) * dx2) + &
                      0.5d0 * (1d0 - eta) * (pn(j+1,k) - pn(j-1,k)) / &
                      (s(j)**2 * delx) - &
                      (pn(j,k+1) - 2d0 * pn(j,k) + pn(j,k-1)) / &
                      (s(j) * dz2)  !ICS based on above fields for vorticity
        end do
      end do
   end if
end if

return
END SUBROUTINE ICS

SUBROUTINE state_restart(u, zn, pn, bn, jn, p)
!Get restart data
use parameters
implicit none
double precision, intent(out) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                 pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                 jn(0:nx,0:nz)
double precision :: dt_prev
integer, intent(out) :: p
integer :: j, k

open (50, file = 'end_state.dat')

read(50, *) p
read(50, *) dt_prev
read(50, *) ((u(j,k), k = 0, nz), j = 0, nx)
read(50, *) ((zn(j,k), k = 0, nz), j = 0, nx)
read(50, *) ((pn(j,k), k = 0, nz), j = 0, nx)
read(50, *) ((bn(j,k), k = 0, nz), j = 0, nx)
read(50, *) ((jn(j,k), k = 0, nz), j = 0, nx)

close (50)   

!u = u
!zn = zn
!pn = pn
!bn = bn
!jn = jn

if (dt_prev /= dt) then
   p = p * dt_prev / dt
end if

return
END SUBROUTINE state_restart

SUBROUTINE u_BCS(u, t)
!Boundary conditions for total azimuthal velocity (including CCF)
use parameters
implicit none
double precision, intent(out) :: u(0:nx,0:nz)
double precision, intent(in) :: t
integer :: k

u(0,:) = Re1 + Re1_mod * dcos(om1 * t) + &
         eps1 * Re1 * (1d0 / eta - 1d0) * dcos(freq1 * z(:))
u(nx,:) = Re2 + Re2_mod * dcos(om2 * t) - &
         eps2 * dcos(freq2 * z(:))

if (tau == 1) then
   u(:,0) = 0d0
   u(:,nz) = 0d0
end if

return
END SUBROUTINE u_BCS

SUBROUTINE z_BCS(zn, pn, t)
!Boundary conditions for azimuthal vorticity
use parameters
implicit none
double precision, intent(out) :: zn(0:nx,0:nz)
double precision, intent(in) :: t, pn(0:nx,0:nz)

zn(0,:) = -(8d0 * pn(1,:) - pn(2,:)) / (2d0 * s(0) * dx2)
zn(nx,:) = -(8d0 * pn(nx1,:) - pn(nx-2,:)) / (2d0 * s(nx) * dx2)

if (tau == 1) then
   zn(:,0) = -(8d0 * pn(:,1) - pn(:,2)) / &
              (2d0 * (s(:)) * dz2)
   zn(:,nz) = -(8d0 * pn(:,nz1) - pn(:,nz-2)) / &
               (2d0 * (s(:)) * dz2)
else
   zn(:,0) = (-tau / (s(:) * (1d0 - tau))) * &
             (0.5d0 * (-pn(:,2) + 4d0 * pn(:,1)) / delz)
   zn(:,nz) = (tau / (s(:) * (1d0 - tau))) * &
             (0.5d0 * (pn(:,nz-2) - 4d0 * pn(:,nz1)) / delz)
end if

return
END SUBROUTINE z_BCS

SUBROUTINE p_BCS(p)
!Boundary conditions for stream-function, psi
use parameters
implicit none
double precision, intent(out) :: p(0:nx,0:nz)

p(0,:) = 0d0
p(nx,:) = 0d0

p(:,0) = 0d0
p(:,nz) = 0d0

return
END SUBROUTINE p_BCS

SUBROUTINE b_BCS(bn)
!Boundary conditions for azimuthal magnetic field
use parameters
implicit none
double precision, intent(out) :: bn(0:nx,0:nz)

if (tau == 0) then
   bn(:,0) = 0d0
   bn(:,nz) = 0d0
end if

return
END SUBROUTINE b_BCS

SUBROUTINE j_BCS(jn)
!Boundary conditions for azimuthal current
use parameters
implicit none
double precision, intent(out) :: jn(0:nx,0:nz)

jn(0,:) = 0d0
jn(nx,:) = 0d0

if (tau == 1) then
   jn(:,0) = 0d0
   jn(:,nz) = 0d0
end if

return
END SUBROUTINE j_BCS

FUNCTION f1(index)
use parameters
implicit none

integer, intent(in) :: index
double precision :: f1

f1 = eps1 * dcos(freq1*z(index))

return
END FUNCTION f1   

FUNCTION f2(index)
use parameters
implicit none

integer, intent(in) :: index
double precision :: f2

f2 = eps2 * dcos(freq2*z(index)-pi)

return
END FUNCTION f2 

END MODULE ic_bc
