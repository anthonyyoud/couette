MODULE linear
implicit none

contains

SUBROUTINE get_rhs_ux(uo, u)
use parameters
use variables
use derivs
use ic_bc
implicit none

double precision, intent(in) :: uo(0:nx,0:nz)
double precision, intent(out) :: u(0:nx,0:nz)
type (deriv) :: du
integer :: j, k

call deriv_x(uo, du%x)
call deriv_xx(uo, du%xx)
call deriv_zz(uo, du%zz)

do k = 1, nz1
   u(1:nx1,k) = uo(1:nx1,k) + (0.5d0 * rxx * du%xx(1:nx1,k)) + &
            (((1d0 - eta) * rx) / (4d0 * s(1:nx1))) * du%x(1:nx1,k) - &
            (((1d0 - eta)**2 * dt) / (2d0 * s(1:nx1)**2)) * &
            uo(1:nx1,k) + 0.5d0 * rzz * du%zz(1:nx1,k)
end do

if (tau /= 1) then
   u(1:nx1,0) = uo(1:nx1,0) + (0.5d0 * rxx * du%xx(1:nx1,0)) + &
                (((1d0 - eta) * rx) / (4d0 * s(1:nx1))) * du%x(1:nx1,0) - &
                (((1d0 - eta)**2 * dt) / (2d0 * s(1:nx1)**2)) * &
                uo(1:nx1,0) + 0.5d0 * rzz * du%zz(1:nx1,0)

   u(1:nx1,nz) = uo(1:nx1,nz) + (0.5d0 * rxx * du%xx(1:nx1,nz)) + &
                 (((1d0 - eta) * rx) / (4d0 * s(1:nx1))) * du%x(1:nx1,nz) - &
                 (((1d0 - eta)**2 * dt) / (2d0 * s(1:nx1)**2)) * &
                 uo(1:nx1,nz) + 0.5d0 * rzz * du%zz(1:nx1,nz)
end if

return
END SUBROUTINE get_rhs_ux

SUBROUTINE get_rhs_Zx(zo, zn)
use parameters
use variables
use derivs
use ic_bc
implicit none
double precision, intent(in) :: zo(0:nx,0:nz)
double precision, intent(out) :: zn(0:nx,0:nz)
type (deriv) :: dz
integer :: j, k

call deriv_x(zo, dz%x)
call deriv_xx(zo, dz%xx)
call deriv_zz(zo, dz%zz)

do k = 1, nz1
   zn(1:nx1,k) = zo(1:nx1,k) + (0.5d0 * rxx * dz%xx(1:nx1,k)) + &
                 (((1d0 - eta) * rx) / (4d0 * s(1:nx1))) * dz%x(1:nx1,k) - &
                 (((1d0 - eta)**2 * dt) / (2d0 * s(1:nx1)**2)) * &
                 zo(1:nx1,k) + 0.5d0 * rzz * dz%zz(1:nx1,k)
end do

return
END SUBROUTINE get_rhs_Zx

END MODULE linear
