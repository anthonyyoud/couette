! $Id: linear.f90 1 2010/01/30 13:54:36 youd $
!----------------------------------------------------------------------------

module linear
  !Routines to get the linear parts of the right-hand side of the azimuthal
  !velocity/vorticity equations.
  implicit none

  private
  public :: get_rhs_ux, get_rhs_Zx

  contains

  subroutine get_rhs_ux(uo, u)
    !Get the linear part for the right-hand side of the tridiagonal system for u
    use parameters
    use variables
    use derivs
    use ic_bc, only : s
    implicit none

    double precision, intent(in) :: uo(0:nx,0:nz)
    double precision, intent(out) :: u(0:nx,0:nz)
    integer :: j, k
    type (deriv) :: du

    call deriv_x(uo, du%x)
    call deriv_xx(uo, du%xx)   !get derivatives
    call deriv_zz(uo, du%zz)

    do k = 1, nz1
      u(1:nx1,k) = uo(1:nx1,k) + 0.5d0 * rxx * du%xx(1:nx1,k) + &
                   one_eta * rx * 0.25d0 * du%x(1:nx1,k) / s(1:nx1) - &
                   one_eta**2 * dt * 0.5d0 * uo(1:nx1,k) / s(1:nx1)**2 + &
                   0.5d0 * rzz * du%zz(1:nx1,k)
    end do

    if (abs(tau - 1d0) > epsilon(tau)) then
      u(1:nx1,0) = uo(1:nx1,0) + 0.5d0 * rxx * du%xx(1:nx1,0) + &
                   one_eta * rx * 0.25d0 * du%x(1:nx1,0) / s(1:nx1) - &
                   one_eta**2 * dt * 0.5d0 * uo(1:nx1,0) / s(1:nx1)**2 + &
                   0.5d0 * rzz * du%zz(1:nx1,0)

      u(1:nx1,nz) = uo(1:nx1,nz) + 0.5d0 * rxx * du%xx(1:nx1,nz) + &
                    one_eta * rx * 0.25d0 * du%x(1:nx1,nz) / s(1:nx1) - &
                    one_eta**2 * dt * 0.5d0 * uo(1:nx1,nz) / s(1:nx1)**2 + &
                    0.5d0 * rzz * du%zz(1:nx1,nz)
    end if

    return
  end subroutine get_rhs_ux

  subroutine get_rhs_Zx(zo, zn)
    !Get linear part for the right-hand side of the tridiagonal system for Z
    use parameters
    use variables
    use derivs
    use ic_bc, only : s
    implicit none

    double precision, intent(in) :: zo(0:nx,0:nz)
    double precision, intent(out) :: zn(0:nx,0:nz)
    integer :: j, k
    type (deriv) :: dz

    call deriv_x(zo, dz%x)
    call deriv_xx(zo, dz%xx)   !get derivatives
    call deriv_zz(zo, dz%zz)

    do k = 1, nz1
      zn(1:nx1,k) = zo(1:nx1,k) + 0.5d0 * rxx * dz%xx(1:nx1,k) + &
                    one_eta * rx * 0.25d0 * dz%x(1:nx1,k) / s(1:nx1) - &
                    one_eta**2 * dt * 0.5d0 * zo(1:nx1,k) / s(1:nx1)**2 + &
                    0.5d0 * rzz * dz%zz(1:nx1,k)
    end do

    return
  end subroutine get_rhs_Zx

end module linear
