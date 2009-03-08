module nonlinear
  !Routines to get the nonlinear parts of the right-hand side of the azimuthal
  !velocity/vorticity equations.
  implicit none

  private
  public :: get_nlin_ux, get_nlin_Zx

  contains

  subroutine get_nlin_ux(uo, uo2, po, po2, bo, bo2, u_nl_n)
    !Get nonlinear part for the right-hand side of the tridiagonal system for u
    use parameters
    use variables
    use derivs
    use ic_bc, only : s 
    implicit none

    real, intent(in)  :: uo(0:nx,0:nz), uo2(0:nx,0:nz), &
                         po(0:nx,0:nz), po2(0:nx,0:nz), &
                         bo(0:nx,0:nz), bo2(0:nx,0:nz)
    real, intent(out) :: u_nl_n(0:nx,0:nz) 
    integer :: j, k  
    type (deriv) :: du, du2, dp, dp2, db, db2

    call deriv_x(uo, du%x)
    call deriv_x(uo2, du2%x)
    call deriv_z(uo, du%z)
    call deriv_z(uo2, du2%z)
    call deriv_x(po, dp%x)
    call deriv_x(po2, dp2%x)   !get derivatives
    call deriv_z(po, dp%z)
    call deriv_z(po2, dp2%z)
    call deriv_z(bo, db%z)
    call deriv_z(bo2, db2%z)

    do k = 1, nz1
      u_nl_n(1:nx1,k) = (-0.125 * rx / (s(1:nx1) * delz)) * &
                        (3.0 * (dp%x(1:nx1,k) * du%z(1:nx1,k) - &
                        dp%z(1:nx1,k) * du%x(1:nx1,k)) - &
                        (dp2%x(1:nx1,k) * du2%z(1:nx1,k) - &
                        dp2%z(1:nx1,k) * du2%x(1:nx1,k))) + &
                        (one_eta * rz * 0.25 / s(1:nx1)**2) * &
                        (3.0 * uo(1:nx1,k) * dp%z(1:nx1,k) - &
                        uo2(1:nx1,k) * dp2%z(1:nx1,k)) + &
                        0.25 * rz * Q * (3.0 * db%z(1:nx1,k) - &
                        db2%z(1:nx1,k))
    end do

    if (abs(tau - 1.0) > epsilon(tau)) then
      u_nl_n(1:nx1,0) = (-0.125 * rx / (s(1:nx1) * delz)) * &
                        (-3.0 * dp%z(1:nx1,0) * du%x(1:nx1,0) + &
                        dp2%z(1:nx1,0) * du2%x(1:nx1,0)) + &
                        (one_eta * rz * 0.25 / s(1:nx1)**2) * &
                        (3.0 * uo(1:nx1,0) * dp%z(1:nx1,0) - &
                        uo2(1:nx1,0) * dp2%z(1:nx1,0)) + &
                        0.25 * rz * Q * (3.0 * db%z(1:nx1,0) - &
                        db2%z(1:nx1,0))

      u_nl_n(1:nx1,nz) = (-0.125 * rx / (s(1:nx1) * delz)) * &
                         (-3.0 * dp%z(1:nx1,nz) * du%x(1:nx1,nz) + &
                         dp2%z(1:nx1,nz) * du2%x(1:nx1,nz)) + &
                         (one_eta * rz * 0.25 / s(1:nx1)**2) * &
                        (3.0 * uo(1:nx1,nz) * dp%z(1:nx1,nz) - &
                         uo2(1:nx1,nz) * dp2%z(1:nx1,nz)) + &
                         0.25 * rz * Q * (3.0 * db%z(1:nx1,nz) - &
                         db2%z(1:nx1,nz))
    end if

    return
  end subroutine get_nlin_ux

  subroutine get_nlin_Zx(t, uo, uo2, po, po2, zo, zo2, jo, jo2, z_nl_n)
    !Get nonlinear part for the right-hand side of the tridiagonal system for Z
    use parameters
    use variables
    use derivs
    use ic_bc, only : s
    implicit none

    real, intent(in)  :: t, uo(0:nx,0:nz), uo2(0:nx,0:nz), &
                         po(0:nx,0:nz), po2(0:nx,0:nz), &
                         zo(0:nx,0:nz), zo2(0:nx,0:nz), &
                         jo(0:nx,0:nz), jo2(0:nx,0:nz)
    real, intent(out) :: z_nl_n(0:nx,0:nz)
    integer :: j, k
    type (deriv) :: du, du2, dp, dp2, dz, dz_2, dj, dj2

    call deriv_z(uo, du%z)
    call deriv_z(uo2, du2%z)
    call deriv_x(po, dp%x)
    call deriv_x(po2, dp2%x)
    call deriv_z(po, dp%z)
    call deriv_z(po2, dp2%z)
    call deriv_x(zo, dz%x)    !get derivatives
    call deriv_x(zo2, dz_2%x)
    call deriv_z(zo, dz%z)
    call deriv_z(zo2, dz_2%z)
    call deriv_z(jo, dj%z)
    call deriv_z(jo2, dj2%z)

    do k = 1, nz1
      z_nl_n(1:nx1,k) = (one_eta * rz * 0.5 / s(1:nx1)) * &
                        (3.0 * uo(1:nx1,k) * du%z(1:nx1,k) - &
                        uo2(1:nx1,k) * du2%z(1:nx1,k)) - &
                        (one_eta * rz * 0.25 / s(1:nx1)**2) * &
                        (3.0 * zo(1:nx1,k) * dp%z(1:nx1,k) - &
                        zo2(1:nx1,k) * dp2%z(1:nx1,k)) - &
                        (0.125 * rx / (s(1:nx1) * delz)) * &
                        ((3.0 * (dp%x(1:nx1,k) * dz%z(1:nx1,k) - &
                        dp%z(1:nx1,k) * dz%x(1:nx1,k))) - &
                        (dp2%x(1:nx1,k) * dz_2%z(1:nx1,k) - &
                        dp2%z(1:nx1,k) * dz_2%x(1:nx1,k))) + &
                        0.25 * rz * Q * (3.0 * dj%z(1:nx1,k) - &
                        dj2%z(1:nx1,k))
    end do

    return
  end subroutine get_nlin_Zx

end module nonlinear
