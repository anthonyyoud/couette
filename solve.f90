! $Id$
!----------------------------------------------------------------------------

module solve
  !Routines to solve the equations for the azimuthal velocity/vorticity.
  implicit none

  private
  public :: solve_ux, solve_uz, solve_Zx, solve_Zz

  contains

  subroutine solve_ux(uo, u, u_nl, t, ux)
    !Solve for the azimuthal velocity field in the x-direction
    use parameters
    use variables
    use ic_bc, only : u_BCS, s
    implicit none

    double precision, intent(in) :: t, u(0:nx,0:nz), u_nl(0:nx,0:nz)
    double precision, intent(inout) :: uo(0:nx,0:nz)
    type (mat_comp), intent(in) :: ux
    integer :: j, k
    double precision :: ux_rhs(nx1)

    call u_BCS(uo, t)

    do k = 1, nz1
      ux_rhs(:) = u(1:nx1,k) + u_nl(1:nx1,k)   !RHS at fixed z, looping over x

      ux_rhs(1) = ux_rhs(1) + 0.5d0 * rxx * uo(0,k) - &
                  one_eta * rx * 0.25d0 * uo(0,k) / s(1)    !BCS
      ux_rhs(nx1) = ux_rhs(nx1) + 0.5d0 * rxx * uo(nx,k) + &
                    one_eta * rx * 0.25d0 * uo(nx,k) / s(nx1)

      call thomas(xlb, nx1, ux%up, ux%di, ux%lo, ux_rhs)   !Thomas algorithm
                                                           !at each z
      uo(1:nx1,k) = ux_rhs(:)   !intermediate solution
    end do

    if (abs(tau - 1d0) > epsilon(tau)) then !extra entries due to BCS
      do k = 0, nz, nz
        ux_rhs(:) = u(1:nx1,k) + u_nl(1:nx1,k)

        ux_rhs(1) = ux_rhs(1) + 0.5d0 * rxx * uo(0,k) - &
                    one_eta * rx * 0.25d0 * uo(0,k) / s(1)
        ux_rhs(nx1) = ux_rhs(nx1) + 0.5d0 * rxx * uo(nx,k) + &
                      one_eta * rx * 0.25d0 * uo(nx,k) / s(nx1)

        call thomas(xlb, nx1, ux%up, ux%di, ux%lo, ux_rhs)

        uo(1:nx1,k) = ux_rhs(:)
      end do
    end if

    return
  end subroutine solve_ux

  subroutine solve_Zx(zo, zn, z_nl, po, t, zx)
    !Solve for the azimuthal vorticity in the x-direction
    use parameters
    use variables
    use ic_bc, only : z_BCS, s
    implicit none

    double precision, intent(in) :: t, zn(0:nx,0:nz), po(0:nx,0:nz), &
                        z_nl(0:nx,0:nz)
    double precision, intent(inout) :: zo(0:nx,0:nz)
    type (mat_comp), intent(in) :: zx
    integer :: j, k
    double precision :: zx_rhs(nx1)

    call z_BCS(zo, po, t)

    do k = 1, nz1
      zx_rhs(:) = zn(1:nx1,k) + z_nl(1:nx1,k)   !RHS at fixed z, looping over x

      zx_rhs(1) = zx_rhs(1) + 0.5d0 * rxx * zo(0,k) - &
                  one_eta * rx * 0.25d0 * zo(0,k) / s(1)   !BCS
      zx_rhs(nx1) = zx_rhs(nx1) + 0.5d0 * rxx * zo(nx,k) + &
                    one_eta * rx * 0.25d0 * zo(nx,k) / s(nx1)

      call thomas(xlb, nx1, zx%up, zx%di, zx%lo, zx_rhs)   !Thomas algorithm

      zo(1:nx1,k) = zx_rhs(:)   !intermediate solution
    end do

    return
  end subroutine solve_Zx

  subroutine solve_uz(uo, u, t, uz)
    !Solve for azimuthal velocity in z-direction to give full solution
    use parameters
    use variables
    use ic_bc, only : u_BCS, s
    implicit none

    double precision, intent(in) :: t, uo(0:nx,0:nz)
    double precision, intent(inout) :: u(0:nx,0:nz)
    type (uz_mat_comp), intent(in) :: uz
    integer :: j, k
    double precision :: uz_rhs(0:nz), uz_rhs_t1(nz1), &
                        up(nz-2), di(nz1), lo(2:nz1)

    call u_BCS(u, t)

    if (abs(tau - 1d0) < epsilon(tau)) then
      up(:) = uz%up(1:nz-2)
      di(:) = uz%di(1:nz1)   !dimensions of upper, lower, diagonal change
      lo(:) = uz%lo(2:nz1)   !for tau=1

      do j = 1, nx1
        uz_rhs_t1(:) = uo(j,1:nz1)   !RHS at fixed x looping over z

        uz_rhs_t1(1) = uz_rhs_t1(1) + 0.5d0 * rzz * u(j,0)
        uz_rhs_t1(nz1) = uz_rhs_t1(nz1) + 0.5d0 * rzz * u(j,nz)   !BCS

        call thomas(zlb+1, nz1, up, di, lo, uz_rhs_t1)   !Thomas algorithm

        u(j,1:nz1) = uz_rhs_t1(:)   !full solution
      end do
    else   !straight-forward if tau /= 1
      do j = 1, nx1
        uz_rhs(:) = uo(j,:)   !RHS at fixed x looping over z

        call thomas(zlb, nz, uz%up, uz%di, uz%lo, uz_rhs)   !Thomas algorithm

        u(j,:) = uz_rhs(:)   !full solution
      end do
    end if

    return
  end subroutine solve_uz

  subroutine solve_Zz(zo, po, zn, t, zz)
    !Solve for azimuthal vorticity in z-direction to give full solution
    use parameters
    use variables
    use ic_bc, only : z_BCS, s
    implicit none

    double precision, intent(in) :: t, zo(0:nx,0:nz), po(0:nx,0:nz)
    double precision, intent(inout) :: zn(0:nx,0:nz)
    type (zz_mat_comp), intent(in) :: zz
    integer :: j, k
    double precision :: Zz_rhs(nz1)

    call z_BCS(zn, po, t)

    do j = 1, nx1
      Zz_rhs(:) = zo(j,1:nz1)   !RHS at fixed x, looping over z

      Zz_rhs(1) = Zz_rhs(1) + 0.5d0 * rzz * zn(j,0)
      Zz_rhs(nz1) = Zz_rhs(nz1) + 0.5d0 * rzz * zn(j,nz)   !BCS

      call thomas(zlb+1, nz1, zz%up, zz%di, zz%lo, Zz_rhs)   !Thomas algorithm

      zn(j,1:nz1) = Zz_rhs(:)   !full solution
    end do

    return
  end subroutine solve_Zz

  subroutine thomas (lb, m, up, di, lo, r)
    !Thomas algorithm for solving a tridiagonal system of equations
    implicit none

    integer, intent(in) :: m, lb  !lb is the vector lower-bound
    double precision, intent(in) :: up(lb:m-1), di(lb:m), lo(lb+1:m)
    double precision, intent(inout) :: r(lb:m)
    double precision :: dnew(lb:m), aa = 0d0
    integer :: j

    dnew = di   !create new diagonal so as not to destroy original
    do j = lb+1, m
      aa = -lo(j) / dnew(j-1)
      dnew(j) = dnew(j) + aa * up(j-1)   !eliminate lower-diagonal
      r(j) = r(j) + aa * r(j-1)
    end do

    r(m) = r(m) / dnew(m)   !first step back-substitution

    do j = m-1, lb, -1
      r(j) = (r(j) - up(j) * r(j+1)) / dnew(j)  !solve by back-substitution
    end do

    return
  end subroutine thomas

end module solve
