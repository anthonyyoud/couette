module current 
  !Routines to do with solving the Poisson equation for the azimuthal current.
  implicit none

  private
  public :: j_poisson, fin_j_poisson

  contains

  subroutine j_poisson(p_mat, jn, j_mat, IPIV)
    !Solve Poisson equation for the azimuthal current when tau/=1
    use parameters
    use variables
    use ic_bc, only : j_BCS, s
    use derivs
    implicit none

    real, intent(in)  :: p_mat(0:nx,0:nz), j_mat(2*nx1+nx1+1,nx1*nzp1)
    real, intent(out) :: jn(0:nx,0:nz)
    integer, intent(in)  :: IPIV(nx1*nzp1)
    integer :: j, k, info
    real :: p_vec(nx1*nzp1)
    type (deriv)             :: dp

    call deriv_z(p_mat, dp%z)
    call deriv_x(dp%z, dp%zx)   !get derivatives for RHS
    call deriv_xx(dp%z, dp%zxx)
    call deriv_zz(dp%z, dp%zzz)

    !distribute RHS as a vector over the process grid, in block column format
    do k = 0, nz
      do j = 1, nx1
        p_vec(nx1*k+j) = dx2 * dz2 * &
                        (0.5 * dp%zzz(j,k) / (s(j) * delz**3) + &
                        0.5 * dp%zxx(j,k) / (s(j) * dx2 * delz) - &
                        0.25 * one_eta * dp%zx(j,k) / &
                        (s(j)**2 * delx * delz))!transform RHS matrix to vector
      end do
    end do

    !Solve Poisson equation using factorised matrix from D/SDBTRF
    call SGBTRS('N', nx1*nzp1, nx1, nx1, 1, j_mat, 2*nx1+nx1+1, &
                 IPIV, p_vec, nx1*nzp1, info)
    if (info /= 0) then
      print*, 'ERROR: SPr solve error cur_inf_SGBTRS, INFO=', info
      stop
    end if

    do k = 0, nz
      do j = 1, nx1
        jn(j,k) = p_vec(nx1*k+j)
      end do
    end do

    call j_BCS(jn)   !update boundary conditions

    return
  end subroutine j_poisson

  subroutine fin_j_poisson(p_mat, jn, j_mat, IPIV)
    !Solve Poisson equation for the azimuthal current when tau=1.
    !Algorithm as above but indices change to reflect different dimensions.
    use parameters
    use variables
    use ic_bc, only : j_BCS, s
    use derivs
    implicit none

    real, intent(in)  :: p_mat(0:nx,0:nz), j_mat(2*nx1+nx1+1,nx1*nz1)
    real, intent(out) :: jn(0:nx,0:nz)
    integer, intent(in)  :: IPIV(nx1*nz1)
    integer :: j, k, info
    real :: p_vec(nx1*nz1)
    type (deriv) :: dp

    call deriv_z(p_mat, dp%z)
    call deriv_x(dp%z, dp%zx)
    call deriv_xx(dp%z, dp%zxx)
    call deriv_zz(dp%z, dp%zzz)

    do k = 1, nz1
      do j = 1, nx1
        p_vec(nx1*(k-1)+j) = dx2 * dz2 * (0.5 * dp%zzz(j,k) / &
                             (s(j) * delz**3) + &
                             0.5 * dp%zxx(j,k) / (s(j) * dx2 * delz) - &
                             0.25 * one_eta * dp%zx(j,k) / &
                             (s(j)**2 * delx * delz))
      end do
    end do

    call SGBTRS('N', nx1*nz1, nx1, nx1, 1, j_mat, 2*nx1+nx1+1, &
                 IPIV, p_vec, nx1*nz1, info)
    if (info /= 0) then
      print*, 'ERROR: SPr solve error cur_fin_SGBTRS, INFO=', info
      stop
    end if

    do k = 1, nz1
      do j = 1, nx1
        jn(j,k) = p_vec(nx1*(k-1)+j)
      end do
    end do

    call j_BCS(jn)

    return
  end subroutine fin_j_poisson

end module current
