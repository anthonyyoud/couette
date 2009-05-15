module magnetic 
  !Routines to do with solving the Poisson equation for the azimuthal magnetic
  !field.
  !Algorithms for the solution of the magnetic Poisson equation are as
  !for the current in current.f90
  implicit none

  private
  public :: b_poisson, fin_b_poisson

  contains

  subroutine b_poisson(u_mat, bn, b_mat, IPIV)
    !Solve Poisson equation for azimuthal magnetic field when tau=0
    use parameters
    use ic_bc, only : b_BCS
    implicit none

    double precision, intent(in) :: u_mat(0:nx,0:nz), &
                                    b_mat(2*nxp1+nxp1+1,0:nxp1*nz1-1)
    double precision, intent(out) :: bn(0:nx,0:nz)
    integer, intent(in) :: IPIV(nxp1*nz1)
    integer :: j, k, info
    double precision :: u_vec(0:nxp1*nz1-1)

    do k = 1, nz1
      do j = 0, nx
        u_vec(nxp1*(k-1)+j) = 0.5d0 * dx2 * delz * &
                              (u_mat(j,k-1) - u_mat(j,k+1))
      end do
    end do

    call DGBTRS('N', nxp1*nz1, nxp1, nxp1, 1, b_mat, 2*nxp1+nxp1+1, &
                 IPIV, u_vec, nxp1*nz1, info)
    if (info /= 0) then
      print*, 'ERROR: Solve error mag_inf_DGBTRS, INFO=', info
      stop
    end if

    do k = 1, nz1
      do j = 0, nx
        bn(j,k) = u_vec(nxp1*(k-1)+j)
      end do
    end do

    call b_BCS(bn)

    return
  end subroutine b_poisson

  subroutine fin_b_poisson(u_mat, bn, b_mat, IPIV)
    !Solve Poisson equation for azimuthal magnetic field when tau/=0
    use parameters
    use derivs, only : deriv_z
    implicit none

    double precision, intent(in) :: u_mat(0:nx,0:nz), &
                                    b_mat(2*nxp1+nxp1+1,0:nxp1*nzp1-1)
    double precision, intent(out) :: bn(0:nx,0:nz)
    integer, intent(in)  :: IPIV(nxp1*nzp1)
    integer :: j, k, info
    double precision :: u_vec(0:nxp1*nzp1-1), u_mat_z(0:nx,0:nz)

    call deriv_z(u_mat, u_mat_z)

    do k = 0, nz
      do j = 0, nx
        u_vec(nxp1*k+j) = -0.5d0 * dx2 * delz * u_mat_z(j,k)
      end do
    end do

    call DGBTRS('N', nxp1*nzp1, nxp1, nxp1, 1, b_mat, 2*nxp1+nxp1+1, &
                 IPIV, u_vec, nxp1*nzp1, info)
    if (info /= 0) then
      print*, 'ERROR: Solve error mag_fin_DGBTRS, INFO=', info
      stop
    end if

    do k = 0, nz
      do j = 0, nx
        bn(j,k) = u_vec(nxp1*k+j)
      end do
    end do

    return
  end subroutine fin_b_poisson

end module magnetic
