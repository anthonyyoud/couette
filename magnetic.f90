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

    real    (r2), intent(in)  :: u_mat(0:nx,0:nz), &
                                 b_mat(2*nxp1+nxp1+1,0:nxp1*nz1-1)
    real    (r2), intent(out) :: bn(0:nx,0:nz)
    integer (i1), intent(in)  :: IPIV(nxp1*nz1)
    integer (i1)              :: j, k, info
    real    (r2)              :: u_vec(0:nxp1*nz1-1)

    do k = 1, nz1
      do j = 0, nx
        u_vec(nxp1*(k-1)+j) = 0.5_r2 * dx2 * delz * &
                              (u_mat(j,k-1) - u_mat(j,k+1))
      end do
    end do

    select case (r2)
      case (SPr)
        call SGBTRS('N', nxp1*nz1, nxp1, nxp1, 1, b_mat, 2*nxp1+nxp1+1, &
                     IPIV, u_vec, nxp1*nz1, info)
        if (info /= 0) then
          print*, 'ERROR: SPr solve error mag_inf_SGBTRS, INFO=', info
          stop
        end if
      case (DPr)
        call DGBTRS('N', nxp1*nz1, nxp1, nxp1, 1, b_mat, 2*nxp1+nxp1+1, &
                     IPIV, u_vec, nxp1*nz1, info)
        if (info /= 0) then
          print*, 'ERROR: DPr solve error mag_inf_DGBTRS, INFO=', info
          stop
        end if
      case default
        stop 'ERROR: Precision selection error - magnetic.f90 infinite *GBTRS'
    end select

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

    real    (r2), intent(in)  :: u_mat(0:nx,0:nz), &
                                 b_mat(2*nxp1+nxp1+1,0:nxp1*nzp1-1)
    real    (r2), intent(out) :: bn(0:nx,0:nz)
    integer (i1), intent(in)  :: IPIV(nxp1*nzp1)
    integer (i1)              :: j, k, info
    real    (r2)              :: u_vec(0:nxp1*nzp1-1), u_mat_z(0:nx,0:nz)

    call deriv_z(u_mat, u_mat_z)

    do k = 0, nz
      do j = 0, nx
        u_vec(nxp1*k+j) = -0.5_r2 * dx2 * delz * u_mat_z(j,k)
      end do
    end do

    select case (r2)
      case (SPr)
        call SGBTRS('N', nxp1*nzp1, nxp1, nxp1, 1, b_mat, 2*nxp1+nxp1+1, &
                     IPIV, u_vec, nxp1*nzp1, info)
        if (info /= 0) then
          print*, 'ERROR: SPr solve error mag_fin_SGBTRS, INFO=', info
          stop
        end if
      case (DPr)
        call DGBTRS('N', nxp1*nzp1, nxp1, nxp1, 1, b_mat, 2*nxp1+nxp1+1, &
                     IPIV, u_vec, nxp1*nzp1, info)
        if (info /= 0) then
          print*, 'ERROR: DPr solve error mag_fin_DGBTRS, INFO=', info
          stop
        end if
      case default
        stop 'ERROR: Precision selection error - magnetic.f90 finite *GBTRS'
    end select

    do k = 0, nz
      do j = 0, nx
        bn(j,k) = u_vec(nxp1*k+j)
      end do
    end do

    return
  end subroutine fin_b_poisson

end module magnetic
