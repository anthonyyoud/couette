module stream
  !Routines to do with solving the Poisson equation for the stream function.
  !Algorithms for the solution of the stream function Poisson equation are as
  !for the current in current.f90
  implicit none

  private
  public :: p_poisson

  contains

  subroutine p_poisson(Z_mat, psi, p_mat, IPIV)
    !Solve Poisson equation for the stream function, psi for all tau
    use parameters
    use ic_bc, only : p_BCS, s
    implicit none

    real    (r2), intent(in)  :: Z_mat(0:nx,0:nz), &
                                 p_mat(2*nx1+nx1+1,nx1*nz1)
    real    (r2), intent(out) :: psi(0:nx,0:nz)
    integer (i1), intent(in)  :: IPIV(nx1*nz1)
    integer (i1)              :: j, k, info
    real    (r2)              :: zvec(nx1*nz1)

    do k = 1, nz1
      do j = 1, nx1
        zvec(nx1*(k-1)+j) = -s(j) * dx2 * dz2 * Z_mat(j,k)
      end do
    end do

    select case (r2)
      case (SPr)
        call SGBTRS('N', nx1*nz1, nx1, nx1, 1, p_mat, 2*nx1+nx1+1, &
                     IPIV, zvec, nx1*nz1, info)
        if (info /= 0) then
          print*, 'ERROR: SPr solve error psi_SGBTRS, INFO=', info
          stop
        end if
      case (DPr)
        call DGBTRS('N', nx1*nz1, nx1, nx1, 1, p_mat, 2*nx1+nx1+1, &
                     IPIV, zvec, nx1*nz1, info)
        if (info /= 0) then
          print*, 'ERROR: DPr solve error psi_DGBTRS, INFO=', info
          stop
        end if
      case default
        stop 'ERROR: Precision selection error - stream.f90 *GBTRS'
    end select

    do k = 1, nz1
      do j = 1, nx1
        psi(j,k) = zvec(nx1*(k-1)+j)
      end do
    end do

    call p_BCS(psi)

    return
  end subroutine p_poisson

end module stream
