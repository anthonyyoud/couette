module stream
  !Routines to do with solving the Poisson equation for the stream function.
  !Algorithms for the solution of the stream function Poisson equation are as
  !for the current in current.f90
  implicit none

  private
  public :: p_poisson

  contains

  subroutine p_poisson(Z_mat, psi, p_mat, desc_p, af)
    !Solve Poisson equation for the stream function, psi for all tau
    use parameters
    use ic_bc, only : p_BCS, s
    implicit none

    integer (i1), intent(in)  :: desc_p(7)
    real    (r2), intent(in)  :: af(laf), Z_mat(0:nx,0:nz), p_mat(p_M,p_N)
    real    (r2), intent(out) :: psi(0:nx,0:nz)
    integer (i1)              :: i, j, k, info, cpcol, desc_z(7)
    real    (r2)              :: zvec(nb), work(lwork_sol)

    desc_z(1) = 502
    desc_z(2) = ictxt
    desc_z(3) = nx1*nz1
    desc_z(4) = nb
    desc_z(5) = 0
    desc_z(6) = nb

    cpcol = 0

    do j = 1, nx1*nz1, nb
      if (mycol == cpcol) then
        do k = 1, min(nb, nx1*nz1-j+1)
          i = k + j - 1
          zvec(k) = -s(modulo(i-1, nx1) + 1) * dx2 * dz2 * &
                     Z_mat(modulo(i-1, nx1) + 1, (i-1)/nx1 + 1)
        end do
      end if
      if (cpcol == npcol) exit
      cpcol = cpcol + 1
    end do

    select case (r2)
      case (SPr)
        call PSDBTRS('N', nx1*nz1, nx1, nx1, 1, p_mat, 1, desc_p, zvec, 1, &
                      desc_z, af, laf, work, lwork_sol, info)
        if (info /= 0) then
          print*, 'ERROR: SPr solve error psi_PSDBTRS, INFO=', info
          stop
        end if
      case (DPr)
        call PDDBTRS('N', nx1*nz1, nx1, nx1, 1, p_mat, 1, desc_p, zvec, 1, &
                      desc_z, af, laf, work, lwork_sol, info)
        if (info /= 0) then
          print*, 'ERROR: DPr solve error psi_PDDBTRS, INFO=', info
          stop
        end if
      case default
        stop 'ERROR: Precision selection error - stream.f90 P*DBTRS'
    end select

    cpcol = 0

    psi = 0.0_r2
    do j = 1, nx1*nz1, nb
      if (mycol == cpcol) then
        do k = 1, min(nb, nx1*nz1-j+1)
          i = k + j - 1
          psi(modulo(i-1, nx1) + 1, (i-1)/nx1 + 1) = zvec(k)
        end do
      end if
      if (cpcol == npcol) exit
      cpcol = cpcol + 1
    end do

    call SLTIMER(6)
    if (npcol > 1) then
      call DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, psi, nxp1, 0, 0)
    end if
    call SLTIMER(6)

    if (mycol == 0) then
      call p_BCS(psi)
    end if

    return
  end subroutine p_poisson

end module stream
