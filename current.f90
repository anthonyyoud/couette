module current 
  !Routines to do with solving the Poisson equation for the azimuthal current.
  implicit none

  private
  public :: j_poisson, fin_j_poisson

  contains

  subroutine j_poisson(p_mat, jn, j_mat, desc_j, af)
    !Solve Poisson equation for the azimuthal current when tau/=1
    use parameters
    use variables
    use ic_bc, only : j_BCS, s
    use derivs
    implicit none

    integer (i1),   intent(in)  :: desc_j(7)
    real    (r2),   intent(in)  :: af(laf), p_mat(0:nx,0:nz), j_mat(j_M,j_N)
    real    (r2),   intent(out) :: jn(0:nx,0:nz)
    real    (r2)                :: p_vec(nb), work(lwork_sol)
    integer (i1)                :: h, i, j, k, l, info, cpcol, desc_rp(7)
    type    (deriv)             :: dp

    desc_rp(1) = 502
    desc_rp(2) = ictxt
    desc_rp(3) = nx1*nzp1   !ScaLAPACK descriptor vector for RHS
    desc_rp(4) = nb
    desc_rp(5) = 0
    desc_rp(6) = nb

    if (mycol == 0) then
      call deriv_z(p_mat, dp%z)
      call deriv_x(dp%z, dp%zx)   !get derivatives for RHS
      call deriv_xx(dp%z, dp%zxx)
      call deriv_zz(dp%z, dp%zzz)
    end if

    call SLTIMER(9)
    if (npcol > 1) then
      !Broadcast RHS to all processes
      call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zx, nxp1, 0, 0) 
      call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zxx, nxp1, 0, 0)
      call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zzz, nxp1, 0, 0)
    end if
    call SLTIMER(9)

    cpcol = 0   !initialise current process column

    !distribute RHS as a vector over the process grid, in block column format
    do j = 1, nx1*nzp1, nb
      if (mycol == cpcol) then  !make sure each process gets right pieces
        do k = 1, min(nb, nx1*nzp1-j+1)
          i = k + j - 1
          h = modulo(i-1, nx1) + 1
          l = (i-1)/nx1
          p_vec(k) = dx2 * dz2 * (0.5_r2 * dp%zzz(h,l) / (s(h) * delz**3) + &
                     0.5_r2 * dp%zxx(h,l) / (s(h) * dx2 * delz) - &
                     0.25_r2 * one_eta * dp%zx(h,l) / &
                     (s(h)**2 * delx * delz))!transform RHS matrix to vector
        end do
      end if
      if (cpcol == npcol) exit  !if last process then exit
      cpcol = cpcol + 1  !next process column
    end do

    !Solve Poisson equation using factorised matrix from PD/PSDBTRF
    select case (r2)
      case (SPr)
        call PSDBTRS('N', nx1*nzp1, nx1, nx1, 1, j_mat, 1, desc_j, p_vec, 1, &
                      desc_rp, af, laf, work, lwork_sol, info)
        if (info /= 0) then
          print*, 'ERROR: SPr solve error cur_inf_PSDBTRS, INFO=', info
          stop
        end if
      case (DPr)
        call PDDBTRS('N', nx1*nzp1, nx1, nx1, 1, j_mat, 1, desc_j, p_vec, 1, &
                      desc_rp, af, laf, work, lwork_sol, info)
        if (info /= 0) then
          print*, 'ERROR: DPr solve error cur_inf_PDDBTRS, INFO=', info
          stop
        end if
      case default
        stop 'ERROR: Precision selection error - current.f90 infinite P*DBTRS'
    end select

    cpcol = 0   !reset current process column

    jn = 0.0_r2   !set matrix to zero on all processes
    do j = 1, nx1*nzp1, nb
      if (mycol == cpcol) then
        do k = 1, min(nb, nx1*nzp1-j+1)   !transform distributed RHS vector
          i = k + j - 1                     !into distributed matrix
          h = modulo(i-1, nx1) + 1
          l = (i-1)/nx1
          jn(h,l) = p_vec(k)
        end do
      end if
      if (cpcol == npcol) exit
      cpcol = cpcol + 1
    end do

    call SLTIMER(7)                                        !collect distributed
    if (npcol > 1) then
      call DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, jn, nxp1, 0, 0) !matrix onto
    end if
    call SLTIMER(7)                                            !master process

    if (mycol == 0) then
      call j_BCS(jn)   !update boundary conditions
    end if

    return
  end subroutine j_poisson

  subroutine fin_j_poisson(p_mat, jn, j_mat, desc_j, af)
    !Solve Poisson equation for the azimuthal current when tau=1.
    !Algorithm as above but indices change to reflect different dimensions.
    use parameters
    use variables
    use ic_bc, only : j_BCS, s
    use derivs
    implicit none

    integer (i1),   intent(in)  :: desc_j(7)
    real    (r2),   intent(in)  :: af(laf), p_mat(0:nx,0:nz), j_mat(j_M,j_N)
    real    (r2),   intent(out) :: jn(0:nx,0:nz)
    integer (i1)                :: h, i, j, k, l, info, cpcol, desc_rp(7)
    real    (r2)                :: work(lwork_sol), p_vec(nb)
    type    (deriv)             :: dp

    desc_rp(1) = 502
    desc_rp(2) = ictxt
    desc_rp(3) = nx1*nz1
    desc_rp(4) = nb
    desc_rp(5) = 0
    desc_rp(6) = nb

    if (mycol == 0) then
      call deriv_z(p_mat, dp%z)
      call deriv_x(dp%z, dp%zx)
      call deriv_xx(dp%z, dp%zxx)
      call deriv_zz(dp%z, dp%zzz)
    end if

    call SLTIMER(9)
    if (npcol > 1) then
      call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zx, nxp1, 0, 0)
      call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zxx, nxp1, 0, 0)
      call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zzz, nxp1, 0, 0)
    end if
    call SLTIMER(9)

    cpcol = 0

    do j = 1, nx1*nz1, nb
      if (mycol == cpcol) then
        do k = 1, min(nb, nx1*nz1-j+1)
           i = k + j - 1
           h = modulo(i-1, nx1) + 1
           l = (i-1)/nx1 + 1
           p_vec(k) = dx2 * dz2 * (0.5_r2 * dp%zzz(h,l) / &
                      (s(h) * delz**3) + &
                      0.5_r2 * dp%zxx(h,l) / (s(h) * dx2 * delz) - &
                      0.25_r2 * one_eta * dp%zx(h,l) / &
                      (s(h)**2 * delx * delz))
        end do
      end if
      if (cpcol == npcol) exit
      cpcol = cpcol + 1
    end do

    select case (r2)
      case (SPr)
        call PSDBTRS('N', nx1*nz1, nx1, nx1, 1, j_mat, 1, desc_j, p_vec, 1, &
                      desc_rp, af, laf, work, lwork_sol, info)
        if (info /= 0) then
          print*, 'ERROR: SPr solve error cur_fin_PSDBTRS, INFO=', info
          stop
        end if
      case (DPr)
        call PDDBTRS('N', nx1*nz1, nx1, nx1, 1, j_mat, 1, desc_j, p_vec, 1, &
                      desc_rp, af, laf, work, lwork_sol, info)
        if (info /= 0) then
          print*, 'ERROR: DPr solve error cur_fin_PDDBTRS, INFO=', info
          stop
        end if
      case default
        stop 'ERROR: Precision selection error - current.f90 finite P*DBTRS'
    end select

    cpcol = 0

    jn = 0.0_r2
    do j = 1, nx1*nz1, nb
      if (mycol == cpcol) then
        do k = 1, min(nb, nx1*nz1-j+1)
           i = k + j - 1
           h = modulo(i-1, nx1) + 1
           l = (i-1)/nx1 + 1
           jn(h,l) = p_vec(k)
        end do
      end if
      if (cpcol == npcol) exit
      cpcol = cpcol + 1
    end do

    call SLTIMER(7)
    if (npcol > 1) then
      call DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, jn, nxp1, 0, 0)
    end if
    call SLTIMER(7)

    if (mycol == 0) then
      call j_BCS(jn)
    end if

    return
  end subroutine fin_j_poisson

end module current
