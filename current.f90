MODULE current 
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: j_poisson, fin_j_poisson

  CONTAINS

  SUBROUTINE j_poisson(p_mat, jn, j_mat, desc_j, af)
    !Solve Poisson equation for the azimuthal current when tau/=1
    USE parameters
    USE variables
    USE ic_bc, ONLY : j_BCS, s
    USE derivs
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN)  :: desc_j(7)
    REAL (r2),    INTENT(IN)  :: af(laf), p_mat(0:nx,0:nz), j_mat(j_M,j_N)
    REAL (r2),    INTENT(OUT) :: jn(0:nx,0:nz)
    REAL (r2)                 :: p_vec(nb), work(lwork_sol)
    INTEGER (i1)              :: h, i, j, k, l, info, cpcol, desc_rp(7)
    TYPE (DERIV)              :: dp

    desc_rp(1) = 502
    desc_rp(2) = ictxt
    desc_rp(3) = nx1*nzp1   !ScaLAPACK descriptor vector for RHS
    desc_rp(4) = nb
    desc_rp(5) = 0
    desc_rp(6) = nb

    IF (mycol == 0) THEN
      CALL deriv_z(p_mat, dp%z)
      CALL deriv_x(dp%z, dp%zx)   !get derivatives for RHS
      CALL deriv_xx(dp%z, dp%zxx)
      CALL deriv_zz(dp%z, dp%zzz)
    END IF

    CALL SLTIMER(9)
    IF (npcol > 1) THEN
      !Broadcast RHS to all processes
      CALL DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zx, nxp1, 0, 0) 
      CALL DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zxx, nxp1, 0, 0)
      CALL DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zzz, nxp1, 0, 0)
    END IF
    CALL SLTIMER(9)

    cpcol = 0   !initialise current process column

    !distribute RHS as a vector over the process grid, in block column format
    DO j = 1, nx1*nzp1, nb
      IF (mycol == cpcol) THEN  !make sure each process gets right pieces
        DO k = 1, MIN(nb, nx1*nzp1-j+1)
          i = k + j - 1
          h = MODULO(i-1, nx1) + 1
          l = (i-1)/nx1
          p_vec(k) = dx2 * dz2 * (0.5_r2 * dp%zzz(h,l) / (s(h) * delz**3) + &
                     0.5_r2 * dp%zxx(h,l) / (s(h) * dx2 * delz) - &
                     0.25_r2 * one_eta * dp%zx(h,l) / &
                     (s(h)**2 * delx * delz))!transform RHS matrix to vector
        END DO
      END IF
      IF (cpcol == npcol) EXIT  !if last process then exit
      cpcol = cpcol + 1  !next process column
    END DO

    !Solve Poisson equation using factorised matrix from PDDBTRF
    CALL PDDBTRS('N', nx1*nzp1, nx1, nx1, 1, j_mat, 1, desc_j, p_vec, 1, &
                  desc_rp, af, laf, work, lwork_sol, info)
    IF (info /= 0) PRINT*, 'j_infinite_PDDBTRS ', info

    cpcol = 0   !reset current process column

    jn = 0_r2   !set matrix to zero on all processes
    DO j = 1, nx1*nzp1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nx1*nzp1-j+1)   !transform distributed RHS vector
          i = k + j - 1                     !into distributed matrix
          h = MODULO(i-1, nx1) + 1
          l = (i-1)/nx1
          jn(h,l) = p_vec(k)
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    CALL SLTIMER(7)                                        !collect distributed
    IF (npcol > 1) THEN
      CALL DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, jn, nxp1, 0, 0) !matrix onto
    END IF
    CALL SLTIMER(7)                                            !master process

    IF (mycol == 0) THEN
      CALL j_BCS(jn)   !update boundary conditions
    END IF

    RETURN
  END SUBROUTINE j_poisson

  SUBROUTINE fin_j_poisson(p_mat, jn, j_mat, desc_j, af)
    !Solve Poisson equation for the azimuthal current when tau=1.
    !Algorithm as above but indices change to reflect different dimensions.
    USE parameters
    USE variables
    USE ic_bc, ONLY : j_BCS, s
    USE derivs
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN)  :: desc_j(7)
    REAL (r2),    INTENT(IN)  :: af(laf), p_mat(0:nx,0:nz), j_mat(j_M,j_N)
    REAL (r2),    INTENT(OUT) :: jn(0:nx,0:nz)
    REAL (r2)                 :: p_vec(nb)
    INTEGER (i1)              :: h, i, j, k, l, info, cpcol, desc_rp(7)
    REAL (r2)                 :: work(lwork_sol)
    TYPE (DERIV)              :: dp

    desc_rp(1) = 502
    desc_rp(2) = ictxt
    desc_rp(3) = nx1*nz1
    desc_rp(4) = nb
    desc_rp(5) = 0
    desc_rp(6) = nb

    IF (mycol == 0) THEN
      CALL deriv_z(p_mat, dp%z)
      CALL deriv_x(dp%z, dp%zx)
      CALL deriv_xx(dp%z, dp%zxx)
      CALL deriv_zz(dp%z, dp%zzz)
    END IF

    CALL SLTIMER(9)
    IF (npcol > 1) THEN
      CALL DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zx, nxp1, 0, 0)
      CALL DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zxx, nxp1, 0, 0)
      CALL DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, dp%zzz, nxp1, 0, 0)
    END IF
    CALL SLTIMER(9)

    cpcol = 0

    DO j = 1, nx1*nz1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nx1*nz1-j+1)
           i = k + j - 1
           h = MODULO(i-1, nx1) + 1
           l = (i-1)/nx1 + 1
           p_vec(k) = dx2 * dz2 * (0.5_r2 * dp%zzz(h,l) / &
                      (s(h) * delz**3) + &
                      0.5_r2 * dp%zxx(h,l) / (s(h) * dx2 * delz) - &
                      0.25_r2 * one_eta * dp%zx(h,l) / &
                      (s(h)**2 * delx * delz))
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    CALL PDDBTRS('N', nx1*nz1, nx1, nx1, 1, j_mat, 1, desc_j, p_vec, 1, &
                  desc_rp, af, laf, work, lwork_sol, info)
    IF (info /= 0) PRINT*, 'j_finite_PDDBTRS ', info

    cpcol = 0

    jn = 0.0_r2
    DO j = 1, nx1*nz1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nx1*nz1-j+1)
           i = k + j - 1
           h = MODULO(i-1, nx1) + 1
           l = (i-1)/nx1 + 1
           jn(h,l) = p_vec(k)
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    CALL SLTIMER(7)
    IF (npcol > 1) THEN
      CALL DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, jn, nxp1, 0, 0)
    END IF
    CALL SLTIMER(7)

    IF (mycol == 0) THEN
      CALL j_BCS(jn)
    END IF

    RETURN
  END SUBROUTINE fin_j_poisson

END MODULE current
