MODULE stream
  !Algorithms for the solution of the stream-function Poisson equation are as
  !for the current in current.f90
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: p_poisson

  CONTAINS

  SUBROUTINE p_poisson(Z_mat, psi, p_mat, desc_p, af)
    !Solve Poisson equation for the stream-function, psi for all tau
    USE parameters
    USE ic_bc, ONLY : p_BCS, s
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN)  :: desc_p(7)
    REAL (r2),    INTENT(IN)  :: af(laf), Z_mat(0:nx,0:nz), p_mat(p_M,p_N)
    REAL (r2),    INTENT(OUT) :: psi(0:nx,0:nz)
    REAL (r2)                 :: zvec(nb), work(lwork_sol)
    INTEGER (i1)              :: i, j, k, info, cpcol, desc_z(7)

    desc_z(1) = 502
    desc_z(2) = ictxt
    desc_z(3) = nx1*nz1
    desc_z(4) = nb
    desc_z(5) = 0
    desc_z(6) = nb

    cpcol = 0

    DO j = 1, nx1*nz1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nx1*nz1-j+1)
          i = k + j - 1
          zvec(k) = -s(MODULO(i-1, nx1) + 1) * dx2 * dz2 * &
                     Z_mat(MODULO(i-1, nx1) + 1, (i-1)/nx1 + 1)
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    CALL PDDBTRS('N', nx1*nz1, nx1, nx1, 1, p_mat, 1, desc_p, zvec, 1, &
                  desc_z, af, laf, work, lwork_sol, info)
    IF (info /= 0) PRINT*, 'psi_PDDBTRS ', info

    cpcol = 0

    psi = 0.0_r2
    DO j = 1, nx1*nz1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nx1*nz1-j+1)
          i = k + j - 1
          psi(MODULO(i-1, nx1) + 1, (i-1)/nx1 + 1) = zvec(k)
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    CALL SLTIMER(6)
    IF (npcol > 1) THEN
      CALL DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, psi, nxp1, 0, 0)
    END IF
    CALL SLTIMER(6)

    IF (mycol == 0) THEN
      CALL p_BCS(psi)
    END IF

    RETURN
  END SUBROUTINE p_poisson

END MODULE stream
