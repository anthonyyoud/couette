MODULE magnetic 
  !Routines to do with solving the Poisson equation for the azimuthal magnetic
  !field.
  !Algorithms for the solution of the magnetic Poisson equation are as
  !for the current in current.f90
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: b_poisson, fin_b_poisson

  CONTAINS

  SUBROUTINE b_poisson(u_mat, bn, b_mat, desc_b, af)
    !Solve Poisson equation for azimuthal magnetic field when tau=0
    USE parameters
    USE ic_bc, ONLY : b_BCS
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN)  :: desc_b(7)
    REAL    (r2), INTENT(IN)  :: af(b_laf), u_mat(0:nx,0:nz), b_mat(b_M,b_N)
    REAL    (r2), INTENT(OUT) :: bn(0:nx,0:nz)
    INTEGER (i1)              :: i, j, k, info, cpcol, desc_u(7)
    REAL    (r2)              :: u_vec(nb), work(lwork_b_sol)

    desc_u(1) = 502
    desc_u(2) = ictxt
    desc_u(3) = nxp1*nz1
    desc_u(4) = nb
    desc_u(5) = 0
    desc_u(6) = nb

    cpcol = 0

    DO j = 0, nxp1*nz1-1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nxp1*nz1-j)
          i = k + j - 1
          u_vec(k) = 0.5_r2 * dx2 * delz * &
                     (u_mat(MODULO(i, nxp1), i/nxp1) - &
                      u_mat(MODULO(i, nxp1), i/nxp1 + 2))
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    SELECT CASE (r2)
      CASE (SPr)
        CALL PSDBTRS('N', nxp1*nz1, nxp1, nxp1, 1, b_mat, 1, desc_b, u_vec, 1, &
                      desc_u, af, b_laf, work, lwork_b_sol, info)
        IF (info /= 0) PRINT*, 'b_infinite_PSDBTRS ', info
      CASE (DPr)
        CALL PDDBTRS('N', nxp1*nz1, nxp1, nxp1, 1, b_mat, 1, desc_b, u_vec, 1, &
                      desc_u, af, b_laf, work, lwork_b_sol, info)
        IF (info /= 0) PRINT*, 'b_infinite_PDDBTRS ', info
      CASE DEFAULT
        STOP 'ERROR: Precision selection error - magnetic.f90 infinite P*DBTRS'
    END SELECT

    cpcol = 0

    bn = 0.0_r2
    DO j = 0, nxp1*nz1-1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nxp1*nz1-j)
          i = k + j - 1
          bn(MODULO(i, nxp1), i/nxp1 + 1) = u_vec(k)
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    CALL SLTIMER(8)
    IF (npcol > 1) THEN
      CALL DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, bn, nxp1, 0, 0)
    END IF
    CALL SLTIMER(8)

    IF (mycol == 0) THEN
      CALL b_BCS(bn)
    END IF

    RETURN
  END SUBROUTINE b_poisson

  SUBROUTINE fin_b_poisson(u_mat, bn, b_mat, desc_b, af)
    !Solve Poisson equation for azimuthal magnetic field when tau/=0
    USE parameters
    USE derivs, ONLY : deriv_z
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN)  :: desc_b(7)
    REAL    (r2), INTENT(IN)  :: af(b_laf), u_mat(0:nx,0:nz), b_mat(b_M,b_N)
    REAL    (r2), INTENT(OUT) :: bn(0:nx,0:nz)
    INTEGER (i1)              :: i, j, k, info, cpcol, desc_u(7)
    REAL    (r2)              :: u_vec(nb), u_mat_z(0:nx,0:nz), &
                                 work(lwork_b_sol)

    IF (mycol == 0) THEN
      CALL deriv_z(u_mat, u_mat_z)
    END IF

    IF (npcol > 1) THEN
      CALL DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, u_mat_z, nxp1, 0, 0)
    END IF

    desc_u(1) = 502
    desc_u(2) = ictxt
    desc_u(3) = nxp1*nzp1
    desc_u(4) = nb
    desc_u(5) = 0
    desc_u(6) = nb

    cpcol = 0

    DO j = 0, nxp1*nzp1-1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nxp1*nzp1-j)
          i = k + j - 1
          u_vec(k) = -0.5_r2 * dx2 * delz * &
                      u_mat_z(MODULO(i, nxp1), i/nxp1)
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    SELECT CASE (r2)
      CASE (SPr)
        CALL PSDBTRS('N', nxp1*nzp1, nxp1, nxp1, 1, b_mat, 1, desc_b, u_vec, 1, &
                      desc_u, af, b_laf, work, lwork_b_sol, info)
        IF (info /= 0) PRINT*, 'b_finite_PSDBTRS ', info
      CASE (DPr)
        CALL PDDBTRS('N', nxp1*nzp1, nxp1, nxp1, 1, b_mat, 1, desc_b, u_vec, 1, &
                      desc_u, af, b_laf, work, lwork_b_sol, info)
        IF (info /= 0) PRINT*, 'b_finite_PDDBTRS ', info
      CASE DEFAULT
        STOP 'ERROR: Precision selection error - magnetic.f90 finite P*DBTRS'
    END SELECT

    cpcol = 0

    bn = 0.0_r2
    DO j = 0, nxp1*nzp1-1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nxp1*nzp1-j)
          i = k + j - 1
          bn(MODULO(i, nxp1), i/nxp1) = u_vec(k)
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    CALL SLTIMER(8)
    IF (npcol > 1) THEN
      CALL DGSUM2D(ictxt, 'A', ' ', nxp1, nzp1, bn, nxp1, 0, 0)
    END IF
    CALL SLTIMER(8)

    RETURN
  END SUBROUTINE fin_b_poisson

END MODULE magnetic
