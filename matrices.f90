MODULE matrices
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: matrix_setup, psi_mat_setup, b_mat_setup, fin_b_mat_setup, &
            j_mat_setup, fin_j_mat_setup

  CONTAINS

  SUBROUTINE matrix_setup(ux, uz, zx, zz)
    !Set up the upper, lower and diagonal parts of the tridiagonal matrices
    !for the solution of the azimuthal velocity and vorticity equations
    USE parameters
    USE variables
    USE ic_bc, ONLY : s
    IMPLICIT NONE

    TYPE (MAT_COMP),    INTENT(OUT) :: ux, zx
    TYPE (UZ_MAT_COMP), INTENT(OUT) :: uz
    TYPE (ZZ_MAT_COMP), INTENT(OUT) :: zz

    ux%di(:) = 1.0_r2 + rxx + one_eta**2 * dt * 0.5_r2 / s(1:nx-1)**2
    ux%lo(:) = -0.5_r2 * rxx + one_eta * rx * 0.25_r2 / s(2:nx-1)
    ux%up(:) = -0.5_r2 * rxx - one_eta * rx * 0.25_r2 / s(1:nx-2)

    uz%di(:) = 1.0_r2 + rzz
    uz%lo(:) = -0.5_r2 * rzz
    uz%up(:) = -0.5_r2 * rzz

    IF (tau /= 1) THEN
      uz%di(0) = 1.0_r2 + rzz + (rz * tau / (1.0_r2 - tau))  !extra entries -
      uz%di(nz) = 1.0_r2 + rzz + (rz * tau / (1.0_r2 - tau)) !Neumann BCS ends
      uz%lo(nz) = -rzz
      uz%up(0) = -rzz
    END IF

    zx%di(:) = 1.0_r2 + rxx + one_eta**2 * dt * 0.5_r2 / s(1:nx-1)**2
    zx%lo(:) = -0.5_r2 * rxx + one_eta * rx * 0.25_r2 / s(2:nx-1)
    zx%up(:) = -0.5_r2 * rxx - one_eta * rx * 0.25_r2 / s(1:nx-2)

    zz%di(:) = 1.0_r2 + rzz
    zz%lo(:) = -0.5_r2 * rzz
    zz%up(:) = -0.5_r2 * rzz

    RETURN
  END SUBROUTINE matrix_setup

  SUBROUTINE psi_mat_setup(p_mat, desc_p, af)
    !Setup of LHS matrix in solution of stream-function Poisson equation
    USE parameters
    USE ic_bc, ONLY : s
    IMPLICIT NONE

    INTEGER (i1), INTENT(OUT) :: desc_p(7)
    REAL (r2),    INTENT(OUT) :: p_mat(p_M,p_N), af(laf)
    REAL (r2)                 :: alp(0:nx), gam(0:nx), beta, delta, &
                                 work(lwork_fac)
    INTEGER (i1)              :: i, j, k, info, cpcol

    alp(:) = dz2 + 0.5_r2 * delx * dz2 * one_eta / s(:)   !coefficients
    gam(:) = dz2 - 0.5_r2 * delx * dz2 * one_eta / s(:)   !in matrix
    beta = -2.0_r2 * (dz2 + dx2)
    delta = dx2

    desc_p(1) = 501
    desc_p(2) = ictxt
    desc_p(3) = nx1*nz1   !ScaLAPACK descriptor vector for LHS matrix
    desc_p(4) = nb
    desc_p(5) = 0
    desc_p(6) = 2*nx1+1

    cpcol = 0

    !distribute the matrix over the process grid using the block-column
    !distribution for banded matrices
    DO j = 1, nx1*nz1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nx1*nz1-j+1)
          i = k + j - 1
          p_mat(nx,k) = beta   !diagonal
          IF (i /= 1) THEN
            p_mat(nx1,k) = gam(MODULO(i-1, nx1))   !upper-diagonal
            IF(MODULO(i-1, nx1) == 0) THEN
              p_mat(nx1,k) = 0.0_r2   !upper-diagonal, BCS
            END IF
          END IF
          IF (i /= nx1*nz1) THEN
            p_mat(nxp1,k) = alp(MODULO(i, nx1) + 1) !lower-diagonal
            IF (MODULO(i, nx1) == 0) THEN
              p_mat(nxp1,k) = 0.0_r2   !lower-diagonal, BCS
            END IF
          END IF
          IF (i > nx1) THEN
            p_mat(1,k) = delta  !upper band
          END IF
          IF (i <= nx1*nz1-nx1) THEN
            p_mat(2*nx1+1,k) = delta   !lower band
          END IF
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    !LU factorisation of the matrix for use in PD/PSDBTRS to solve
    SELECT CASE (r2)
      CASE (SPr)
        CALL PSDBTRF(nx1*nz1, nx1, nx1, p_mat, 1, desc_p, af, laf, &
                     work, lwork_fac, info)
        IF (info /= 0) PRINT*, 'psi_PSDBTRF ', info
      CASE (DPr)
        CALL PDDBTRF(nx1*nz1, nx1, nx1, p_mat, 1, desc_p, af, laf, &
                     work, lwork_fac, info)
        IF (info /= 0) PRINT*, 'psi_PDDBTRF ', info
      CASE DEFAULT
        STOP 'ERROR: Precision selection error - matrices.f90 psi_P*DBTRF'
    END SELECT

    RETURN
  END SUBROUTINE psi_mat_setup

  SUBROUTINE b_mat_setup(b_mat, desc_b, af)
    !Setup of LHS matrix in solution of magnetic Poisson equation.
    !Algorithm as for stream-function above.
    USE parameters
    USE ic_bc, ONLY : s
    IMPLICIT NONE

    INTEGER (i1), INTENT(OUT) :: desc_b(7)
    REAL (r2),    INTENT(OUT) :: b_mat(b_M,b_N), af(b_laf)
    REAL (r2)                 :: alp(0:nx), beta(0:nx), gam(0:nx), delta, &
                                 work(lwork_b_fac)
    INTEGER (i1)              :: i, j, k, info, cpcol

    alp(:) = dz2 - 0.5_r2 * delx * dz2 * one_eta / s(:)
    beta(:) = -2.0_r2 * (dz2 + dx2) - dx2 * dz2 * one_eta**2 / s(:)**2
    gam(:) = dz2 + 0.5_r2 * delx * dz2 * one_eta / s(:)
    delta = dx2

    desc_b(1) = 501
    desc_b(2) = ictxt
    desc_b(3) = nxp1*nz1
    desc_b(4) = nb
    desc_b(5) = 0
    desc_b(6) = 2*nxp1+1

    cpcol = 0

    DO j = 0, nxp1*nz1-1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nxp1*nz1-j)
          i = k + j - 1
          b_mat(nx+2,k) = beta(MODULO(i, nxp1))   !diagonal
          IF (MODULO(i, nxp1) == 0) THEN
            b_mat(nx+2,k) = (2.0_r2 * alp(0) * delx * one_eta / &
                             s(0)) + beta(0)   !diagonal, BCS
          END IF
          IF (MODULO(i+1, nxp1) == 0) THEN
            b_mat(nx+2,k) = (-2.0_r2 * gam(nx) * delx * one_eta / &
                             s(nx)) + beta(nx) !diagonal, BCS
          END IF
          IF (i /= 0) THEN
            b_mat(nxp1,k) = gam(MODULO(i-1, nxp1))   !upper-diagonal
            IF(MODULO(i, nxp1) == 0) THEN
              b_mat(nxp1,k) = 0.0_r2   !upper-diagonal, BCS
            END IF
            IF(MODULO(i-1, nxp1) == 0) THEN
              b_mat(nxp1,k) = alp(0) + gam(0)   !upper-diagonal, BCS
            END IF
          END IF
          IF (i /= nxp1*nz1-1) THEN
            IF (MODULO(i+1, nxp1) /= 0) THEN
              b_mat(nx+3,k) = alp(MODULO(i+1,nxp1))   !lower-diagonal
            END IF
            IF (MODULO(i+1, nxp1) == 0) THEN
              b_mat(nx+3,k) = 0.0_r2   !lower-diagonal, BCS
            END IF
            IF (MODULO(i+2, nxp1) == 0) THEN
              b_mat(nx+3,k) = alp(nx) + gam(nx)   !lower-diagonal, BCS
            END IF
          END IF
          IF (i > nx) THEN
            b_mat(1,k) = delta   !upper band
          END IF
          IF (i < nxp1*nz1-nxp1) THEN
            b_mat(2*nxp1+1,k) = delta   !lower band
          END IF
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    SELECT CASE (r2)
      CASE (SPr)
        CALL PSDBTRF(nxp1*nz1, nxp1, nxp1, b_mat, 1, desc_b, af, b_laf, &
                     work, lwork_b_fac, info)
        IF (info /= 0) PRINT*, 'b_infinite_PSDBTRF ', info
      CASE (DPr)
        CALL PDDBTRF(nxp1*nz1, nxp1, nxp1, b_mat, 1, desc_b, af, b_laf, &
                     work, lwork_b_fac, info)
        IF (info /= 0) PRINT*, 'b_infinite_PDDBTRF ', info
      CASE DEFAULT
        STOP 'ERROR: Precision selection error - matrices.f90 b_inf_P*DBTRF'
    END SELECT

    RETURN
  END SUBROUTINE b_mat_setup

  SUBROUTINE fin_b_mat_setup(b_mat, desc_b, af)
    USE parameters
    USE ic_bc, ONLY : s
    IMPLICIT NONE

    INTEGER (i1), INTENT(OUT) :: desc_b(7)
    REAL (r2),    INTENT(OUT) :: b_mat(b_M,b_N), af(b_laf)
    REAL (r2)                 :: alp(0:nx), beta(0:nx), gam(0:nx), delta, &
                                 work(lwork_b_fac)
    INTEGER (i1)              :: i, j, k, info, cpcol

    alp(:) = dz2 - 0.5_r2 * delx * dz2 * one_eta / s(:)
    beta(:) = -2.0_r2 * (dz2 + dx2) - dx2 * dz2 * one_eta**2 / s(:)**2
    gam(:) = dz2 + 0.5_r2 * delx * dz2 * one_eta / s(:)
    delta = dx2

    desc_b(1) = 501
    desc_b(2) = ictxt
    desc_b(3) = nxp1*nzp1
    desc_b(4) = nb
    desc_b(5) = 0
    desc_b(6) = 2*nxp1+1

    cpcol = 0

    DO j = 0, nxp1*nzp1-1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nxp1*nzp1-j)
          i = k + j - 1
          !diagonal
          b_mat(nx+2,k) = beta(MODULO(i, nxp1))
          !diagonal, j=0
          IF (MODULO(i, nxp1) == 0) THEN
            b_mat(nx+2,k) = (2.0_r2 * alp(0) * delx * one_eta / &
                             s(0)) + beta(0)
          END IF
          !diagonal, k=0
          IF (i < nxp1) THEN
            b_mat(nx+2,k) = beta(MODULO(i, nxp1)) - &
                            2.0_r2 * delta * delz * (1.0_r2 - tau) / tau
          END IF
          !diagonal, j=nx
          IF (MODULO(i+1, nxp1) == 0) THEN
            b_mat(nx+2,k) = (-2.0_r2 * gam(nx) * delx * one_eta / &
                             s(nx)) + beta(nx)
          END IF
          !diagonal, k=nz
          IF (i > nxp1*nzp1-nxp1) THEN
            b_mat(nx+2,k) = beta(MODULO(i, nxp1)) - &
                            2.0_r2 * delta * delz * (1.0_r2 - tau) / tau
          END IF
          !diagonal, j=0, k=0
          IF (i == 0) THEN
            b_mat(nx+2,k) = beta(0) - &
                            2.0_r2 * delta * delz * (1.0_r2 - tau) / tau + &
                           (2.0_r2 * alp(0) * delx * one_eta / s(0))
          END IF
          !diagonal, j=0, k=nz
          IF (i == nxp1*nzp1-nxp1) THEN
            b_mat(nx+2,k) = beta(0) - &
                            2.0_r2 * delta * delz * (1.0_r2 - tau) / tau + &
                           (2.0_r2 * alp(0) * delx * one_eta / s(0))
          END IF
          !diagonal, j=nx, k=0
          IF (i == nx) THEN
            b_mat(nx+2,k) = beta(nx) - &
                            2.0_r2 * delta * delz * (1.0_r2 - tau) / tau - &
                           (2.0_r2 * gam(nx) * delx * one_eta / s(nx))
          END IF
          !diagonal, j=nx, k=nz
          IF (i == nxp1*nzp1-1) THEN
            b_mat(nx+2,k) = beta(nx) - &
                            2.0_r2 * delta * delz * (1.0_r2 - tau) / tau - &
                           (2.0_r2 * gam(nx) * delx * one_eta / s(nx))
          END IF
          IF (i /= 0) THEN
            !super-diagonal
            b_mat(nxp1,k) = gam(MODULO(i-1, nxp1))
            !super-diagonal, j=nx, so j=0 not present
            IF(MODULO(i, nxp1) == 0) THEN
              b_mat(nxp1,k) = 0.0_r2
            END IF
            !super-diagonal, j=0
            IF(MODULO(i-1, nxp1) == 0) THEN
              b_mat(nxp1,k) = alp(0) + gam(0)
            END IF
          END IF
          IF (i /= nxp1*nzp1-1) THEN
            !sub-diagonal
            IF (MODULO(i+1, nxp1) /= 0) THEN
              b_mat(nx+3,k) = alp(MODULO(i+1, nxp1))
            END IF
            !sub-diagonal, j=0, so j=nx not present
            IF (MODULO(i+1, nxp1) == 0) THEN
              b_mat(nx+3,k) = 0.0_r2
            END IF
            !sub-diagonal, j=nx
            IF (MODULO(i+2, nxp1) == 0) THEN
              b_mat(nx+3,k) = alp(nx) + gam(nx)
            END IF
          END IF
          !upper branch
          IF (i > nx) THEN
            b_mat(1,k) = delta
            IF (i < 2*nxp1) THEN
              b_mat(1,k) = 2.0_r2 * delta
            END IF
          END IF
          !lower branch
          IF (i < nxp1*nzp1-nxp1) THEN
            b_mat(2*nxp1+1,k) = delta
            IF (i >= nxp1*nzp1-2*nxp1) THEN
              b_mat(2*nxp1+1,k) = 2.0_r2 * delta
            END IF
          END IF
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    SELECT CASE (r2)
      CASE (SPr)
        CALL PSDBTRF(nxp1*nzp1, nxp1, nxp1, b_mat, 1, desc_b, af, b_laf, &
                     work, lwork_b_fac, info)
        IF (info /= 0) PRINT*, 'b_finite_PSDBTRF ', info
      CASE (DPr)
        CALL PDDBTRF(nxp1*nzp1, nxp1, nxp1, b_mat, 1, desc_b, af, b_laf, &
                     work, lwork_b_fac, info)
        IF (info /= 0) PRINT*, 'b_finite_PDDBTRF ', info
      CASE DEFAULT
        STOP 'ERROR: Precision selection error - matrices.f90 b_fin_P*DBTRF'
    END SELECT

    RETURN
  END SUBROUTINE fin_b_mat_setup

  SUBROUTINE j_mat_setup(j_mat, desc_j, af)
    USE parameters
    USE ic_bc, ONLY : s
    IMPLICIT NONE

    INTEGER (i1), INTENT(OUT) :: desc_j(7)
    REAL (r2),    INTENT(OUT) :: j_mat(j_M,j_N), af(laf)
    REAL (r2)                 :: alp(0:nx), beta(0:nx), gam(0:nx), delta, &
                                 work(lwork_fac)
    INTEGER (i1)              :: i, j, k, info, cpcol

    alp(:) = dz2 - 0.5_r2 * delx * dz2 * one_eta / s(:)
    beta(:) = -2.0_r2 * (dz2 + dx2) - dx2 * dz2 * one_eta**2 / s(:)**2
    gam(:) = dz2 + 0.5_r2 * delx * dz2 * one_eta / s(:)
    delta = dx2

    desc_j(1) = 501
    desc_j(2) = ictxt
    desc_j(3) = nx1*nzp1
    desc_j(4) = nb
    desc_j(5) = 0
    desc_j(6) = 2*nx1+1

    cpcol = 0

    DO j = 1, nx1*nzp1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nx1*nzp1-j+1)
          i = k + j - 1
          j_mat(nx,k) = beta(MODULO(i-1, nx1) + 1)   !diagonal
          IF (i <= nx1) THEN
            j_mat(nx,k) = beta(MODULO(i-1, nx1) + 1) - &  !extra diagonal due
                          2.0_r2 * delta * delz * tau / (1.0_r2 - tau) !to tau
          END IF
          IF (i >= nx1*nzp1-nx1+1) THEN
            j_mat(nx,k) = beta(MODULO(i-1, nx1) + 1) - &  !extra diagonal due
                          2.0_r2 * delta * delz * tau / (1.0_r2 - tau) !to tau
          END IF
          IF (i /= 1) THEN
            j_mat(nx1,k) = gam(MODULO(i-1, nx1))   !upper-diagonal
            IF(MODULO(i-1, nx1) == 0) THEN
              j_mat(nx1,k) = 0.0_r2  !upper-diagonal, BCS
            END IF
          END IF
          IF (i /= nx1*nzp1) THEN
            j_mat(nxp1,k) = alp(MODULO(i,nx1) + 1)   !lower-diagonal
            IF (MODULO(i, nx1) == 0) THEN
              j_mat(nxp1,k) = 0.0_r2   !lower-diagonal, BCS
            END IF
          END IF
          IF (i > nx1) THEN
            j_mat(1,k) = delta   !upper band
            IF (i <= 2*nx1) THEN
              j_mat(1,k) = 2.0_r2 * delta   !upper band, BCS
            END IF
          END IF
          IF (i <= nx1*nzp1-nx1) THEN
            j_mat(2*nx1+1,k) = delta   !lower band
            IF (i > nx1*nzp1-2*nx1) THEN
              j_mat(2*nx1+1,k) = 2.0_r2 * delta   !lower band, BCS
            END IF
          END IF
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    SELECT CASE (r2)
      CASE (SPr)
        CALL PSDBTRF(nx1*nzp1, nx1, nx1, j_mat, 1, desc_j, af, laf, &
                     work, lwork_fac, info)
        IF (info /= 0) PRINT*, 'j_infinite_PSDBTRF ', info
      CASE (DPr)
        CALL PDDBTRF(nx1*nzp1, nx1, nx1, j_mat, 1, desc_j, af, laf, &
                     work, lwork_fac, info)
        IF (info /= 0) PRINT*, 'j_infinite_PDDBTRF ', info
      CASE DEFAULT
        STOP 'ERROR: Precision selection error - matrices.f90 j_inf_P*DBTRF'
    END SELECT

    RETURN
  END SUBROUTINE j_mat_setup

  SUBROUTINE fin_j_mat_setup(j_mat, desc_j, af)
    USE parameters
    USE ic_bc, ONLY : s
    IMPLICIT NONE

    INTEGER (i1), INTENT(OUT) :: desc_j(7)
    REAL (r2),    INTENT(OUT) :: j_mat(j_M,j_N), af(laf)
    REAL (r2)                 :: alp(0:nx), beta(0:nx), gam(0:nx), delta, &
                                 work(lwork_fac)
    INTEGER (i1)              :: i, j, k, info, cpcol

    alp(:) = dz2 - 0.5_r2 * delx * dz2 * one_eta / s(:)
    beta(:) = -2.0_r2 * (dz2 + dx2) - dx2 * dz2 * one_eta**2 / s(:)**2
    gam(:) = dz2 + 0.5_r2 * delx * dz2 * one_eta / s(:)
    delta = dx2

    desc_j(1) = 501
    desc_j(2) = ictxt
    desc_j(3) = nx1*nz1
    desc_j(4) = nb
    desc_j(5) = 0
    desc_j(6) = 2*nx1+1

    cpcol = 0

    DO j = 1, nx1*nz1, nb
      IF (mycol == cpcol) THEN
        DO k = 1, MIN(nb, nx1*nz1-j+1)
          i = k + j - 1
          j_mat(nx,k) = beta(MODULO(i-1, nx1) + 1)   !diagonal
          IF (i /= 1) THEN
            j_mat(nx1,k) = gam(MODULO(i-1, nx1))   !upper-diagonal
            IF(MODULO(i-1, nx1) == 0) THEN
              j_mat(nx1,k) = 0.0_r2   !upper-diagonal, BCS
            END IF
          END IF
          IF (i /= nx1*nz1) THEN
            j_mat(nxp1,k) = alp(MODULO(i, nx1) + 1)   !lower-diagonal
            IF (MODULO(i, nx1) == 0) THEN
              j_mat(nxp1,k) = 0.0_r2   !lower-diagonal, BCS
            END IF
          END IF
          IF (i > nx1) THEN
            j_mat(1,k) = delta   !upper band
          END IF
          IF (i <= nx1*nz1-nx1) THEN
            j_mat(2*nx1+1,k) = delta   !lower band
          END IF
        END DO
      END IF
      IF (cpcol == npcol) EXIT
      cpcol = cpcol + 1
    END DO

    SELECT CASE (r2)
      CASE (SPr)
        CALL PSDBTRF(nx1*nz1, nx1, nx1, j_mat, 1, desc_j, af, laf, &
                     work, lwork_fac, info)
        IF (info /= 0) PRINT*, 'j_infinite_PSDBTRF ', info
      CASE (DPr)
        CALL PDDBTRF(nx1*nz1, nx1, nx1, j_mat, 1, desc_j, af, laf, &
                     work, lwork_fac, info)
        IF (info /= 0) PRINT*, 'j_infinite_PDDBTRF ', info
      CASE DEFAULT
        STOP 'ERROR: Precision selection error - matrices.f90 j_inf_P*DBTRF'
    END SELECT

    RETURN
  END SUBROUTINE fin_j_mat_setup

END MODULE matrices
