MODULE linear
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_rhs_ux, get_rhs_Zx

  CONTAINS

  SUBROUTINE get_rhs_ux(uo, u)
    !Get the linear part for the right-hand side of the tridiagonal system for u
    USE parameters
    USE variables
    USE derivs
    USE ic_bc, ONLY : s
    IMPLICIT NONE

    REAL    (r2),   INTENT(IN)  :: uo(0:nx,0:nz)
    REAL    (r2),   INTENT(OUT) :: u(0:nx,0:nz)
    INTEGER (i1)                :: j, k
    TYPE    (DERIV)             :: du

    CALL deriv_x(uo, du%x)
    CALL deriv_xx(uo, du%xx)   !get derivatives
    CALL deriv_zz(uo, du%zz)

    DO k = 1, nz1
      u(1:nx1,k) = uo(1:nx1,k) + 0.5_r2 * rxx * du%xx(1:nx1,k) + &
                   one_eta * rx * 0.25_r2 * du%x(1:nx1,k) / s(1:nx1) - &
                   one_eta**2 * dt * 0.5_r2 * uo(1:nx1,k) / s(1:nx1)**2 + &
                   0.5_r2 * rzz * du%zz(1:nx1,k)
    END DO

    IF (ABS(tau - 1.0_r2) > EPSILON(tau)) THEN
      u(1:nx1,0) = uo(1:nx1,0) + 0.5_r2 * rxx * du%xx(1:nx1,0) + &
                   one_eta * rx * 0.25_r2 * du%x(1:nx1,0) / s(1:nx1) - &
                   one_eta**2 * dt * 0.5_r2 * uo(1:nx1,0) / s(1:nx1)**2 + &
                   0.5_r2 * rzz * du%zz(1:nx1,0)

      u(1:nx1,nz) = uo(1:nx1,nz) + 0.5_r2 * rxx * du%xx(1:nx1,nz) + &
                    one_eta * rx * 0.25_r2 * du%x(1:nx1,nz) / s(1:nx1) - &
                    one_eta**2 * dt * 0.5_r2 * uo(1:nx1,nz) / s(1:nx1)**2 + &
                    0.5_r2 * rzz * du%zz(1:nx1,nz)
    END IF

    RETURN
  END SUBROUTINE get_rhs_ux

  SUBROUTINE get_rhs_Zx(zo, zn)
    !Get linear part for the right-hand side of the tridiagonal system for Z
    USE parameters
    USE variables
    USE derivs
    USE ic_bc, ONLY : s
    IMPLICIT NONE

    REAL    (r2),   INTENT(IN)  :: zo(0:nx,0:nz)
    REAL    (r2),   INTENT(OUT) :: zn(0:nx,0:nz)
    INTEGER (i1)                :: j, k
    TYPE    (DERIV)             :: dz

    CALL deriv_x(zo, dz%x)
    CALL deriv_xx(zo, dz%xx)   !get derivatives
    CALL deriv_zz(zo, dz%zz)

    DO k = 1, nz1
      zn(1:nx1,k) = zo(1:nx1,k) + 0.5_r2 * rxx * dz%xx(1:nx1,k) + &
                    one_eta * rx * 0.25_r2 * dz%x(1:nx1,k) / s(1:nx1) - &
                    one_eta**2 * dt * 0.5_r2 * zo(1:nx1,k) / s(1:nx1)**2 + &
                    0.5_r2 * rzz * dz%zz(1:nx1,k)
    END DO

    RETURN
  END SUBROUTINE get_rhs_Zx

END MODULE linear
