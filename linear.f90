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

    REAL (r2), INTENT(IN)  :: uo(0:nx,0:nz)
    REAL (r2), INTENT(OUT) :: u(0:nx,0:nz)
    TYPE (DERIV)           :: du
    INTEGER (i1)           :: j, k

    CALL deriv_x(uo, du%x)
    CALL deriv_xx(uo, du%xx)   !get derivatives
    CALL deriv_zz(uo, du%zz)

    DO k = 1, nz1
      u(1:nx1,k) = uo(1:nx1,k) + (0.5_r2 * rxx * du%xx(1:nx1,k)) + &
                   (((1.0_r2 - eta) * rx) / (4.0_r2 * s(1:nx1))) * du%x(1:nx1,k) - &
                   (((1.0_r2 - eta)**2 * dt) / (2.0_r2 * s(1:nx1)**2)) * &
                   uo(1:nx1,k) + 0.5_r2 * rzz * du%zz(1:nx1,k)
    END DO

    IF (tau /= 1) THEN
      u(1:nx1,0) = uo(1:nx1,0) + (0.5_r2 * rxx * du%xx(1:nx1,0)) + &
                   (((1.0_r2 - eta) * rx) / (4.0_r2 * s(1:nx1))) * du%x(1:nx1,0) - &
                   (((1.0_r2 - eta)**2 * dt) / (2.0_r2 * s(1:nx1)**2)) * &
                   uo(1:nx1,0) + 0.5_r2 * rzz * du%zz(1:nx1,0)

      u(1:nx1,nz) = uo(1:nx1,nz) + (0.5_r2 * rxx * du%xx(1:nx1,nz)) + &
                    (((1.0_r2 - eta) * rx) / (4.0_r2 * s(1:nx1))) * du%x(1:nx1,nz) - &
                    (((1.0_r2 - eta)**2 * dt) / (2.0_r2 * s(1:nx1)**2)) * &
                    uo(1:nx1,nz) + 0.5_r2 * rzz * du%zz(1:nx1,nz)
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

    REAL (r2), INTENT(IN)  :: zo(0:nx,0:nz)
    REAL (r2), INTENT(OUT) :: zn(0:nx,0:nz)
    TYPE (DERIV)           :: dz
    INTEGER (i1)           :: j, k

    CALL deriv_x(zo, dz%x)
    CALL deriv_xx(zo, dz%xx)   !get derivatives
    CALL deriv_zz(zo, dz%zz)

    DO k = 1, nz1
      zn(1:nx1,k) = zo(1:nx1,k) + (0.5_r2 * rxx * dz%xx(1:nx1,k)) + &
                    (((1.0_r2 - eta) * rx) / (4.0_r2 * s(1:nx1))) * dz%x(1:nx1,k) - &
                    (((1.0_r2 - eta)**2 * dt) / (2.0_r2 * s(1:nx1)**2)) * &
                    zo(1:nx1,k) + 0.5_r2 * rzz * dz%zz(1:nx1,k)
    END DO

    RETURN
  END SUBROUTINE get_rhs_Zx

END MODULE linear
