MODULE nonlinear
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_nlin_ux, get_nlin_Zx

  CONTAINS

  SUBROUTINE get_nlin_ux(uo, uo2, po, po2, bo, bo2, u_nl_n)
    !Get nonlinear part for the right-hand side of the tridiagonal system for u
    USE parameters
    USE variables
    USE derivs
    USE ic_bc, ONLY : s 
    IMPLICIT NONE

    REAL (r2), INTENT(IN)  :: uo(0:nx,0:nz), uo2(0:nx,0:nz), &
                              po(0:nx,0:nz), po2(0:nx,0:nz), &
                              bo(0:nx,0:nz), bo2(0:nx,0:nz)
    REAL (r2), INTENT(OUT) :: u_nl_n(0:nx,0:nz) 
    TYPE (DERIV)           :: du, du2, dp, dp2, db, db2
    INTEGER (i1)           :: j, k  

    CALL deriv_x(uo, du%x)
    CALL deriv_x(uo2, du2%x)
    CALL deriv_z(uo, du%z)
    CALL deriv_z(uo2, du2%z)
    CALL deriv_x(po, dp%x)
    CALL deriv_x(po2, dp2%x)   !get derivatives
    CALL deriv_z(po, dp%z)
    CALL deriv_z(po2, dp2%z)
    CALL deriv_z(bo, db%z)
    CALL deriv_z(bo2, db2%z)

    DO k = 1, nz1
      u_nl_n(1:nx1,k) = (-0.125_r2 * rx / (s(1:nx1) * delz)) * &
                        (3.0_r2 * (dp%x(1:nx1,k) * du%z(1:nx1,k) - &
                        dp%z(1:nx1,k) * du%x(1:nx1,k)) - &
                        (dp2%x(1:nx1,k) * du2%z(1:nx1,k) - &
                        dp2%z(1:nx1,k) * du2%x(1:nx1,k))) + &
                        (one_eta * rz * 0.25_r2 / s(1:nx1)**2) * &
                        (3.0_r2 * uo(1:nx1,k) * dp%z(1:nx1,k) - &
                        uo2(1:nx1,k) * dp2%z(1:nx1,k)) + &
                        0.25_r2 * rz * Q * (3.0_r2 * db%z(1:nx1,k) - &
                        db2%z(1:nx1,k))
    END DO

    IF (ABS(tau - 1.0_r2) > EPSILON(tau)) THEN
      u_nl_n(1:nx1,0) = (-0.125_r2 * rx / (s(1:nx1) * delz)) * &
                        (-3.0_r2 * dp%z(1:nx1,0) * du%x(1:nx1,0) + &
                        dp2%z(1:nx1,0) * du2%x(1:nx1,0)) + &
                        (one_eta * rz * 0.25_r2 / s(1:nx1)**2) * &
                        (3.0_r2 * uo(1:nx1,0) * dp%z(1:nx1,0) - &
                        uo2(1:nx1,0) * dp2%z(1:nx1,0)) + &
                        0.25_r2 * rz * Q * (3.0_r2 * db%z(1:nx1,0) - &
                        db2%z(1:nx1,0))

      u_nl_n(1:nx1,nz) = (-0.125_r2 * rx / (s(1:nx1) * delz)) * &
                         (-3.0_r2 * dp%z(1:nx1,nz) * du%x(1:nx1,nz) + &
                         dp2%z(1:nx1,nz) * du2%x(1:nx1,nz)) + &
                         (one_eta * rz * 0.25_r2 / s(1:nx1)**2) * &
                        (3.0_r2 * uo(1:nx1,nz) * dp%z(1:nx1,nz) - &
                         uo2(1:nx1,nz) * dp2%z(1:nx1,nz)) + &
                         0.25_r2 * rz * Q * (3.0_r2 * db%z(1:nx1,nz) - &
                         db2%z(1:nx1,nz))
    END IF

    RETURN
  END SUBROUTINE get_nlin_ux

  SUBROUTINE get_nlin_Zx(t, uo, uo2, po, po2, zo, zo2, jo, jo2, z_nl_n)
    !Get nonlinear part for the right-hand side of the tridiagonal system for Z
    USE parameters
    USE variables
    USE derivs
    USE ic_bc, ONLY : s
    IMPLICIT NONE

    REAL (r2), INTENT(IN)  :: t, uo(0:nx,0:nz), uo2(0:nx,0:nz), &
                              po(0:nx,0:nz), po2(0:nx,0:nz), &
                              zo(0:nx,0:nz), zo2(0:nx,0:nz), &
                              jo(0:nx,0:nz), jo2(0:nx,0:nz)
    REAL (r2), INTENT(OUT) :: z_nl_n(0:nx,0:nz)
    TYPE (DERIV)           :: du, du2, dp, dp2, dz, dz_2, dj, dj2
    INTEGER (i1)           :: j, k

    CALL deriv_z(uo, du%z)
    CALL deriv_z(uo2, du2%z)
    CALL deriv_x(po, dp%x)
    CALL deriv_x(po2, dp2%x)
    CALL deriv_z(po, dp%z)
    CALL deriv_z(po2, dp2%z)
    CALL deriv_x(zo, dz%x)    !get derivatives
    CALL deriv_x(zo2, dz_2%x)
    CALL deriv_z(zo, dz%z)
    CALL deriv_z(zo2, dz_2%z)
    CALL deriv_z(jo, dj%z)
    CALL deriv_z(jo2, dj2%z)

    DO k = 1, nz1
      z_nl_n(1:nx1,k) = (one_eta * rz * 0.5_r2 / s(1:nx1)) * &
                        (3.0_r2 * uo(1:nx1,k) * du%z(1:nx1,k) - &
                        uo2(1:nx1,k) * du2%z(1:nx1,k)) - &
                        (one_eta * rz * 0.25_r2 / s(1:nx1)**2) * &
                        (3.0_r2 * zo(1:nx1,k) * dp%z(1:nx1,k) - &
                        zo2(1:nx1,k) * dp2%z(1:nx1,k)) - &
                        (0.125_r2 * rx / (s(1:nx1) * delz)) * &
                        ((3.0_r2 * (dp%x(1:nx1,k) * dz%z(1:nx1,k) - &
                        dp%z(1:nx1,k) * dz%x(1:nx1,k))) - &
                        (dp2%x(1:nx1,k) * dz_2%z(1:nx1,k) - &
                        dp2%z(1:nx1,k) * dz_2%x(1:nx1,k))) + &
                        0.25_r2 * rz * Q * (3.0_r2 * dj%z(1:nx1,k) - &
                        dj2%z(1:nx1,k))
    END DO

    RETURN
  END SUBROUTINE get_nlin_Zx

END MODULE nonlinear
