MODULE solve
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: solve_ux, solve_uz, solve_Zx, solve_Zz

  CONTAINS

  SUBROUTINE solve_ux(uo, u, u_nl, t, ux)
    !Solve for the azimuthal velocity field in the x-direction
    USE parameters
    USE variables
    USE ic_bc, ONLY : u_BCS, s
    IMPLICIT NONE

    REAL (r2),       INTENT(IN)    :: t, u(0:nx,0:nz), u_nl(0:nx,0:nz)
    TYPE (MAT_COMP), INTENT(IN)    :: ux
    REAL (r2),       INTENT(INOUT) :: uo(0:nx,0:nz)
    REAL (r2)                      :: ux_rhs(nx1)
    INTEGER (i1)                   :: j, k

    CALL u_BCS(uo, t)

    DO k = 1, nz1
      ux_rhs(:) = u(1:nx1,k) + u_nl(1:nx1,k)   !RHS at fixed z, looping over x

      ux_rhs(1) = ux_rhs(1) + 0.5_r2 * rxx * uo(0,k) - &
                  one_eta * rx * 0.25_r2 * uo(0,k) / s(1)    !BCS
      ux_rhs(nx1) = ux_rhs(nx1) + 0.5_r2 * rxx * uo(nx,k) + &
                    one_eta * rx * 0.25_r2 * uo(nx,k) / s(nx1)

      CALL thomas(xlb, nx1, ux%up, ux%di, ux%lo, ux_rhs)   !Thomas algorithm
                                                           !at each z
      uo(1:nx1,k) = ux_rhs(:)   !intermediate solution
    END DO

    IF (ABS(tau - 1.0_r2) > EPSILON(tau)) THEN !extra entries due to BCS
      DO k = 0, nz, nz
        ux_rhs(:) = u(1:nx1,k) + u_nl(1:nx1,k)

        ux_rhs(1) = ux_rhs(1) + 0.5_r2 * rxx * uo(0,k) - &
                    one_eta * rx * 0.25_r2 * uo(0,k) / s(1)
        ux_rhs(nx1) = ux_rhs(nx1) + 0.5_r2 * rxx * uo(nx,k) + &
                      one_eta * rx * 0.25_r2 * uo(nx,k) / s(nx1)

        CALL thomas(xlb, nx1, ux%up, ux%di, ux%lo, ux_rhs)

        uo(1:nx1,k) = ux_rhs(:)
      END DO
    END IF

    RETURN
  END SUBROUTINE solve_ux

  SUBROUTINE solve_Zx(zo, zn, z_nl, po, t, zx)
    !Solve for the azimuthal vorticity in the x-direction
    USE parameters
    USE variables
    USE ic_bc, ONLY : z_BCS, s
    IMPLICIT NONE

    REAL (r2),       INTENT(IN)    :: t, zn(0:nx,0:nz), po(0:nx,0:nz), &
                                      z_nl(0:nx,0:nz)
    TYPE (MAT_COMP), INTENT(IN)    :: zx
    REAL (r2),       INTENT(INOUT) :: zo(0:nx,0:nz)
    REAL (r2)                      :: zx_rhs(nx1)
    INTEGER (i1)                   :: j, k

    CALL z_BCS(zo, po, t)

    DO k = 1, nz1
      zx_rhs(:) = zn(1:nx1,k) + z_nl(1:nx1,k)   !RHS at fixed z, looping over x

      zx_rhs(1) = zx_rhs(1) + 0.5_r2 * rxx * zo(0,k) - &
                  one_eta * rx * 0.25_r2 * zo(0,k) / s(1)   !BCS
      zx_rhs(nx1) = zx_rhs(nx1) + 0.5_r2 * rxx * zo(nx,k) + &
                    one_eta * rx * 0.25_r2 * zo(nx,k) / s(nx1)

      CALL thomas(xlb, nx1, zx%up, zx%di, zx%lo, zx_rhs)   !Thomas algorithm

      zo(1:nx1,k) = zx_rhs(:)   !intermediate solution
    END DO

    RETURN
  END SUBROUTINE solve_Zx

  SUBROUTINE solve_uz(uo, u, t, uz)
    !Solve for azimuthal velocity in z-direction to give full solution
    USE parameters
    USE variables
    USE ic_bc, ONLY : u_BCS, s
    IMPLICIT NONE

    REAL (r2),          INTENT(IN)    :: t, uo(0:nx,0:nz)
    TYPE (UZ_MAT_COMP), INTENT(IN)    :: uz
    REAL (r2),          INTENT(INOUT) :: u(0:nx,0:nz)
    REAL (r2)                         :: uz_rhs(0:nz), uz_rhs_t1(nz1), &
                                         up(nz-2), di(nz1), lo(2:nz1)
    INTEGER (i1)                      :: j, k

    CALL u_BCS(u, t)

    IF (ABS(tau - 1.0_r2) < EPSILON(tau)) THEN
      up(:) = uz%up(1:nz-2)
      di(:) = uz%di(1:nz1)   !dimensions of upper, lower, diagonal change
      lo(:) = uz%lo(2:nz1)   !for tau=1

      DO j = 1, nx1
        uz_rhs_t1(:) = uo(j,1:nz1)   !RHS at fixed x looping over z

        uz_rhs_t1(1) = uz_rhs_t1(1) + 0.5_r2 * rzz * u(j,0)
        uz_rhs_t1(nz1) = uz_rhs_t1(nz1) + 0.5_r2 * rzz * u(j,nz)   !BCS

        CALL thomas(zlb+1, nz1, up, di, lo, uz_rhs_t1)   !Thomas algorithm

        u(j,1:nz1) = uz_rhs_t1(:)   !full solution
      END DO
    ELSE   !straight-forward if tau /= 1
      DO j = 1, nx1
        uz_rhs(:) = uo(j,:)   !RHS at fixed x looping over z

        CALL thomas(zlb, nz, uz%up, uz%di, uz%lo, uz_rhs)   !Thomas algorithm

        u(j,:) = uz_rhs(:)   !full solution
      END DO
    END IF

    RETURN
  END SUBROUTINE solve_uz

  SUBROUTINE solve_Zz(zo, po, zn, t, zz)
    !Solve for azimuthal vorticity in z-direction to give full solution
    USE parameters
    USE variables
    USE ic_bc, ONLY : z_BCS, s
    IMPLICIT NONE

    REAL (r2),          INTENT(IN)    :: t, zo(0:nx,0:nz), po(0:nx,0:nz)
    TYPE (ZZ_MAT_COMP), INTENT(IN)    :: zz
    REAL (r2),          INTENT(INOUT) :: zn(0:nx,0:nz)
    REAL (r2)                         :: Zz_rhs(nz1)
    INTEGER (i1)                      :: j, k

    CALL z_BCS(zn, po, t)

    DO j = 1, nx1
      Zz_rhs(:) = zo(j,1:nz1)   !RHS at fixed x, looping over z

      Zz_rhs(1) = Zz_rhs(1) + 0.5_r2 * rzz * zn(j,0)
      Zz_rhs(nz1) = Zz_rhs(nz1) + 0.5_r2 * rzz * zn(j,nz)   !BCS

      CALL thomas(zlb+1, nz1, zz%up, zz%di, zz%lo, Zz_rhs)   !Thomas algorithm

      zn(j,1:nz1) = Zz_rhs(:)   !full solution
    END DO

    RETURN
  END SUBROUTINE solve_Zz

  SUBROUTINE thomas (lb, m, up, di, lo, r)
    !Thomas algorithm for solving a tridiagonal system of equations
    USE parameters, ONLY : i1, r2
    IMPLICIT NONE

    INTEGER (i1)                :: j
    INTEGER (i1), INTENT(IN)    :: m, lb  !lb is the vector lower-bound
    REAL (r2),    INTENT(IN)    :: up(lb:m-1), di(lb:m), lo(lb+1:m)
    REAL (r2),    INTENT(INOUT) :: r(lb:m)
    REAL (r2)                   :: dnew(lb:m), aa = 0.0_r2

    dnew = di   !create new diagonal so as not to destroy original
    DO j = lb+1, m
      aa = -lo(j) / dnew(j-1)
      dnew(j) = dnew(j) + aa * up(j-1)   !eliminate lower-diagonal
      r(j) = r(j) + aa * r(j-1)
    END DO

    r(m) = r(m) / dnew(m)   !first step back-substitution

    DO j = m-1, lb, -1
      r(j) = (r(j) - up(j) * r(j+1)) / dnew(j)  !solve by back-substitution
    END DO

    RETURN
  END SUBROUTINE thomas

END MODULE solve
