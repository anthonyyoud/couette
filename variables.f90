MODULE variables
  USE parameters
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: copy_var, vr_vz

  TYPE, PUBLIC :: VAR
    !variables
    REAL (r2) :: new(0:nx,0:nz)
    REAL (r2) :: old(0:nx,0:nz)
    REAL (r2) :: old2(0:nx,0:nz)
    REAL (r2) :: inter(0:nx,0:nz)
    REAL (r2) :: nlin_new(0:nx,0:nz)
    REAL (r2) :: nlin_old(0:nx,0:nz)
  END TYPE VAR

  TYPE, PUBLIC :: DERIV
    !derivatives
    REAL (r2) :: x(0:nx,0:nz)
    REAL (r2) :: xx(0:nx,0:nz)
    REAL (r2) :: z(0:nx,0:nz)
    REAL (r2) :: zz(0:nx,0:nz)
    REAL (r2) :: zx(0:nx,0:nz)
    REAL (r2) :: zxx(0:nx,0:nz)
    REAL (r2) :: zzz(0:nx,0:nz)
  END TYPE DERIV

  TYPE, PUBLIC :: MAT_COMP
    !matrix components
    REAL (r2) :: lo(2:nx1)
    REAL (r2) :: di(nx1)
    REAL (r2) :: up(nx-2)
  END TYPE MAT_COMP

  TYPE, PUBLIC :: UZ_MAT_COMP
    !matrix components for u in the z-direction
    REAL (r2) :: lo(nz)
    REAL (r2) :: di(0:nz)
    REAL (r2) :: up(0:nz1)
  END TYPE UZ_MAT_COMP

  TYPE, PUBLIC :: ZZ_MAT_COMP
    !matrix components for Z in the z-direction
    REAL (r2) :: lo(2:nz1)
    REAL (r2) :: di(nz1)
    REAL (r2) :: up(nz-2)
  END TYPE ZZ_MAT_COMP

  TYPE (VAR), PUBLIC, SAVE :: ut, zt, psi, bt, jt
  REAL (r2),  PUBLIC, SAVE :: vr(0:nx,0:nz), vz(0:nx,0:nz), &
                              vrold(0:nx,0:nz) = 0.0_r2, &
                              vzold(0:nx,0:nz) = 0.0_r2

  CONTAINS

  SUBROUTINE copy_var(var_out, var_in)
    !Copy a variable
    USE parameters
    IMPLICIT NONE

    REAL (r2), INTENT(IN)  :: var_in(0:nx,0:nz)
    REAL (r2), INTENT(OUT) :: var_out(0:nx,0:nz)

    var_out = var_in

    RETURN
  END SUBROUTINE copy_var

  SUBROUTINE vr_vz(p, vr, vz)
    !Calculate radial and axial velocity components from the stream-function
    USE parameters
    USE ic_bc, ONLY : s
    IMPLICIT NONE

    REAL    (r2), INTENT(IN)  :: p(0:nx,0:nz)
    REAL    (r2), INTENT(OUT) :: vr(0:nx,0:nz), vz(0:nx,0:nz)
    INTEGER (i1)              :: j, k

    DO k = 1, nz1
      vr(:,k) = (-1.0_r2 / (2.0_r2 * s(:) * delz)) * (p(:,k+1) - p(:,k-1))
    END DO

    IF (ABS(tau - 1.0_r2) > EPSILON(tau)) THEN
      vr(:,0) = (-1.0_r2 / (2.0_r2 * s(:) * delz)) * &
                (-3.0_r2 * p(:,0) + 4.0_r2 * p(:,1) - p(:,2))
      vr(:,nz) = (-1.0_r2 / (2.0_r2 * s(:) * delz)) * &
                 (3.0_r2 * p(:,nz) - 4.0_r2 * p(:,nz1) + p(:,nz-2))
    ELSE
      vr(:,0) = 0.0_r2
      vr(:,nz) = 0.0_r2
    END IF

    DO j = 1, nx1
      vz(j,:) = (1.0_r2 / (2.0_r2 * s(j) * delx)) * (p(j+1,:) - p(j-1,:))
    END DO

    vz(0,:) = 0.0_r2 !(1.0_r2 / (2.0_r2 * s(0) * delx)) * &
                  !(-3.0_r2 * p(0,k) + 4.0_r2 * p(1,k) - p(2,k))
    vz(nx,:) = 0.0_r2 !(1.0_r2 / (2.0_r2 * s(nx) * delx)) * &
                  !(3.0_r2 * p(nx,k) - 4.0_r2 * p(nx1,k) + p(nx-2,k))

    RETURN
  END SUBROUTINE vr_vz

END MODULE variables
