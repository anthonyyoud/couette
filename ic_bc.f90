MODULE ic_bc
  !Routines to do with grid setup, initial and boundary conditions.
  USE parameters
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_xzs, ICS, u_BCS, z_BCS, p_BCS, b_BCS, j_BCS

  REAL (r2), PUBLIC :: x(0:nx), x_(0:nx), th(0:nt), z(0:nz), s(0:nx) 
                                                !finite-difference mesh
                                                !s=eta+(1-eta)*x
  CONTAINS

  SUBROUTINE get_xzs()
    !Finite-difference mesh
    USE parameters
    IMPLICIT NONE

    INTEGER (i1) :: j, k, l

    DO k = 0, nz
      z(k) = REAL(k,r2) * delz
    END DO

    DO l = 0, nt
      th(l) = REAL(l,r2) * delt
    END DO

    DO j = 0, nx
      x(j) = REAL(j,r2) * delx
    END DO

    x_ = x + 1.0_r2                 !shift radial coordinate for OpenDX
    s = eta + one_eta * x

    RETURN
  END SUBROUTINE get_xzs

  SUBROUTINE ICS(u, zn, pn, bn, jn, p)
    !Initial conditions
    USE parameters
    IMPLICIT NONE

    INTEGER (i1), INTENT(OUT)   :: p
    REAL    (r2), INTENT(INOUT) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                   pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                   jn(0:nx,0:nz)
    INTEGER (i1)                :: j, k
    LOGICAL                     :: state_exist

    IF (restart) THEN
      INQUIRE(FILE='end_state.dat', EXIST=state_exist)
      !exit if doing restart but end_state.dat does not exist
      IF (.NOT. state_exist) STOP 'ERROR: restart=.TRUE.&
                                   &but end_state.dat does not exist.'
      PRINT*, 'Getting restart conditions'
      CALL state_restart(u, zn, pn, bn, jn, p)  !get saved data if restart
    ELSE                   !put in estimate of eigen-function shape
      IF (ABS(tau - 1.0_r2) < EPSILON(tau)) THEN !which satisfies BCS
                                                 !multiplied by small seed
        DO k = 0, nz
          u(:,k) = seed * SIN(2.0_r2*pi*z(k)/gamma) * SIN(pi*x(:))
          pn(:,k) = seed * SIN(2.0_r2*pi*z(k)/gamma) * SIN(pi*x(:))
          bn(:,k) = seed * COS(2.0_r2*pi*z(k)/gamma) / s(:)
          jn(:,k) = seed * SIN(pi*x(:)) * SIN(2.0_r2*pi*z(k)/gamma)
        END DO
      ELSE
        DO k = 0, nz
          u(:,k) = seed * SIN(pi*x(:)) * COS(2.0_r2*pi*z(k)/gamma)
          pn(:,k) = seed * SIN(pi*x(:)) * SIN(2.0_r2*pi*z(k)/gamma)
          bn(:,k) = seed * SIN(2.0_r2*pi*z(k)/gamma) / s(:)
          jn(:,k) = seed * SIN(pi*x(:)) * COS(2.0_r2*pi*z(k)/gamma)
        END DO

        DO k = 1, nz1
          DO j = 1, nx1
            zn(j,k) = -(pn(j+1,k) - 2.0_r2 * pn(j,k) + pn(j-1,k)) / &
                       (s(j) * dx2) + &
                       0.5_r2 * one_eta * (pn(j+1,k) - pn(j-1,k)) / &
                       (s(j)**2 * delx) - &
                       (pn(j,k+1) - 2.0_r2 * pn(j,k) + pn(j,k-1)) / &
                       (s(j) * dz2)  !ICS based on above fields for vorticity
          END DO
        END DO
      END IF
    END IF

    RETURN
  END SUBROUTINE ICS

  SUBROUTINE state_restart(u, zn, pn, bn, jn, p)
    !Get restart data
    USE parameters
    IMPLICIT NONE

    INTEGER (i1), INTENT(OUT) :: p
    REAL    (r2), INTENT(OUT) :: u(0:nx,0:nz), zn(0:nx,0:nz), pn(0:nx,0:nz), &
                                 bn(0:nx,0:nz), jn(0:nx,0:nz)
    INTEGER (i1)              :: j, k, nx_prev, nz_prev, alloc_err
    REAL    (r2)              :: dt_prev
    REAL    (r2), ALLOCATABLE :: u_prev(:,:), z_prev(:,:), p_prev(:,:), &
                                 b_prev(:,:), j_prev(:,:)

    OPEN (50, FILE = 'end_state.dat', FORM='unformatted')

    READ (50) nx_prev
    READ (50) nz_prev
    READ (50) p
    READ (50) dt_prev

    IF ((nx_prev /= nx) .OR. (nz_prev /= nz)) THEN  !interpolate onto new grid
      PRINT*, 'Interpolating onto new grid...'
      ALLOCATE(u_prev(0:nx_prev,0:nz_prev), z_prev(0:nx_prev,0:nz_prev), &
               p_prev(0:nx_prev,0:nz_prev), b_prev(0:nx_prev,0:nz_prev) , &
               j_prev(0:nx_prev,0:nz_prev), STAT=alloc_err)
      IF (alloc_err /= 0) STOP 'ERROR: Interpolating allocation error' 
      READ (50) u_prev
      READ (50) z_prev
      READ (50) p_prev
      READ (50) b_prev
      READ (50) j_prev
      CALL inter(u_prev, z_prev, p_prev, b_prev, j_prev, nx_prev, nz_prev, &
                 u, zn, pn, bn, jn)
      IF (ALLOCATED(u_prev)) DEALLOCATE(u_prev, STAT=alloc_err)
      IF (alloc_err /= 0) STOP 'ERROR: Interpolating deallocation error'
      IF (ALLOCATED(z_prev)) DEALLOCATE(z_prev, STAT=alloc_err)
      IF (alloc_err /= 0) STOP 'ERROR: Interpolating deallocation error'
      IF (ALLOCATED(p_prev)) DEALLOCATE(p_prev, STAT=alloc_err)
      IF (alloc_err /= 0) STOP 'ERROR: Interpolating deallocation error'
      IF (ALLOCATED(b_prev)) DEALLOCATE(b_prev, STAT=alloc_err)
      IF (alloc_err /= 0) STOP 'ERROR: Interpolating deallocation error'
      IF (ALLOCATED(j_prev)) DEALLOCATE(j_prev, STAT=alloc_err)
      IF (alloc_err /= 0) STOP 'ERROR: Interpolating deallocation error'
    ELSE     !just read the data
      READ (50) u
      READ (50) zn
      READ (50) pn
      READ (50) bn
      READ (50) jn
    END IF

    CLOSE (50)   

    IF (dt_prev /= dt) THEN  !if restart tstep /= old tstep adjust index 'p'
      p = p * dt_prev / dt
    END IF

    RETURN
  END SUBROUTINE state_restart

  SUBROUTINE inter(u_, z_, p_, b_, j_, nxp, nzp, &
                   u, zn, pn, bn, jn)
    !Setup for interpolating fields from smaller grid onto larger grid
    USE parameters
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN)  :: nxp, nzp
    REAL    (r2), INTENT(IN)  :: u_(0:nxp,0:nzp), z_(0:nxp,0:nzp), &
                                 p_(0:nxp,0:nzp), b_(0:nxp,0:nzp), &
                                 j_(0:nxp,0:nzp)
    REAL    (r2), INTENT(OUT) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                 pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                 jn(0:nx,0:nz)
    INTEGER (i1)              :: j, k
    REAL    (r2)              :: dx_prev, dz_prev, x_prev(0:nxp), z_prev(0:nzp)

    dx_prev = 1.0_r2 / nxp       !previous space mesh
    dz_prev = gamma / nzp

    DO j = 0, nxp
      x_prev(j) = REAL(j,r2) * dx_prev
    END DO
                                  !previous coordinates
    DO k = 0, nzp
      z_prev(k) = REAL(k,r2) * dz_prev
    END DO

    CALL inter_var(u_, x_prev, z_prev, nxp, nzp, u)
    CALL inter_var(z_, x_prev, z_prev, nxp, nzp, zn)
    CALL inter_var(p_, x_prev, z_prev, nxp, nzp, pn)    !interpolate
    CALL inter_var(b_, x_prev, z_prev, nxp, nzp, bn)
    CALL inter_var(j_, x_prev, z_prev, nxp, nzp, jn)

    RETURN
  END SUBROUTINE inter

  SUBROUTINE inter_var(in_var, x_prev, z_prev, nxp, nzp, out_var)
    !Bilinearly interpolate the actual fields onto a new grid
    !See Numerical Recipes in Fortran77 Chap. 3.6 p.116
    USE parameters
    IMPLICIT NONE
    
    INTEGER (i1), INTENT(IN)  :: nxp, nzp
    REAL    (r2), INTENT(IN)  :: in_var(0:nxp,0:nzp), x_prev(0:nxp), &
                                 z_prev(0:nzp)
    REAL    (r2), INTENT(OUT) :: out_var(0:nx,0:nz)
    INTEGER (i1)              :: j, k, j2, k2
    REAL    (r2)              :: int1, int2

    DO k = 0, nz   !new 'z' index
      k2 = INT(nzp*k/nz,i1)   !old 'z' index
      DO j = 0, nx   !new 'x' index
        j2 = INT(nxp*j/nx,i1)   !old 'x' index
        int1 = (x(j) - x_prev(j2)) / (x_prev(j2+1) - x_prev(j2))  !interpolating
        int2 = (z(k) - z_prev(k2)) / (z_prev(k2+1) - z_prev(k2))  !constants
        out_var(j,k) = (1.0_r2 - int1) * (1.0_r2 - int2) * in_var(j2,k2) + &
                       int1 * (1.0_r2 - int2) * in_var(j2+1,k2) + &
                       int1 * int2 * in_var(j2+1,k2+1) + &   !bilinear
                       (1.0_r2 - int1) * int2 * in_var(j2,k2+1) !interpolation
      END DO
    END DO

    RETURN
  END SUBROUTINE inter_var

  SUBROUTINE u_BCS(u, t)
    !Boundary conditions for total azimuthal velocity (including CCF)
    USE parameters
    IMPLICIT NONE

    REAL    (r2), INTENT(IN)  :: t
    REAL    (r2), INTENT(OUT) :: u(0:nx,0:nz)
    INTEGER (i1)              :: k

    u(0,:) = Re1 + Re1_mod * COS(om1 * t) + &
             eps1 * Re1 * (1.0_r2 / eta - 1.0_r2) * COS(freq1 * z(:))
    u(nx,:) = Re2 + Re2_mod * COS(om2 * t) - &
             eps2 * COS(freq2 * z(:))

    IF (ABS(tau - 1.0_r2) < EPSILON(tau)) THEN
      u(:,0) = 0.0_r2
      u(:,nz) = 0.0_r2
    END IF

    RETURN
  END SUBROUTINE u_BCS

  SUBROUTINE z_BCS(zn, pn, t)
    !Boundary conditions for azimuthal vorticity
    USE parameters
    IMPLICIT NONE

    REAL (r2), INTENT(IN)  :: t, pn(0:nx,0:nz)
    REAL (r2), INTENT(OUT) :: zn(0:nx,0:nz)

    zn(0,:) = -(8.0_r2 * pn(1,:) - pn(2,:)) / (2.0_r2 * s(0) * dx2)
    zn(nx,:) = -(8.0_r2 * pn(nx1,:) - pn(nx-2,:)) / (2.0_r2 * s(nx) * dx2)

    IF (ABS(tau - 1.0_r2) < EPSILON(tau)) THEN
      zn(:,0) = -(8.0_r2 * pn(:,1) - pn(:,2)) / &
                 (2.0_r2 * (s(:)) * dz2)
      zn(:,nz) = -(8.0_r2 * pn(:,nz1) - pn(:,nz-2)) / &
                  (2.0_r2 * (s(:)) * dz2)
    ELSE
      zn(:,0) = (-tau / (s(:) * (1.0_r2 - tau))) * &
                (0.5_r2 * (-pn(:,2) + 4.0_r2 * pn(:,1)) / delz)
      zn(:,nz) = (tau / (s(:) * (1.0_r2 - tau))) * &
                (0.5_r2 * (pn(:,nz-2) - 4.0_r2 * pn(:,nz1)) / delz)
    END IF

    RETURN
  END SUBROUTINE z_BCS

  SUBROUTINE p_BCS(p)
    !Boundary conditions for stream-function, psi
    USE parameters
    IMPLICIT NONE

    REAL (r2), INTENT(OUT) :: p(0:nx,0:nz)

    p(0,:) = 0.0_r2
    p(nx,:) = 0.0_r2

    p(:,0) = 0.0_r2
    p(:,nz) = 0.0_r2

    RETURN
  END SUBROUTINE p_BCS

  SUBROUTINE b_BCS(bn)
    !Boundary conditions for azimuthal magnetic field
    USE parameters
    IMPLICIT NONE

    REAL (r2), INTENT(OUT) :: bn(0:nx,0:nz)

    IF (ABS(tau - 0.0_r2) < EPSILON(tau)) THEN
      bn(:,0) = 0.0_r2
      bn(:,nz) = 0.0_r2
    END IF

    RETURN
  END SUBROUTINE b_BCS

  SUBROUTINE j_BCS(jn)
    !Boundary conditions for azimuthal current
    USE parameters
    IMPLICIT NONE

    REAL (r2), INTENT(OUT) :: jn(0:nx,0:nz)

    jn(0,:) = 0.0_r2
    jn(nx,:) = 0.0_r2

    IF (ABS(tau - 1.0_r2) < EPSILON(tau)) THEN
      jn(:,0) = 0.0_r2
      jn(:,nz) = 0.0_r2
    END IF

    RETURN
  END SUBROUTINE j_BCS

  FUNCTION f1(index)
    USE parameters
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN) :: index
    REAL    (r2)             :: f1

    f1 = eps1 * COS(freq1*z(index))

    RETURN
  END FUNCTION f1   

  FUNCTION f2(index)
    USE parameters
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN) :: index
    REAL    (r2)             :: f2

    f2 = eps2 * COS(freq2*z(index)-pi)

    RETURN
  END FUNCTION f2 

END MODULE ic_bc
