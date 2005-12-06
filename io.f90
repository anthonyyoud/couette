MODULE io
  !Routines to do with input/output to/from files.
  USE parameters
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: open_files, close_files, save_xsect, save_3d, &
            save_surface, write_data, terminate, save_run, &
            end_state

  CONTAINS

  FUNCTION itos(n)
    !Function to convert an integer into a string of length 7
    USE parameters, ONLY : i1
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN) :: n
    INTEGER (i1)             :: i, n_, d(7)
    CHARACTER(7)             :: itos
    CHARACTER                :: c(0:9) = (/'0','1','2','3','4','5',&
                                           '6','7','8','9'/)

    n_ = n
    DO i = 7, 1, -1
      d(i) = MOD(n_,10)
      n_ = n_ / 10
    END DO

    itos = c(d(1))//c(d(2))//c(d(3))//c(d(4))//c(d(5))//c(d(6))//c(d(7))

    RETURN
  END FUNCTION itos

  SUBROUTINE open_files()
    !Open runtime files
    USE parameters
    IMPLICIT NONE

    OPEN (20, STATUS = 'unknown', FILE = 'u_growth.dat')
    OPEN (22, STATUS = 'unknown', FILE = 'max_psi.dat')
    OPEN (24, STATUS = 'unknown', FILE = 'particle.dat')
    OPEN (33, STATUS = 'unknown', FILE = 'torque.dat')
    IF (save3d) OPEN (35, STATUS = 'unknown', FILE = 'isosurface.dat')
    OPEN (36, STATUS = 'unknown', FILE = 'energy.dat')
    IF (auto_Re .OR. hyst_Re) THEN
      OPEN (37, STATUS = 'unknown', FILE = 'auto_Re.dat')
    END IF
    OPEN (51, FILE = 'time_tau.dat')
    IF (divergence) THEN
      OPEN (97, FILE = 'divergence.dat')
    END IF
    OPEN (99, FILE = 'RUNNING')
    CLOSE (99)

    RETURN
  END SUBROUTINE open_files

  SUBROUTINE close_files()
    !Close runtime files
    USE parameters
    IMPLICIT NONE

    CLOSE (20)
    CLOSE (22)
    CLOSE (24)
    CLOSE (33)
    IF (save3d) CLOSE (35)
    CLOSE (36)
    IF (auto_Re .OR. hyst_Re) THEN
      CLOSE (37)
    END IF
    CLOSE (51)
    IF (divergence) THEN
      CLOSE (97)
    END IF

    RETURN
  END SUBROUTINE close_files

  SUBROUTINE save_growth(t, ur, ur_prev, uz, uz_prev, pn, v, zn, bn, jn, &
                         growth, growth_vz)
    !Save fields at particular points in (x,z)-plane
    USE parameters
    IMPLICIT NONE

    REAL    (r2), INTENT(IN)  :: t, ur(0:nx,0:nz), uz(0:nx,0:nz), &
                                 pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                 jn(0:nx,0:nz), v(0:nx,0:nz), &
                                 zn(0:nx,0:nz), ur_prev(0:nx,0:nz), &
                                 uz_prev(0:nx,0:nz)
    REAL    (r2), INTENT(OUT) :: growth, growth_vz
    INTEGER (i1)              :: zpos, xpos
    REAL    (r2), SAVE        :: min_p, max_p, min_ur, max_ur, min_uz, max_uz

    !growth rate of vortices
    growth = LOG(ABS(ur(nx/2,nz/2)/ur_prev(nx/2,nz/2))) / (dt * save_rate)
    growth_vz = LOG(ABS(uz(nx/2,nz/2)/uz_prev(nx/2,nz/2))) / (dt * save_rate)

    xpos = nx/2
    zpos = nz/2   !position at which to save fields

    WRITE(20, '(12e17.9)') t, ur(nx/2,nz/2), ur(nx/2,0), &
                           growth, growth_vz, &
                           uz(xpos,zpos), &
                           pn(nx/4,3*nz/4), v(nx/2,nz/2), &
                           zn(nx/2,nz/4), bn(nx/2,nz/4), jn(nx/2,nz/2), &
                           Re1 + Re1_mod * COS(om1 * t)

    IF (MAXVAL(ur) > max_ur) THEN
      max_ur = MAXVAL(ur)
    END IF
    IF (MINVAL(ur) < min_ur) THEN
      min_ur = MINVAL(ur)
    END IF

    IF (MAXVAL(uz) > max_uz) THEN        !calculate maximum/minimum
      max_uz = MAXVAL(uz)               !values of fields
    END IF
    IF (MINVAL(uz) < min_uz) THEN
      min_uz = MINVAL(uz)
    END IF

    IF (MAXVAL(pn) > max_p) THEN
      max_p = MAXVAL(pn)
    END IF
    IF (MINVAL(pn) < min_p) THEN
      min_p = MINVAL(pn)
    END IF

    !WRITE (22, '(7e17.9)') t, max_p, min_p, max_ur, min_ur, max_uz, min_uz

    RETURN
  END SUBROUTINE save_growth

  SUBROUTINE save_torque(t, v)
    !Save torques on inner and outer cylinders as well as torque due to CCF
    USE parameters
    USE variables, ONLY : integrate_z
    USE ic_bc, ONLY : Re_1, Re_2
    IMPLICIT NONE

    REAL    (r2), INTENT(IN) :: t, v(0:nx,0:nz)
    REAL    (r2)             :: A1G1, A1G2, A2G1, A2G2, G1, G2, Gc, &
                                uc(0:nx), var1(0:nz), var2(0:nz), &
                                z_int1, z_int2

    CALL get_CCF(t, uc)
    
    A1G1 = (2.0_r2 * pi * eta**2) / one_eta**2
    A1G2 = (2.0_r2 * pi) / one_eta**2

    A2G1 = Re_1(t) * one_eta * gamma / eta
    A2G2 = Re_2(t) * one_eta * gamma

    var1(:) = (0.5_r2 * (-3.0_r2 * v(0,:) + 4.0_r2 * v(1,:) - v(2,:))) / delx
    var2(:) = (0.5_r2 * (v(nx-2,:) - 4.0_r2 * v(nx1,:) + &
               3.0_r2 * v(nx,:))) / delx

    CALL integrate_z(var1, z_int1)
    CALL integrate_z(var2, z_int2)

    G1 = A1G1 * (A2G1 - z_int1)
    G2 = A1G2 * (A2G2 - z_int2)
    Gc = 4.0_r2 * pi * eta * gamma * uc(0) / &
         ((1.0_r2 - eta**2) * one_eta)
    
    WRITE (33, '(5e17.9)') t, G1, G2, Gc, G1 / Gc

    RETURN
  END SUBROUTINE save_torque

  SUBROUTINE get_CCF(t, uc)
    !(Dimensionless) Circular Couette flow
    USE parameters
    USE ic_bc, ONLY : s, Re_1, Re_2
    IMPLICIT NONE

    REAL (r2), INTENT(IN)  :: t
    REAL (r2), INTENT(OUT) :: uc(0:nx)
    REAL (r2)              :: A, B
  
    A = (Re_2(t) - eta * Re_1(t)) / (1.0_r2 + eta)
    B = eta * (Re_1(t) - eta * Re_2(t)) / ((1.0_r2 + eta) * one_eta**2)

    uc = A * s / one_eta + B * one_eta / s

    RETURN
  END SUBROUTINE get_CCF
                                              
  SUBROUTINE save_xsect(ur, uz, pn, ut, zt, bt, jt, t, p)
    !Save cross-sections of fields for use in IDL
    USE parameters
    USE ic_bc, ONLY : x, z
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN) :: p
    REAL    (r2), INTENT(IN) :: ur(0:nx,0:nz), uz(0:nx,0:nz), &
                                pn(0:nx,0:nz), ut(0:nx,0:nz), &
                                zt(0:nx,0:nz), bt(0:nx,0:nz), &
                                jt(0:nx,0:nz), t
    INTEGER (i1)             :: j, k

    OPEN (32, STATUS = 'unknown', FILE = 'xsect'//itos(p)//'.dat', &
          FORM = 'unformatted')

    !WRITE (32, '(e17.9)') t
    !WRITE (32, '(2i5)') nx, nz
    !WRITE (32, '(e17.9)') ((ur(j,k), j = 0, nx), k = 0, nz)
    !WRITE (32, '(e17.9)') ((uz(j,k), j = 0, nx), k = 0, nz)
    !WRITE (32, '(e17.9)') ((pn(j,k), j = 0, nx), k = 0, nz)
    !WRITE (32, '(e17.9)') ((ut(j,k), j = 0, nx), k = 0, nz)
    !WRITE (32, '(e17.9)') ((zt(j,k), j = 0, nx), k = 0, nz)
    !WRITE (32, '(e17.9)') ((bt(j,k), j = 0, nx), k = 0, nz)
    !WRITE (32, '(e17.9)') ((jt(j,k), j = 0, nx), k = 0, nz)
    !WRITE (32, '(e17.9)') (x(j), j = 0, nx)
    !WRITE (32, '(e17.9)') (z(k), k = 0, nz)
    
    WRITE (32) t
    WRITE (32) nx, nz
    WRITE (32) ur
    WRITE (32) uz
    WRITE (32) pn
    WRITE (32) ut
    WRITE (32) zt
    WRITE (32) bt
    WRITE (32) jt
    WRITE (32) x
    WRITE (32) z

    CLOSE (32)

    RETURN
  END SUBROUTINE save_xsect

  SUBROUTINE save_3d(u_r, u_t, u_z, pn, p)
    !Save 3D isosurface for use in OpenDX
    USE parameters
    USE ic_bc, ONLY : x_, th, z
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN) :: p
    REAL    (r2), INTENT(IN) :: u_r(0:nx,0:nz), u_t(0:nx,0:nz), &
                                u_z(0:nx,0:nz), pn(0:nx,0:nz)
    INTEGER (i1)             :: j, k, l
    REAL    (r2)             :: hel(0:nx,0:nz)

    !OPEN (35, STATUS = 'unknown', FILE = 'p3d'//itos(p)//'.dat')

    IF (iso_hel) CALL helicity(u_r, u_t, u_z, hel)

    WRITE(35,*) 'nx, nt, nz = ', nx, nt, nz, p

    IF (iso_hel) THEN            !save helicity u.(curl u)
      DO j = 0, nx
        DO l = 0, nt
          DO k = 0, nz
            WRITE(35,'(4e11.3)') x_(j) * COS(th(l)), x_(j) * SIN(th(l)), &
                                 z(k), hel(j,k)
          END DO
        END DO
      END DO
    ELSE
      DO j = 0, nx
        DO l = 0, nt
          DO k = 0, nz
            WRITE(35,'(4e11.3)') x_(j) * COS(th(l)), x_(j) * SIN(th(l)), &
                                 z(k), pn(j,k)
          END DO
        END DO
      END DO
    END IF

    WRITE(35,*)

    !CLOSE(35)

    RETURN
  END SUBROUTINE save_3d

  SUBROUTINE helicity(u_r, u_t, u_z, hel)
    !Calculate helicity u.(curl u)
    USE parameters
    USE ic_bc, ONLY : s
    USE derivs, ONLY : deriv_x, deriv_z
    IMPLICIT NONE

    REAL    (r2), INTENT(IN)  :: u_r(0:nx,0:nz), u_t(0:nx,0:nz), u_z(0:nx,0:nz)
    REAL    (r2), INTENT(OUT) :: hel(0:nx, 0:nz)
    INTEGER (i1)              :: k
    REAL    (r2)              :: u_r_z(0:nx,0:nz), u_t_z(0:nx,0:nz), &
                                 u_t_x(0:nx,0:nz), u_z_x(0:nx,0:nz)

    CALL deriv_z(u_r, u_r_z)
    CALL deriv_z(u_t, u_t_z)
    CALL deriv_x(u_t, u_t_x)
    CALL deriv_x(u_z, u_z_x)

    DO k = 0, nz
      hel(:,k) = 0.5_r2 * delx * (-u_r(:,k) * u_t_z(:,k) + &
                                   u_t(:,k) * u_r_z(:,k) - &
                                   u_t(:,k) * u_z_x(:,k) + &
                                   u_z(:,k) * u_t_x(:,k)) + &
                 (s(:) * u_t(:,k) * u_z(:,k)) / one_eta
    END DO

    RETURN
  END SUBROUTINE helicity

  SUBROUTINE save_surface(pn, v, zn, ur, uz, bn, jn, p, t)
    !Save surfaces of fields for use in gnuplot
    USE parameters
    USE ic_bc, ONLY : x, z
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN) :: p
    REAL    (r2), INTENT(IN) :: t, pn(0:nx,0:nz), v(0:nx,0:nz), &
                                zn(0:nx,0:nz), bn(0:nx,0:nz), &
                                jn(0:nx,0:nz), ur(0:nx,0:nz), &
                                uz(0:nx,0:nz)
    INTEGER (i1)             :: j, k

    OPEN (19, STATUS = 'unknown', FILE = 'p'//itos(p)//'.dat')
    !OPEN (70, STATUS = 'unknown', FILE = 'sin.dat')
    OPEN (21, STATUS = 'unknown', FILE = 'z'//itos(p)//'.dat')
    OPEN (23, STATUS = 'unknown', FILE = 'u'//itos(p)//'.dat')
    OPEN (30, STATUS = 'unknown', FILE = 'vr'//itos(p)//'.dat')
    OPEN (31, STATUS = 'unknown', FILE = 'vz'//itos(p)//'.dat')
    OPEN (34, STATUS = 'unknown', FILE = 'b'//itos(p)//'.dat')
    OPEN (35, STATUS = 'unknown', FILE = 'j'//itos(p)//'.dat')

    WRITE(30, '(2A, i10, e19.7)') '#', 'p=', p, t
    WRITE(31, '(2A, i10, e19.7)') '#', 'p=', p, t
    WRITE(34, '(2A, i10, e19.7)') '#', 'p=', p, t
    WRITE(35, '(2A, i10, e19.7)') '#', 'p=', p, t
    WRITE(23, '(2A, i10, e19.7)') '#', 'p=', p, t
    WRITE(19, '(2A, i10, e19.7)') '#', 'p=', p, t
    WRITE(21, '(2A, i10, e19.7)') '#', 'p=', p, t
    !WRITE(70, '(2A, i10, e19.7)') '#', 'p=', p, t

    DO j = 0, nx
      WRITE(30, '(3e19.7)') (x(j), z(k), ur(j,k), k = 0, nz)
      WRITE(30, *)
      WRITE(31, '(3e19.7)') (x(j), z(k), uz(j,k), k = 0, nz)
      WRITE(31, *)
      WRITE(34, '(3e19.7)') (x(j), z(k), bn(j,k), k = 0, nz)
      WRITE(34, *)
      WRITE(35, '(3e19.7)') (x(j), z(k), jn(j,k), k = 0, nz)
      WRITE(35, *)
      WRITE(23, '(3e19.7)') (x(j), z(k), v(j,k), k = 0, nz)
      WRITE(23, *)
      WRITE(19, '(3e19.7)') (x(j), z(k), pn(j,k), k = 0, nz)
      WRITE(19, *)
      WRITE(21, '(3e19.7)') (x(j), z(k), zn(j,k), k = 0, nz)
      WRITE(21, *)
     !WRITE(70, '(3e19.7)') (x(j), z(k), SIN(pi*x(j))*&
     !                      SIN(alpha*z(k)), k = 0, nz)
     !WRITE(70, *)
    END DO

    CLOSE (19)
    CLOSE (21)
    CLOSE (23)
    CLOSE (30)
    CLOSE (31)
    CLOSE (34)
    CLOSE (35)
    !CLOSE (70)

    RETURN
  END SUBROUTINE save_surface

  SUBROUTINE write_data(p, p_start, t)
    !Calculations of various quantities; save fields
    USE parameters
    USE variables
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN) :: p, p_start
    REAL    (r2), INTENT(IN) :: t
    REAL    (r2)             :: growth_rate, growth_rate_vz, growth_diff, &
                                growth_error
    REAL    (r2), SAVE       :: growth_rate_vz_old

    CALL vr_vz(psi%old, vr, vz)   !get radial, axial velocities
    IF (save_part) THEN
      CALL particle(vr, vrold, vz, vzold) !save particle
    END IF                                              !path
    CALL save_torque(t, ut%new)
    IF ((p /= p_start) .AND. ((p - p_start) > save_rate)) THEN
      CALL save_growth(t, vr, vrold, vz, vzold, psi%old, ut%new, zt%new, &
                       bt%old, jt%old, growth_rate, growth_rate_vz)
      growth_diff = ABS(growth_rate_vz - growth_rate_vz_old)
      growth_error = growth_diff / ABS(growth_rate_vz_old)
      !print*, Re1, growth_diff, growth_error
      IF (auto_Re .AND. ((growth_diff < growth_tol) .OR. &
          (dec_Re .AND. zero_Re .AND. (vz(nx/2,nz/2) < 1E-6_r2)))) THEN
        CALL increment_Re(growth_rate_vz, t)
      END IF
      !WRITE (37, '(4e17.9)') t, growth_rate_vz_old, growth_rate_vz, growth_diff
      growth_rate_vz_old = growth_rate_vz
      CALL save_energy(vr, ut%new, vz, t)
      IF ((ABS(om1 - 0.0_r2) < EPSILON(om1)) .AND. &
          (ABS(om2 - 0.0_r2) < EPSILON(om2))) THEN
        IF (hyst_Re) THEN
          IF (((ABS(growth_rate_vz) < 1E-6_r2) .AND. &
               (ABS(vz(nx/2,nz/2)) > 1E-3_r2))) THEN
            CALL increment_hysteresis(.TRUE., t)
          ELSE IF ((ABS(growth_rate_vz) > 1E-6_r2) .AND. &
                       (growth_error < 1E-8_r2)) THEN
            CALL increment_hysteresis(.FALSE., t)
          END IF
        END IF
        IF ((ABS(growth_rate_vz) < 1E-8_r2) .AND. &  !if vr saturated
            (ABS(vr(nx/2, nz/2)) > 1E-3_r2)) THEN
          saturated = saturated + 1
          IF (saturated > 4) THEN
            IF ((.NOT. auto_tau) .OR. (tau == tau_end)) THEN  !if tau not auto
              CALL save_time_tau(tau, t)                     !or tau at end
              CALL end_state(ut%old, zt%old, psi%old, bt%old, jt%old, p) !finish
            ELSE IF (tau < 1.0_r2) THEN
              CALL save_time_tau(tau, t)
              tau = tau + tau_step   !increment tau
              PRINT*, 'tau = ', tau
              CALL save_xsect(vr, vz, psi%old, ut%new, zt%new, &
                              bt%old, jt%old, t, p)
              CALL save_surface(psi%old, ut%new, zt%new, vr, vz, &
                                bt%old, jt%old, p, t)
            END IF
          END IF
        ELSE
          saturated = 0
        END IF
      END IF
    END IF

    RETURN
  END SUBROUTINE write_data

  SUBROUTINE terminate(p, t)
    !Terminate run if file 'RUNNING' does not exist in run directory
    USE parameters
    USE variables
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN) :: p
    REAL    (r2), INTENT(IN) :: t
    LOGICAL                  :: run_exist

    IF (mycol == 0) THEN
      INQUIRE(FILE='RUNNING', EXIST=run_exist)  !does 'RUNNING' exist?
      IF (.NOT. run_exist) THEN  !if not then finish
        PRINT*
        WRITE(6, '(A33, i2)') 'Stop requested.  Ending process ', mycol
        PRINT*, 'Saving end state...'
        CALL end_state(ut%old, zt%old, psi%old, bt%old, jt%old, p)
        CALL save_xsect(vr, vz, psi%old, ut%new, zt%new, &
                        bt%old, jt%old, t, p)
        CALL save_surface(psi%old, ut%old, zt%old, &
                          vr, vz, bt%old, jt%old, p, t)
        end_proc = 1  !flag to send to all processes
      END IF
    END IF

    CALL SLTIMER(3)                                         !tell all other
    CALL IGEBR2D(ictxt, 'A', ' ', 1, 1, end_proc, 1, 0, 0)  !processes to finish
    CALL SLTIMER(3)                                         !if end_proc=1

    RETURN
  END SUBROUTINE terminate

  SUBROUTINE save_run(p, t)
    !Save cross-sections, surfaces if file 'SAVE' exists in run directory
    USE parameters
    USE variables
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN) :: p
    REAL    (r2), INTENT(IN) :: t
    LOGICAL                  :: save_exist

    INQUIRE(FILE='SAVE', EXIST=save_exist)   !does 'SAVE' exist?
    IF (save_exist) THEN   !if so then save fields
      CALL save_xsect(vr, vz, psi%old, ut%new, zt%new, &
                      bt%old, jt%old, t, p)
      CALL save_surface(psi%old, ut%old, zt%old, &
                        vr, vz, bt%old, jt%old, p, t)
      OPEN (98, FILE = 'SAVE')
      CLOSE (98, STATUS = 'delete')
    END IF

    RETURN
  END SUBROUTINE save_run

  SUBROUTINE end_state(u, zn, pn, bn, jn, p)
    !Save variables for use in a restarted run
    USE parameters
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN) :: p
    REAL    (r2), INTENT(IN) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                jn(0:nx,0:nz)
    INTEGER (i1)             :: j, k

    OPEN (50, FILE = 'end_state.dat', FORM='unformatted')

    WRITE(50) nx
    WRITE(50) nz
    WRITE(50) p
    WRITE(50) dt
    WRITE(50) u
    WRITE(50) zn
    WRITE(50) pn
    WRITE(50) bn
    WRITE(50) jn

    CLOSE (50)

    OPEN (99, FILE = 'RUNNING')  !delete 'RUNNING' to finish run
    CLOSE (99, STATUS = 'delete')

    RETURN
  END SUBROUTINE end_state

  SUBROUTINE save_time_tau (tau, t)
    !Save the times at which tau is incremented
    USE parameters, ONLY : r2
    IMPLICIT NONE

    REAL (r2), INTENT(IN) :: tau, t

    WRITE(51, '(2e17.9)') t, tau

    RETURN
  END SUBROUTINE save_time_tau

  SUBROUTINE particle (vr, vrold, vz, vzold)
    !Calculate and save a particle path
    USE parameters
    USE variables, ONLY : xold, zold
    IMPLICIT NONE

    REAL    (r2), INTENT(IN)    :: vr(0:nx, 0:nz), vz(0:nx, 0:nz), &
                                   vrold(0:nx, 0:nz), vzold(0:nx, 0:nz)
    INTEGER (i1)                :: xmin(num_pars), xplu(num_pars), &
                                   zmin(num_pars), zplu(num_pars), i, j
    REAL    (r2)                :: cc1(num_pars), cc2(num_pars), &
                                   rvel(num_pars), zvel(num_pars), &
                                   xnew(num_pars), znew(num_pars), del_t

    del_t = dt / 1.0_r2

    DO j = 1, 1
      xmin = INT(xold,i1)
      xplu = INT(xold + 1.0_r2,i1)
      zmin = INT(zold,i1)
      zplu = INT(zold + 1.0_r2,i1)

      cc1 = (xold - xmin) / (xplu - xmin)
      cc2 = (zold - zmin) / (zplu - zmin)

      DO i = 1, num_pars
        rvel(i) = (1.0_r2 - cc1(i)) * (1.0_r2 - cc2(i)) * &
                  vrold(xmin(i), zmin(i)) + &
                  cc1(i)* (1.0_r2 - cc2(i)) * vrold(xplu(i), zmin(i)) + &
                  cc1(i) * cc2(i) * vrold(xplu(i), zplu(i)) + &
                  (1.0_r2 - cc1(i)) * cc2(i) * vrold(xmin(i), zplu(i))
      END DO

      DO i = 1, num_pars
        zvel(i) = (1.0_r2 - cc1(i)) * (1.0_r2 - cc2(i)) * &
                  vzold(xmin(i), zmin(i)) + &
                  cc1(i) * (1.0_r2 - cc2(i)) * vzold(xplu(i), zmin(i)) + &
                  cc1(i) * cc2(i) * vzold(xplu(i), zplu(i)) + &
                  (1.0_r2 - cc1(i)) * cc2(i) * vzold(xmin(i), zplu(i))
      END DO

      xnew = xold + del_t * rvel
      znew = zold + del_t * zvel
      xold = xnew
      zold = znew
    END DO

    WRITE (24, '(10e17.9)') xnew(1) / nx, znew(1) * gamma / nz, &
                            xnew(2) / nx, znew(2) * gamma / nz, &
                            xnew(3) / nx, znew(3) * gamma / nz, &
                            xnew(4) / nx, znew(4) * gamma / nz, &
                            xnew(5) / nx, znew(5) * gamma / nz

    RETURN
  END SUBROUTINE particle

  SUBROUTINE energy_CCF(t, Eccf)
    !Explicitly integrated (dimensionless) form of the kinetic energy in CCF
    USE parameters
    USE ic_bc, ONLY : Re_1, Re_2
    IMPLICIT NONE

    REAL (r2), INTENT(IN)  :: t
    REAL (r2), INTENT(OUT) :: Eccf
    REAL (r2)              :: a, b, c

    a = 0.25_r2 * (Re_2(t) - eta * Re_1(t))**2 * (1.0_r2 - eta**4)
    b = (Re_2(t) - eta * Re_1(t)) * (Re_1(t) - eta * Re_2(t)) * &
        eta * (1.0_r2 - eta**2)
    c = eta**2 * (Re_1(t) - eta * Re_2(t))**2 * LOG(1.0_r2 / eta)

    Eccf = gamma * pi * (a + b + c) / &
          ((1.0_r2 + eta)**2 * (one_eta)**4)
    
    RETURN
  END SUBROUTINE energy_CCF

  SUBROUTINE save_energy(ur, ut, uz, t)
    !Total kinetic energy in flow (including CCF)
    USE parameters
    USE variables, ONLY : integrate_r, integrate_z
    IMPLICIT NONE

    REAL (r2), INTENT(IN)  :: ur(0:nx,0:nz), ut(0:nx,0:nz), uz(0:nx,0:nz), t
    REAL (r2)              :: u2(0:nx,0:nz), int_r(0:nz), int_z, E, Eccf

    u2 = ur**2 + ut**2 + uz**2

    CALL energy_CCF(t, Eccf)
    CALL integrate_r(u2, int_r)
    CALL integrate_z(int_r, int_z)

    E = pi * int_z

    WRITE (36, '(3e17.9)') t, Eccf, E

    RETURN
  END SUBROUTINE save_energy

  SUBROUTINE increment_Re(growth_rate, t)
    !Automatically increment/decrement the Reynolds number to find critical
    !values based on the sign of the growth rate
    USE parameters
    IMPLICIT NONE

    REAL (r2), INTENT(IN) :: growth_rate, t
    REAL (r2), SAVE       :: Re_minus, Re_plus, growth_minus, growth_plus

    SELECT CASE (dec_RE)
      CASE (.FALSE.)
        IF (init_Re) THEN
          IF (growth_rate < 0.0_r2) THEN
            growth_minus = growth_rate
            Re_minus = Re1
            Re1 = Re1 + Re_incr
            gm_set = .TRUE.
          ELSE IF (growth_rate > 0.0_r2) THEN
            growth_plus = growth_rate
            Re_plus = Re1
            Re1 = Re1 - Re_incr
            gp_set = .TRUE.
          END IF
          init_Re = .FALSE.
        ELSE IF (.NOT. init_Re) THEN
          IF ((.NOT. gm_set) .OR. (.NOT. gp_set)) THEN
            IF (.NOT. gp_set) THEN
              IF (growth_rate < 0.0_r2) THEN
                growth_minus = growth_rate
                Re_minus = Re1
                Re1 = Re1 + Re_incr
              ELSE IF (growth_rate > 0.0_r2) THEN
                growth_plus = growth_rate
                Re_plus = Re1
                Re1 = 0.5_r2 * (Re1 + Re_minus)
                gp_set = .TRUE.
              END IF
            ELSE IF (.NOT. gm_set) THEN
              IF (growth_rate > 0.0_r2) THEN
                growth_plus = growth_rate
                Re_plus = Re1
                Re1 = Re1 - Re_incr
              ELSE IF (growth_rate < 0.0_r2) THEN
                growth_minus = growth_rate
                Re_minus = Re1
                Re1 = 0.5_r2 * (Re1 + Re_plus)
                gm_set = .TRUE.
              END IF
            END IF
          ELSE
            IF (growth_rate < 0.0_r2) THEN
              growth_minus = growth_rate
              Re_minus = Re1
              Re1 = 0.5_r2 * (Re1 + Re_plus)
            ELSE IF (growth_rate > 0.0_r2) THEN
              growth_plus = growth_rate
              Re_plus = Re1
              Re1 = 0.5_r2 * (Re1 + Re_minus)
            END IF
          END IF
        END IF
      CASE (.TRUE.)
        IF (init_Re) THEN
          IF (growth_rate < 0.0_r2) THEN
            growth_minus = growth_rate
            Re_minus = Re1
            Re1 = Re1 - Re_incr
            gm_set = .TRUE.
          ELSE IF (growth_rate > 0.0_r2) THEN
            STOP 'Decreasing Re boundary only but growth rate > 0'
          END IF
          init_Re = .FALSE.
        ELSE IF (.NOT. init_Re) THEN
          IF (.NOT. gp_set) THEN
            IF (growth_rate < 0.0_r2) THEN
              growth_minus = growth_rate
              Re_minus = Re1
              Re1 = Re1 - Re_incr
            ELSE IF (growth_rate > 0.0_r2) THEN
              growth_plus = growth_rate
              Re_plus = Re1
              Re1 = 0.0_r2
              zero_Re = .TRUE.
              gp_set = .TRUE.
            END IF
          ELSE
            IF (zero_Re) THEN
              Re1 = 0.5_r2 * (Re_minus + Re_plus)
              zero_Re = .FALSE.
            ELSE IF (.NOT. zero_Re) THEN
              IF (growth_rate < 0.0_r2) THEN
                growth_minus = growth_rate
                Re_minus = Re1
                Re1 = 0.5_r2 * (Re1 + Re_plus)
              ELSE IF (growth_rate > 0.0_r2) THEN
                growth_plus = growth_rate
                Re_plus = Re1
                Re1 = 0.0_r2
                zero_Re = .TRUE.
              END IF
            END IF
          END IF
        END IF
    END SELECT
    
    WRITE (37, '(6e17.9)') t, Re1, Re_minus, Re_plus, growth_minus, growth_plus
    
    IF (ABS(Re1 - Re_minus) < 1E-6_r2) THEN
      OPEN (99, FILE = 'RUNNING')
      CLOSE (99, STATUS = 'delete')
    END IF
    
    RETURN
  END SUBROUTINE
  
  SUBROUTINE increment_hysteresis(test, t)
    !Automatically increment/decrement the Reynolds number to find critical
    !values in a hysteresis region
    USE parameters
    USE variables
    IMPLICIT NONE

    REAL     (r2), INTENT(IN) :: t
    LOGICAL,       INTENT(IN) :: test
    REAL     (r2), SAVE       :: Re_minus, Re_plus

    IF (init_Re) THEN
      IF (test) THEN
        Re_plus = Re1
        Re1 = Re1 - Re_incr
      END IF
      init_Re = .FALSE.
    ELSE IF (.NOT. init_Re) THEN
      IF (test) THEN
        IF (.NOT. gm_set) THEN
          Re_plus = Re1
          Re1 = Re1 - Re_incr
        ELSE IF (gm_set) THEN
          Re_plus = Re1
          Re1 = 0.5_r2 * (Re_minus + Re_plus)
        END IF
      ELSE IF (.NOT. test) THEN
        Re_minus = Re1
        Re1 = 0.5_r2 * (Re_minus + Re_plus)
        gm_set = .TRUE.
        OPEN (50, FILE = 'end_state.dat', FORM = 'unformatted')
      
        READ (50)
        READ (50)
        READ (50)
        READ (50)
        READ (50) ut%new
        READ (50) zt%new
        READ (50) psi%new
        READ (50) bt%new
        READ (50) jt%new

        CLOSE (50)
      END IF
    END IF

    WRITE (37, '(4e17.9)') t, Re1, Re_minus, Re_plus
    
    IF (ABS(Re1 - Re_minus) < 1E-6_r2) THEN
      OPEN (99, FILE = 'RUNNING')
      CLOSE (99, STATUS = 'delete')
    END IF

    RETURN
  END SUBROUTINE increment_hysteresis
  
  SUBROUTINE get_timestep()
    !Obsolete: from attempt at adaptive timestep
    USE parameters
    IMPLICIT NONE

    REAL (r2) :: timestep

    IF (Re1 * Re2 >= 0.0_r2) THEN
      timestep = MAX(1.0_r2, ABS(Re1) + ABS(Re1_mod), &
                 ABS(Re2) + ABS(Re2_mod))
    ELSE
      timestep = MAX(1.0_r2, ABS(Re1) + ABS(Re1_mod) + &
                 ABS(Re2) + ABS(Re2_mod))
    END IF

    timestep = 1e-4_r2 * eta * 2.0_r2 * pi / (one_eta * timestep)

    PRINT*, timestep

    RETURN
  END SUBROUTINE get_timestep

END MODULE io
