MODULE io
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

    CHARACTER(7)             :: itos
    INTEGER (i1), INTENT(IN) :: n
    INTEGER (i1)             :: i, n_, d(7)
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
    OPEN (51, FILE = 'time_tau.dat')
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
    CLOSE (51)

    RETURN
  END SUBROUTINE close_files

  SUBROUTINE save_growth(t, ur, ur_prev, uz, uz_prev, pn, v, zn, bn, jn, &
                         growth, growth_vz)
    !Save fields at particular points in (x,z)-plane
    USE parameters
    IMPLICIT NONE

    REAL (r2), INTENT(IN)  :: t, ur(0:nx,0:nz), uz(0:nx,0:nz), &
                              pn(0:nx,0:nz), bn(0:nx,0:nz), &
                              jn(0:nx,0:nz), v(0:nx,0:nz), &
                              zn(0:nx,0:nz), ur_prev(0:nx,0:nz), &
                              uz_prev(0:nx,0:nz)
    REAL (r2), INTENT(OUT) :: growth, growth_vz
    REAL (r2), SAVE        :: min_p, max_p, min_ur, max_ur, min_uz, max_uz
    INTEGER (i1)           :: zpos, xpos

    !growth rate of vortices
    growth = LOG(ABS(ur(nx/2,nz/2)/ur_prev(nx/2,nz/2))) / (dt * save_rate)
    growth_vz = LOG(ABS(uz(nx/2,nz/2)/uz_prev(nx/2,nz/2))) / (dt * save_rate)

    xpos = nx/2
    zpos = nz/2   !position at which to save fields

    WRITE(20, '(12e17.9)') t, ur(nx/2,nz/2), ur(nx/2,nz*2/gamma), &
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
    !Save torques
    USE parameters
    IMPLICIT NONE

    REAL (r2), INTENT(IN) :: t, v(0:nx,0:nz)
    INTEGER (i1)          :: k
    REAL (r2)             :: xi, C1, C2, G1(0:nz), G2(0:nz), G1_, G2_

    xi = Re1_mod * COS(om1 * t) - eta * Re2_mod * COS(om2 * t)
    xi = xi / (Re1 - eta * Re2)

    C1 = (-2.0_r2 * (1.0_r2 + xi)) / (eta * (1.0_r2 + eta))

    C2 = 1.0_r2 / (Re1 - eta * Re2)

    G1(:) = C1 + C2 * (0.5_r2 * (4.0_r2 * v(1,:) - v(2,:))) / delx
    G2(:) = C1 + (C2 / eta**2) * (0.5_r2 * (v(nx-2,:) - 4.0_r2 * &
            v(nx-1,:))) / delx

    G1_ = sum(G1(:))
    G2_ = sum(G2(:))
    !PRINT*,C1,C2
    WRITE (33, '(3e17.9)') t, G1_, G2_ !G1(nz/4), G2(nz/4)

    RETURN
  END SUBROUTINE save_torque

  SUBROUTINE save_xsect(ur, uz, pn, ut, zt, bt, jt, t, p)
    !Save cross-sections of fields for use in IDL
    USE parameters
    USE ic_bc, ONLY : x, z
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN) :: p
    REAL (r2),    INTENT(IN) :: ur(0:nx,0:nz), uz(0:nx,0:nz), &
                                pn(0:nx,0:nz), ut(0:nx,0:nz), &
                                zt(0:nx,0:nz), bt(0:nx,0:nz), &
                                jt(0:nx,0:nz), t
    INTEGER (i1)             :: j, k

    OPEN (32, STATUS = 'unknown', FILE = 'xsect'//itos(p)//'.dat')

    WRITE (32, '(e17.9)') t
    WRITE (32, '(2i5)') nx, nz
    WRITE (32, '(e17.9)') ((ur(j,k), j = 0, nx), k = 0, nz)
    WRITE (32, '(e17.9)') ((uz(j,k), j = 0, nx), k = 0, nz)
    WRITE (32, '(e17.9)') ((pn(j,k), j = 0, nx), k = 0, nz)
    WRITE (32, '(e17.9)') ((ut(j,k), j = 0, nx), k = 0, nz)
    WRITE (32, '(e17.9)') ((zt(j,k), j = 0, nx), k = 0, nz)
    WRITE (32, '(e17.9)') ((bt(j,k), j = 0, nx), k = 0, nz)
    WRITE (32, '(e17.9)') ((jt(j,k), j = 0, nx), k = 0, nz)
    WRITE (32, '(e17.9)') (x(j), j = 0, nx)
    WRITE (32, '(e17.9)') (z(k), k = 0, nz)

    CLOSE (32)

    RETURN
  END SUBROUTINE save_xsect

  SUBROUTINE save_3d(u_r, u_t, u_z, pn, p)
    !Save 3D isosurface for use in OpenDX
    USE parameters
    USE ic_bc, ONLY : x_, th, z
    IMPLICIT NONE

    INTEGER (i1), INTENT(IN) :: p
    REAL (r2),    INTENT(IN) :: u_r(0:nx,0:nz), u_t(0:nx,0:nz), &
                                u_z(0:nx,0:nz), pn(0:nx,0:nz)
    REAL (r2)                :: hel(0:nx,0:nz)
    INTEGER (i1)             :: j, k, l

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

    REAL (r2), INTENT(IN)  :: u_r(0:nx,0:nz), u_t(0:nx,0:nz), u_z(0:nx,0:nz)
    REAL (r2), INTENT(OUT) :: hel(0:nx, 0:nz)
    REAL (r2)              :: u_r_z(0:nx,0:nz), u_t_z(0:nx,0:nz), &
                              u_t_x(0:nx,0:nz), u_z_x(0:nx,0:nz)
    INTEGER (i1)           :: k

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
    REAL (r2),    INTENT(IN) :: t, pn(0:nx,0:nz), v(0:nx,0:nz), &
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
    REAL (r2),    INTENT(IN) :: t
    REAL (r2)                :: growth_rate, growth_rate_vz

    CALL vr_vz(psi%old, vr, vz)   !get radial, axial velocities
    IF (save_part) THEN
      CALL particle(vr, vrold, vz, vzold, x_pos, z_pos) !save particle
    END IF                                              !path
   !IF ((Re1 /= 0.0_r2) .OR. (Re2 /= 0.0_r2)) THEN
   !  CALL save_torque(t, unew)
   !END IF
    IF ((p /= p_start) .AND. ((p - p_start) > save_rate)) THEN
      CALL save_growth(t, vr, vrold, vz, vzold, psi%old, ut%new, zt%new, &
                       bt%old, jt%old, growth_rate, growth_rate_vz)
      IF ((om1 == 0.0_r2) .AND. (om2 == 0.0_r2)) THEN
        IF ((ABS(growth_rate_vz) < 1e-8_r2) .AND. &  !if vr saturated
            (ABS(vr(nx/2, nz/2)) > 1e-3_r2)) THEN
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
    REAL (r2),    INTENT(IN) :: t
    LOGICAL                  :: run_exist

    IF (mycol == 0) THEN
      INQUIRE(FILE='RUNNING', EXIST=run_exist)  !does 'RUNNING' exist?
      IF (.NOT. run_exist) THEN  !if not then finish
        PRINT*
        WRITE(6, '(A33, i2)'), 'Stop requested.  Ending process ', mycol
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
    REAL (r2),    INTENT(IN) :: t
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

    REAL (r2),    INTENT(IN) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                jn(0:nx,0:nz)
    INTEGER (i1), INTENT(IN) :: p
    INTEGER (i1)             :: j, k

    OPEN (50, FILE = 'end_state.dat')

    WRITE(50, *) nx
    WRITE(50, *) nz
    WRITE(50, '(i8)') p
    WRITE(50, '(e19.7)') dt
    WRITE(50, '(e19.7)') ((u(j,k), k = 0, nz), j = 0, nx)
    WRITE(50, '(e19.7)') ((zn(j,k), k = 0, nz), j = 0, nx)
    WRITE(50, '(e19.7)') ((pn(j,k), k = 0, nz), j = 0, nx)
    WRITE(50, '(e19.7)') ((bn(j,k), k = 0, nz), j = 0, nx)
    WRITE(50, '(e19.7)') ((jn(j,k), k = 0, nz), j = 0, nx)

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

  SUBROUTINE particle (vr, vrold, vz, vzold, xold, zold)
    !Calculate and save a particle path
    USE parameters
    IMPLICIT NONE

    REAL (r2), INTENT(IN)    :: vr(0:nx, 0:nz), vz(0:nx, 0:nz), &
                                vrold(0:nx, 0:nz), vzold(0:nx, 0:nz)
    REAL (r2), INTENT(INOUT) :: xold, zold
    INTEGER (i1)             :: xmin, xplu, zmin, zplu, j
    REAL (r2)                :: c1, c2, rvel, zvel, xnew, znew, del_t

    del_t = dt / 1.0_r2

    DO j = 1, 1
      xmin = INT(xold,i1)
      xplu = INT(xold + 1.0_r2,i1)
      zmin = INT(zold,i1)
      zplu = INT(zold + 1.0_r2,i1)

      c1 = (xold - xmin) / (xplu - xmin)
      c2 = (zold - zmin) / (zplu - zmin)

      rvel = (1.0_r2 - c1) * (1.0_r2 - c2) * vrold(xmin, zmin) + &
             c1 * (1.0_r2 - c2) * vrold(xplu, zmin) + &
             c1 * c2 * vrold(xplu, zplu) + &
             (1.0_r2 - c1) * c2 * vrold(xmin, zplu)

      zvel = (1.0_r2 - c1) * (1.0_r2 - c2) * vzold(xmin, zmin) + &
             c1 * (1.0_r2 - c2) * vzold(xplu, zmin) + &
             c1 * c2 * vzold(xplu, zplu) + &
             (1.0_r2 - c1) * c2 * vzold(xmin, zplu)

      xnew = xold + del_t * rvel
      znew = zold + del_t * zvel
      xold = xnew
      zold = znew
    END DO

    WRITE (24, '(2e17.9)') xnew / nx, znew * gamma / nz

    RETURN
  END SUBROUTINE particle

  SUBROUTINE get_timestep()
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
