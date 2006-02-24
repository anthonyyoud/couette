module io
  !Routines to do with input/output to/from files.
  use parameters
  implicit none

  private
  public :: open_files, close_files, save_xsect, save_3d, &
            save_surface, write_data, terminate, save_run, &
            end_state

  contains

  FUNCTION itos(n)
    !Function to convert an integer into a string of length 7
    use parameters, only : i1
    implicit none

    integer (i1), intent(in) :: n
    integer (i1)             :: i, n_, d(7)
    character(7)             :: itos
    character                :: c(0:9) = (/'0','1','2','3','4','5',&
                                           '6','7','8','9'/)

    n_ = n
    do i = 7, 1, -1
      d(i) = MOD(n_,10)
      n_ = n_ / 10
    end do

    itos = c(d(1))//c(d(2))//c(d(3))//c(d(4))//c(d(5))//c(d(6))//c(d(7))

    return
  end function itos

  subroutine open_files()
    !Open runtime files
    use parameters
    implicit none

    open (20, status = 'unknown', file = 'u_growth.dat')
    open (22, status = 'unknown', file = 'max_psi.dat')
    open (24, status = 'unknown', file = 'particle.dat')
    open (33, status = 'unknown', file = 'torque.dat')
    if (save3d) open (35, status = 'unknown', file = 'isosurface.dat')
    open (36, status = 'unknown', file = 'energy.dat')
    if (auto_Re .or. hyst_Re) then
      open (37, status = 'unknown', file = 'auto_Re.dat')
    end if
    open (51, file = 'time_tau.dat')
    if (divergence) then
      open (97, file = 'divergence.dat')
    end if
    open (99, file = 'RUNNING')
    close (99)

    return
  end subroutine open_files

  subroutine close_files()
    !Close runtime files
    use parameters
    implicit none

    close (20)
    close (22)
    close (24)
    close (33)
    if (save3d) close (35)
    close (36)
    if (auto_Re .or. hyst_Re) then
      close (37)
    end if
    close (51)
    if (divergence) then
      close (97)
    end if

    return
  end subroutine close_files

  subroutine save_growth(t, ur, ur_prev, uz, uz_prev, pn, v, zn, bn, jn, &
                         growth, growth_vz)
    !Save fields at particular points in (x,z)-plane
    use parameters
    implicit none

    real    (r2), intent(in)  :: t, ur(0:nx,0:nz), uz(0:nx,0:nz), &
                                 pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                 jn(0:nx,0:nz), v(0:nx,0:nz), &
                                 zn(0:nx,0:nz), ur_prev(0:nx,0:nz), &
                                 uz_prev(0:nx,0:nz)
    real    (r2), intent(out) :: growth, growth_vz
    integer (i1)              :: zpos, xpos
    real    (r2), save        :: min_p, max_p, min_ur, max_ur, min_uz, max_uz

    !growth rate of vortices
    growth = log(abs(ur(nx/2,nz/2)/ur_prev(nx/2,nz/2))) / (dt * save_rate)
    growth_vz = log(abs(uz(nx/2,nz/2)/uz_prev(nx/2,nz/2))) / (dt * save_rate)

    xpos = nx/2
    zpos = nz/2   !position at which to save fields

    write(20, '(12e17.9)') t, ur(nx/2,nz/2), ur(nx/2,0), &
                           growth, growth_vz, &
                           uz(xpos,zpos), &
                           pn(nx/4,3*nz/4), v(nx/2,nz/2), &
                           zn(nx/2,nz/4), bn(nx/2,nz/4), jn(nx/2,nz/2), &
                           Re1 + Re1_mod * cos(om1 * t)

    if (maxval(ur) > max_ur) then
      max_ur = maxval(ur)
    end if
    if (minval(ur) < min_ur) then
      min_ur = minval(ur)
    end if

    if (maxval(uz) > max_uz) then        !calculate maximum/minimum
      max_uz = maxval(uz)               !values of fields
    end if
    if (minval(uz) < min_uz) then
      min_uz = minval(uz)
    end if

    if (maxval(pn) > max_p) then
      max_p = maxval(pn)
    end if
    if (minval(pn) < min_p) then
      min_p = minval(pn)
    end if

    !write (22, '(7e17.9)') t, max_p, min_p, max_ur, min_ur, max_uz, min_uz

    return
  end subroutine save_growth

  subroutine save_torque(t, v)
    !Save torques on inner and outer cylinders as well as torque due to CCF
    use parameters
    use variables, only : integrate_z
    use ic_bc, only : Re_1, Re_2
    implicit none

    real    (r2), intent(in) :: t, v(0:nx,0:nz)
    real    (r2)             :: A1G1, A1G2, A2G1, A2G2, G1, G2, Gc, &
                                uc(0:nx), var1(0:nz), var2(0:nz), &
                                z_int1, z_int2

    call get_CCF(t, uc)
    
    A1G1 = (2.0_r2 * pi * eta**2) / one_eta**2
    A1G2 = (2.0_r2 * pi) / one_eta**2

    A2G1 = Re_1(t) * one_eta * gamma / eta
    A2G2 = Re_2(t) * one_eta * gamma

    var1(:) = (0.5_r2 * (-3.0_r2 * v(0,:) + 4.0_r2 * v(1,:) - v(2,:))) / delx
    var2(:) = (0.5_r2 * (v(nx-2,:) - 4.0_r2 * v(nx1,:) + &
               3.0_r2 * v(nx,:))) / delx

    call integrate_z(var1, z_int1)
    call integrate_z(var2, z_int2)

    G1 = A1G1 * (A2G1 - z_int1)
    G2 = A1G2 * (A2G2 - z_int2)
    Gc = 4.0_r2 * pi * eta * gamma * uc(0) / &
         ((1.0_r2 - eta**2) * one_eta)
    
    write (33, '(5e17.9)') t, G1, G2, Gc, G1 / Gc

    return
  end subroutine save_torque

  subroutine get_CCF(t, uc)
    !(Dimensionless) Circular Couette flow
    use parameters
    use ic_bc, only : s, Re_1, Re_2
    implicit none

    real (r2), intent(in)  :: t
    real (r2), intent(out) :: uc(0:nx)
    real (r2)              :: A, B
  
    A = (Re_2(t) - eta * Re_1(t)) / (1.0_r2 + eta)
    B = eta * (Re_1(t) - eta * Re_2(t)) / ((1.0_r2 + eta) * one_eta**2)

    uc = A * s / one_eta + B * one_eta / s

    return
  end subroutine get_CCF
                                              
  subroutine save_xsect(ur, uz, pn, ut, zt, bt, jt, t, p)
    !Save cross-sections of fields for use in IDL
    use parameters
    use ic_bc, only : x, z
    implicit none

    integer (i1), intent(in) :: p
    real    (r2), intent(in) :: ur(0:nx,0:nz), uz(0:nx,0:nz), &
                                pn(0:nx,0:nz), ut(0:nx,0:nz), &
                                zt(0:nx,0:nz), bt(0:nx,0:nz), &
                                jt(0:nx,0:nz), t
    integer (i1)             :: j, k

    open (32, status = 'unknown', file = 'xsect'//itos(p)//'.dat', &
          form = 'unformatted')

    !write (32, '(e17.9)') t
    !write (32, '(2i5)') nx, nz
    !write (32, '(e17.9)') ((ur(j,k), j = 0, nx), k = 0, nz)
    !write (32, '(e17.9)') ((uz(j,k), j = 0, nx), k = 0, nz)
    !write (32, '(e17.9)') ((pn(j,k), j = 0, nx), k = 0, nz)
    !write (32, '(e17.9)') ((ut(j,k), j = 0, nx), k = 0, nz)
    !write (32, '(e17.9)') ((zt(j,k), j = 0, nx), k = 0, nz)
    !write (32, '(e17.9)') ((bt(j,k), j = 0, nx), k = 0, nz)
    !write (32, '(e17.9)') ((jt(j,k), j = 0, nx), k = 0, nz)
    !write (32, '(e17.9)') (x(j), j = 0, nx)
    !write (32, '(e17.9)') (z(k), k = 0, nz)
    
    write (32) t
    write (32) nx, nz
    write (32) ur
    write (32) uz
    write (32) pn
    write (32) ut
    write (32) zt
    write (32) bt
    write (32) jt
    write (32) x
    write (32) z

    close (32)

    return
  end subroutine save_xsect

  subroutine save_3d(u_r, u_t, u_z, pn, p)
    !Save 3D isosurface for use in OpenDX
    use parameters
    use ic_bc, only : x_, th, z
    implicit none

    integer (i1), intent(in) :: p
    real    (r2), intent(in) :: u_r(0:nx,0:nz), u_t(0:nx,0:nz), &
                                u_z(0:nx,0:nz), pn(0:nx,0:nz)
    integer (i1)             :: j, k, l
    real    (r2)             :: hel(0:nx,0:nz)

    !open (35, status = 'unknown', file = 'p3d'//itos(p)//'.dat')

    if (iso_hel) call helicity(u_r, u_t, u_z, hel)

    write(35,*) 'nx, nt, nz = ', nx, nt, nz, p

    if (iso_hel) then            !save helicity u.(curl u)
      do j = 0, nx
        do l = 0, nt
          do k = 0, nz
            write(35,'(4e11.3)') x_(j) * cos(th(l)), x_(j) * sin(th(l)), &
                                 z(k), hel(j,k)
          end do
        end do
      end do
    else
      do j = 0, nx
        do l = 0, nt
          do k = 0, nz
            write(35,'(4e11.3)') x_(j) * cos(th(l)), x_(j) * sin(th(l)), &
                                 z(k), pn(j,k)
          end do
        end do
      end do
    end if

    write(35,*)

    !close(35)

    return
  end subroutine save_3d

  subroutine helicity(u_r, u_t, u_z, hel)
    !Calculate helicity u.(curl u)
    use parameters
    use ic_bc, only : s
    use derivs, only : deriv_x, deriv_z
    implicit none

    real    (r2), intent(in)  :: u_r(0:nx,0:nz), u_t(0:nx,0:nz), u_z(0:nx,0:nz)
    real    (r2), intent(out) :: hel(0:nx, 0:nz)
    integer (i1)              :: k
    real    (r2)              :: u_r_z(0:nx,0:nz), u_t_z(0:nx,0:nz), &
                                 u_t_x(0:nx,0:nz), u_z_x(0:nx,0:nz)

    call deriv_z(u_r, u_r_z)
    call deriv_z(u_t, u_t_z)
    call deriv_x(u_t, u_t_x)
    call deriv_x(u_z, u_z_x)

    do k = 0, nz
      hel(:,k) = 0.5_r2 * delx * (-u_r(:,k) * u_t_z(:,k) + &
                                   u_t(:,k) * u_r_z(:,k) - &
                                   u_t(:,k) * u_z_x(:,k) + &
                                   u_z(:,k) * u_t_x(:,k)) + &
                 (s(:) * u_t(:,k) * u_z(:,k)) / one_eta
    end do

    return
  end subroutine helicity

  subroutine save_surface(pn, v, zn, ur, uz, bn, jn, p, t)
    !Save surfaces of fields for use in gnuplot
    use parameters
    use ic_bc, only : x, z
    implicit none

    integer (i1), intent(in) :: p
    real    (r2), intent(in) :: t, pn(0:nx,0:nz), v(0:nx,0:nz), &
                                zn(0:nx,0:nz), bn(0:nx,0:nz), &
                                jn(0:nx,0:nz), ur(0:nx,0:nz), &
                                uz(0:nx,0:nz)
    integer (i1)             :: j, k

    open (19, status = 'unknown', file = 'p'//itos(p)//'.dat')
    !open (70, status = 'unknown', file = 'sin.dat')
    open (21, status = 'unknown', file = 'z'//itos(p)//'.dat')
    open (23, status = 'unknown', file = 'u'//itos(p)//'.dat')
    open (30, status = 'unknown', file = 'vr'//itos(p)//'.dat')
    open (31, status = 'unknown', file = 'vz'//itos(p)//'.dat')
    open (34, status = 'unknown', file = 'b'//itos(p)//'.dat')
    open (35, status = 'unknown', file = 'j'//itos(p)//'.dat')

    write(30, '(2A, i10, e19.7)') '#', 'p=', p, t
    write(31, '(2A, i10, e19.7)') '#', 'p=', p, t
    write(34, '(2A, i10, e19.7)') '#', 'p=', p, t
    write(35, '(2A, i10, e19.7)') '#', 'p=', p, t
    write(23, '(2A, i10, e19.7)') '#', 'p=', p, t
    write(19, '(2A, i10, e19.7)') '#', 'p=', p, t
    write(21, '(2A, i10, e19.7)') '#', 'p=', p, t
    !write(70, '(2A, i10, e19.7)') '#', 'p=', p, t

    do j = 0, nx
      write(30, '(3e19.7)') (x(j), z(k), ur(j,k), k = 0, nz)
      write(30, *)
      write(31, '(3e19.7)') (x(j), z(k), uz(j,k), k = 0, nz)
      write(31, *)
      write(34, '(3e19.7)') (x(j), z(k), bn(j,k), k = 0, nz)
      write(34, *)
      write(35, '(3e19.7)') (x(j), z(k), jn(j,k), k = 0, nz)
      write(35, *)
      write(23, '(3e19.7)') (x(j), z(k), v(j,k), k = 0, nz)
      write(23, *)
      write(19, '(3e19.7)') (x(j), z(k), pn(j,k), k = 0, nz)
      write(19, *)
      write(21, '(3e19.7)') (x(j), z(k), zn(j,k), k = 0, nz)
      write(21, *)
     !write(70, '(3e19.7)') (x(j), z(k), sin(pi*x(j))*&
     !                      sin(alpha*z(k)), k = 0, nz)
     !write(70, *)
    end do

    close (19)
    close (21)
    close (23)
    close (30)
    close (31)
    close (34)
    close (35)
    !close (70)

    return
  end subroutine save_surface

  subroutine write_data(p, p_start, t)
    !Calculations of various quantities; save fields
    use parameters
    use variables
    implicit none

    integer (i1), intent(in) :: p, p_start
    real    (r2), intent(in) :: t
    real    (r2)             :: growth_rate, growth_rate_vz, growth_diff, &
                                growth_error
    real    (r2), save       :: growth_rate_vz_old

    call vr_vz(psi%old, vr, vz)   !get radial, axial velocities
    if (save_part) then
      call particle(vr, vrold, vz, vzold) !save particle
    end if                                              !path
    call save_torque(t, ut%new)
    if ((p /= p_start) .and. ((p - p_start) > save_rate)) then
      call save_growth(t, vr, vrold, vz, vzold, psi%old, ut%new, zt%new, &
                       bt%old, jt%old, growth_rate, growth_rate_vz)
      growth_diff = abs(growth_rate_vz - growth_rate_vz_old)
      growth_error = growth_diff / abs(growth_rate_vz_old)
      !print*, Re1, growth_diff, growth_error
      if (auto_Re .and. ((growth_diff < growth_tol) .or. &
          (dec_Re .and. zero_Re .and. (vz(nx/2,nz/2) < 1E-6_r2)))) then
        call increment_Re(growth_rate_vz, t)
      end if
      !write (37, '(4e17.9)') t, growth_rate_vz_old, growth_rate_vz, growth_diff
      growth_rate_vz_old = growth_rate_vz
      call save_energy(vr, ut%new, vz, t)
      if ((abs(om1 - 0.0_r2) < epsilon(om1)) .and. &
          (abs(om2 - 0.0_r2) < epsilon(om2))) then
        if (hyst_Re) then
          if (((abs(growth_rate_vz) < 1E-6_r2) .and. &
               (abs(vz(nx/2,nz/2)) > 1E-3_r2))) then
            call increment_hysteresis(.true., t)
          else if ((abs(growth_rate_vz) > 1E-6_r2) .and. &
                       (growth_error < 1E-8_r2)) then
            call increment_hysteresis(.false., t)
          end if
        end if
        if ((abs(growth_rate_vz) < 1E-8_r2) .and. &  !if vr saturated
            (abs(vr(nx/2, nz/2)) > 1E-3_r2)) then
          saturated = saturated + 1
          if (saturated > 4) then
            if ((.not. auto_tau) .or. (tau == tau_end)) then  !if tau not auto
              call save_time_tau(tau, t)                     !or tau at end
              call end_state(ut%old, zt%old, psi%old, bt%old, jt%old, p) !finish
            else if (tau < 1.0_r2) then
              call save_time_tau(tau, t)
              tau = tau + tau_step   !increment tau
              print*, 'tau = ', tau
              call save_xsect(vr, vz, psi%old, ut%new, zt%new, &
                              bt%old, jt%old, t, p)
              call save_surface(psi%old, ut%new, zt%new, vr, vz, &
                                bt%old, jt%old, p, t)
            end if
          end if
        else
          saturated = 0
        end if
      end if
    end if

    return
  end subroutine write_data

  subroutine terminate(p, t)
    !Terminate run if file 'RUNNING' does not exist in run directory
    use parameters
    use variables
    implicit none

    integer (i1), intent(in) :: p
    real    (r2), intent(in) :: t
    logical                  :: run_exist

    if (mycol == 0) then
      inquire(file='RUNNING', exist=run_exist)  !does 'RUNNING' exist?
      if (.not. run_exist) then  !if not then finish
        print*
        write(6, '(A33, i2)') 'Stop requested.  Ending process ', mycol
        print*, 'Saving end state...'
        call end_state(ut%old, zt%old, psi%old, bt%old, jt%old, p)
        call save_xsect(vr, vz, psi%old, ut%new, zt%new, &
                        bt%old, jt%old, t, p)
        call save_surface(psi%old, ut%old, zt%old, &
                          vr, vz, bt%old, jt%old, p, t)
        end_proc = 1  !flag to send to all processes
      end if
    end if

    call SLTIMER(3)                                         !tell all other
    call IGEBR2D(ictxt, 'A', ' ', 1, 1, end_proc, 1, 0, 0)  !processes to finish
    call SLTIMER(3)                                         !if end_proc=1

    return
  end subroutine terminate

  subroutine save_run(p, t)
    !Save cross-sections, surfaces if file 'save' exists in run directory
    use parameters
    use variables
    implicit none

    integer (i1), intent(in) :: p
    real    (r2), intent(in) :: t
    logical                  :: save_exist

    inquire(file='save', exist=save_exist)   !does 'SAVE' exist?
    if (save_exist) then   !if so then save fields
      call save_xsect(vr, vz, psi%old, ut%new, zt%new, &
                      bt%old, jt%old, t, p)
      call save_surface(psi%old, ut%old, zt%old, &
                        vr, vz, bt%old, jt%old, p, t)
      open (98, file = 'save')
      close (98, status = 'delete')
    end if

    return
  end subroutine save_run

  subroutine end_state(u, zn, pn, bn, jn, p)
    !Save variables for use in a restarted run
    use parameters
    implicit none

    integer (i1), intent(in) :: p
    real    (r2), intent(in) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                jn(0:nx,0:nz)
    integer (i1)             :: j, k

    open (50, file = 'end_state.dat', form='unformatted')

    write(50) nx
    write(50) nz
    write(50) p
    write(50) dt
    write(50) u
    write(50) zn
    write(50) pn
    write(50) bn
    write(50) jn

    close (50)

    open (99, file = 'RUNNING')  !delete 'RUNNING' to finish run
    close (99, status = 'delete')

    return
  end subroutine end_state

  subroutine save_time_tau (tau, t)
    !Save the times at which tau is incremented
    use parameters, only : r2
    implicit none

    real (r2), intent(in) :: tau, t

    write(51, '(2e17.9)') t, tau

    return
  end subroutine save_time_tau

  subroutine particle (vr, vrold, vz, vzold)
    !Calculate and save a particle path
    use parameters
    use variables, only : xold, zold
    implicit none

    real    (r2), intent(in)    :: vr(0:nx, 0:nz), vz(0:nx, 0:nz), &
                                   vrold(0:nx, 0:nz), vzold(0:nx, 0:nz)
    integer (i1)                :: xmin(num_pars), xplu(num_pars), &
                                   zmin(num_pars), zplu(num_pars), i, j
    real    (r2)                :: cc1(num_pars), cc2(num_pars), &
                                   rvel(num_pars), zvel(num_pars), &
                                   xnew(num_pars), znew(num_pars), del_t

    del_t = dt / 1.0_r2

    do j = 1, 1
      xmin = INT(xold,i1)
      xplu = INT(xold + 1.0_r2,i1)
      zmin = INT(zold,i1)
      zplu = INT(zold + 1.0_r2,i1)

      cc1 = (xold - xmin) / (xplu - xmin)
      cc2 = (zold - zmin) / (zplu - zmin)

      do i = 1, num_pars
        rvel(i) = (1.0_r2 - cc1(i)) * (1.0_r2 - cc2(i)) * &
                  vrold(xmin(i), zmin(i)) + &
                  cc1(i)* (1.0_r2 - cc2(i)) * vrold(xplu(i), zmin(i)) + &
                  cc1(i) * cc2(i) * vrold(xplu(i), zplu(i)) + &
                  (1.0_r2 - cc1(i)) * cc2(i) * vrold(xmin(i), zplu(i))
      end do

      do i = 1, num_pars
        zvel(i) = (1.0_r2 - cc1(i)) * (1.0_r2 - cc2(i)) * &
                  vzold(xmin(i), zmin(i)) + &
                  cc1(i) * (1.0_r2 - cc2(i)) * vzold(xplu(i), zmin(i)) + &
                  cc1(i) * cc2(i) * vzold(xplu(i), zplu(i)) + &
                  (1.0_r2 - cc1(i)) * cc2(i) * vzold(xmin(i), zplu(i))
      end do

      xnew = xold + del_t * rvel
      znew = zold + del_t * zvel
      xold = xnew
      zold = znew
    end do

    write (24, '(10e17.9)') xnew(1) / nx, znew(1) * gamma / nz, &
                            xnew(2) / nx, znew(2) * gamma / nz, &
                            xnew(3) / nx, znew(3) * gamma / nz, &
                            xnew(4) / nx, znew(4) * gamma / nz, &
                            xnew(5) / nx, znew(5) * gamma / nz

    return
  end subroutine particle

  subroutine energy_CCF(t, Eccf)
    !Explicitly integrated (dimensionless) form of the kinetic energy in CCF
    use parameters
    use ic_bc, only : Re_1, Re_2
    implicit none

    real (r2), intent(in)  :: t
    real (r2), intent(out) :: Eccf
    real (r2)              :: a, b, c

    a = 0.25_r2 * (Re_2(t) - eta * Re_1(t))**2 * (1.0_r2 - eta**4)
    b = (Re_2(t) - eta * Re_1(t)) * (Re_1(t) - eta * Re_2(t)) * &
        eta * (1.0_r2 - eta**2)
    c = eta**2 * (Re_1(t) - eta * Re_2(t))**2 * log(1.0_r2 / eta)

    Eccf = gamma * pi * (a + b + c) / &
          ((1.0_r2 + eta)**2 * (one_eta)**4)
    
    return
  end subroutine energy_CCF

  subroutine save_energy(ur, ut, uz, t)
    !Total kinetic energy in flow (including CCF)
    use parameters
    use variables, only : integrate_r, integrate_z
    implicit none

    real (r2), intent(in)  :: ur(0:nx,0:nz), ut(0:nx,0:nz), uz(0:nx,0:nz), t
    real (r2)              :: u2(0:nx,0:nz), int_r(0:nz), int_z, E, Eccf

    u2 = ur**2 + ut**2 + uz**2

    call energy_CCF(t, Eccf)
    call integrate_r(u2, int_r)
    call integrate_z(int_r, int_z)

    E = pi * int_z

    write (36, '(3e17.9)') t, Eccf, E

    return
  end subroutine save_energy

  subroutine increment_Re(growth_rate, t)
    !Automatically increment/decrement the Reynolds number to find critical
    !values based on the sign of the growth rate
    use parameters
    implicit none

    real (r2), intent(in) :: growth_rate, t
    real (r2), save       :: Re_minus, Re_plus, growth_minus, growth_plus

    select case (dec_RE)
      case (.false.)
        if (init_Re) then
          if (growth_rate < 0.0_r2) then
            growth_minus = growth_rate
            Re_minus = Re1
            Re1 = Re1 + Re_incr
            gm_set = .true.
          else if (growth_rate > 0.0_r2) then
            growth_plus = growth_rate
            Re_plus = Re1
            Re1 = Re1 - Re_incr
            gp_set = .true.
          end if
          init_Re = .false.
        else if (.not. init_Re) then
          if ((.not. gm_set) .or. (.not. gp_set)) then
            if (.not. gp_set) then
              if (growth_rate < 0.0_r2) then
                growth_minus = growth_rate
                Re_minus = Re1
                Re1 = Re1 + Re_incr
              else if (growth_rate > 0.0_r2) then
                growth_plus = growth_rate
                Re_plus = Re1
                Re1 = 0.5_r2 * (Re1 + Re_minus)
                gp_set = .true.
              end if
            else if (.not. gm_set) then
              if (growth_rate > 0.0_r2) then
                growth_plus = growth_rate
                Re_plus = Re1
                Re1 = Re1 - Re_incr
              else if (growth_rate < 0.0_r2) then
                growth_minus = growth_rate
                Re_minus = Re1
                Re1 = 0.5_r2 * (Re1 + Re_plus)
                gm_set = .true.
              end if
            end if
          else
            if (growth_rate < 0.0_r2) then
              growth_minus = growth_rate
              Re_minus = Re1
              Re1 = 0.5_r2 * (Re1 + Re_plus)
            else if (growth_rate > 0.0_r2) then
              growth_plus = growth_rate
              Re_plus = Re1
              Re1 = 0.5_r2 * (Re1 + Re_minus)
            end if
          end if
        end if
      case (.true.)
        if (init_Re) then
          if (growth_rate < 0.0_r2) then
            growth_minus = growth_rate
            Re_minus = Re1
            Re1 = Re1 - Re_incr
            gm_set = .true.
          else if (growth_rate > 0.0_r2) then
            stop 'Decreasing Re boundary only but growth rate > 0'
          end if
          init_Re = .false.
        else if (.not. init_Re) then
          if (.not. gp_set) then
            if (growth_rate < 0.0_r2) then
              growth_minus = growth_rate
              Re_minus = Re1
              Re1 = Re1 - Re_incr
            else if (growth_rate > 0.0_r2) then
              growth_plus = growth_rate
              Re_plus = Re1
              Re1 = 0.0_r2
              zero_Re = .true.
              gp_set = .true.
            end if
          else
            if (zero_Re) then
              Re1 = 0.5_r2 * (Re_minus + Re_plus)
              zero_Re = .false.
            else if (.not. zero_Re) then
              if (growth_rate < 0.0_r2) then
                growth_minus = growth_rate
                Re_minus = Re1
                Re1 = 0.5_r2 * (Re1 + Re_plus)
              else if (growth_rate > 0.0_r2) then
                growth_plus = growth_rate
                Re_plus = Re1
                Re1 = 0.0_r2
                zero_Re = .true.
              end if
            end if
          end if
        end if
    end select
    
    write (37, '(6e17.9)') t, Re1, Re_minus, Re_plus, growth_minus, growth_plus
    
    if (abs(Re1 - Re_minus) < 1E-6_r2) then
      open (99, file = 'RUNNING')
      close (99, status = 'delete')
    end if
    
    return
  end subroutine
  
  subroutine increment_hysteresis(test, t)
    !Automatically increment/decrement the Reynolds number to find critical
    !values in a hysteresis region
    use parameters
    use variables
    implicit none

    real     (r2), intent(in) :: t
    logical,       intent(in) :: test
    real     (r2), save       :: Re_minus, Re_plus

    if (init_Re) then
      if (test) then
        Re_plus = Re1
        Re1 = Re1 - Re_incr
      end if
      init_Re = .false.
    else if (.not. init_Re) then
      if (test) then
        if (.not. gm_set) then
          Re_plus = Re1
          Re1 = Re1 - Re_incr
        else if (gm_set) then
          Re_plus = Re1
          Re1 = 0.5_r2 * (Re_minus + Re_plus)
        end if
      else if (.not. test) then
        Re_minus = Re1
        Re1 = 0.5_r2 * (Re_minus + Re_plus)
        gm_set = .true.
        open (50, file = 'end_state.dat', form = 'unformatted')
      
        read (50)
        read (50)
        read (50)
        read (50)
        read (50) ut%new
        read (50) zt%new
        read (50) psi%new
        read (50) bt%new
        read (50) jt%new

        close (50)
      end if
    end if

    write (37, '(4e17.9)') t, Re1, Re_minus, Re_plus
    
    if (abs(Re1 - Re_minus) < 1E-6_r2) then
      open (99, file = 'RUNNING')
      close (99, status = 'delete')
    end if

    return
  end subroutine increment_hysteresis
  
  subroutine get_timestep()
    !Obsolete: from attempt at adaptive timestep
    use parameters
    implicit none

    real (r2) :: timestep

    if (Re1 * Re2 >= 0.0_r2) then
      timestep = max(1.0_r2, abs(Re1) + abs(Re1_mod), &
                 abs(Re2) + abs(Re2_mod))
    else
      timestep = max(1.0_r2, abs(Re1) + abs(Re1_mod) + &
                 abs(Re2) + abs(Re2_mod))
    end if

    timestep = 1e-4_r2 * eta * 2.0_r2 * pi / (one_eta * timestep)

    print*, timestep

    return
  end subroutine get_timestep

end module io
