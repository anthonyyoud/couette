module ic_bc
  !Routines to do with grid setup, initial and boundary conditions.
  use parameters
  implicit none

  private
  public :: get_xzs, ICS, u_BCS, z_BCS, p_BCS, b_BCS, j_BCS, Re_1, Re_2, &
            state_restart

  real (r2), public :: x(0:nx), x_(0:nx), th(0:nt), z(0:nz), s(0:nx) 
                                                !finite-difference mesh
                                                !s=eta+(1-eta)*x
  contains

  subroutine get_xzs()
    !Finite-difference mesh
    use parameters
    implicit none

    integer (i1) :: j, k, l

    do k = 0, nz
      z(k) = real(k,r2) * delz
    end do

    do l = 0, nt
      th(l) = real(l,r2) * delt
    end do

    do j = 0, nx
      x(j) = real(j,r2) * delx
    end do

    x_ = x + 1.0_r2                 !shift radial coordinate for OpenDX
    s = eta + one_eta * x

    return
  end subroutine get_xzs

  subroutine ICS(u, zn, pn, bn, jn, p)
    !Initial conditions
    use parameters
    implicit none

    integer (i1), intent(out)   :: p
    real    (r2), intent(inout) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                   pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                   jn(0:nx,0:nz)
    integer (i1)                :: j, k
    logical                     :: state_exist

    if (restart) then
      inquire(file='end_state.dat', exist=state_exist)
      !exit if doing restart but end_state.dat does not exist
      if (.not. state_exist) stop 'ERROR: restart=.true.&
                                   &but end_state.dat does not exist.'
      print*, 'Getting restart conditions'
      call state_restart(u, zn, pn, bn, jn, p)  !get saved data if restart
    else                   !put in estimate of eigen-function shape
      if (abs(tau - 1.0_r2) < epsilon(tau)) then !which satisfies BCS
                                                 !multiplied by small seed
        do k = 0, nz
          u(:,k) = seed * sin(2.0_r2*pi*z(k)/gamma) * sin(pi*x(:))
          pn(:,k) = seed * sin(2.0_r2*pi*z(k)/gamma) * sin(pi*x(:))
          bn(:,k) = seed * cos(2.0_r2*pi*z(k)/gamma) / s(:)
          jn(:,k) = seed * sin(pi*x(:)) * sin(2.0_r2*pi*z(k)/gamma)
        end do
      else
        do k = 0, nz
          u(:,k) = seed * sin(pi*x(:)) * cos(2.0_r2*pi*z(k)/gamma)
          pn(:,k) = seed * sin(pi*x(:)) * sin(2.0_r2*pi*z(k)/gamma)
          bn(:,k) = seed * sin(2.0_r2*pi*z(k)/gamma) / s(:)
          jn(:,k) = seed * sin(pi*x(:)) * cos(2.0_r2*pi*z(k)/gamma)
        end do

        do k = 1, nz1
          do j = 1, nx1
            zn(j,k) = -(pn(j+1,k) - 2.0_r2 * pn(j,k) + pn(j-1,k)) / &
                       (s(j) * dx2) + &
                       0.5_r2 * one_eta * (pn(j+1,k) - pn(j-1,k)) / &
                       (s(j)**2 * delx) - &
                       (pn(j,k+1) - 2.0_r2 * pn(j,k) + pn(j,k-1)) / &
                       (s(j) * dz2)  !ICS based on above fields for vorticity
          end do
        end do
      end if
    end if

    return
  end subroutine ICS

  subroutine state_restart(u, zn, pn, bn, jn, p)
    !Get restart data
    use parameters
    implicit none

    integer (i1), intent(out) :: p
    real    (r2), intent(out) :: u(0:nx,0:nz), zn(0:nx,0:nz), pn(0:nx,0:nz), &
                                 bn(0:nx,0:nz), jn(0:nx,0:nz)
    integer (i1)              :: j, k, nx_prev, nz_prev, alloc_err
    real    (r2)              :: dt_prev
    real    (r2), allocatable :: u_prev(:,:), z_prev(:,:), p_prev(:,:), &
                                 b_prev(:,:), j_prev(:,:)

    open (50, file = 'end_state.dat', form='unformatted')

    read (50) nx_prev
    read (50) nz_prev
    read (50) p
    read (50) dt_prev

    if ((nx_prev /= nx) .or. (nz_prev /= nz)) then  !interpolate onto new grid
      print*, 'Interpolating onto new grid...'
      allocate(u_prev(0:nx_prev,0:nz_prev), z_prev(0:nx_prev,0:nz_prev), &
               p_prev(0:nx_prev,0:nz_prev), b_prev(0:nx_prev,0:nz_prev) , &
               j_prev(0:nx_prev,0:nz_prev), stat=alloc_err)
      if (alloc_err /= 0) stop 'ERROR: Interpolating allocation error' 
      read (50) u_prev
      read (50) z_prev
      read (50) p_prev
      read (50) b_prev
      read (50) j_prev
      call inter(u_prev, z_prev, p_prev, b_prev, j_prev, nx_prev, nz_prev, &
                 u, zn, pn, bn, jn)
      if (allocated(u_prev)) deallocate(u_prev, stat=alloc_err)
      if (alloc_err /= 0) stop 'ERROR: Interpolating deallocation error'
      if (allocated(z_prev)) deallocate(z_prev, stat=alloc_err)
      if (alloc_err /= 0) stop 'ERROR: Interpolating deallocation error'
      if (allocated(p_prev)) deallocate(p_prev, stat=alloc_err)
      if (alloc_err /= 0) stop 'ERROR: Interpolating deallocation error'
      if (allocated(b_prev)) deallocate(b_prev, stat=alloc_err)
      if (alloc_err /= 0) stop 'ERROR: Interpolating deallocation error'
      if (allocated(j_prev)) deallocate(j_prev, stat=alloc_err)
      if (alloc_err /= 0) stop 'ERROR: Interpolating deallocation error'
    else     !just read the data
      read (50) u
      read (50) zn
      read (50) pn
      read (50) bn
      read (50) jn
    end if

    close (50)   

    if (dt_prev /= dt) then  !if restart tstep /= old tstep adjust index 'p'
      p = p * dt_prev / dt
    end if

    return
  end subroutine state_restart

  subroutine inter(u_, z_, p_, b_, j_, nxp, nzp, &
                   u, zn, pn, bn, jn)
    !Setup for interpolating fields from smaller grid onto larger grid
    use parameters
    implicit none

    integer (i1), intent(in)  :: nxp, nzp
    real    (r2), intent(in)  :: u_(0:nxp,0:nzp), z_(0:nxp,0:nzp), &
                                 p_(0:nxp,0:nzp), b_(0:nxp,0:nzp), &
                                 j_(0:nxp,0:nzp)
    real    (r2), intent(out) :: u(0:nx,0:nz), zn(0:nx,0:nz), &
                                 pn(0:nx,0:nz), bn(0:nx,0:nz), &
                                 jn(0:nx,0:nz)
    integer (i1)              :: j, k
    real    (r2)              :: dx_prev, dz_prev, x_prev(0:nxp), z_prev(0:nzp)

    dx_prev = 1.0_r2 / nxp       !previous space mesh
    dz_prev = gamma / nzp

    do j = 0, nxp
      x_prev(j) = real(j,r2) * dx_prev
    end do
                                  !previous coordinates
    do k = 0, nzp
      z_prev(k) = real(k,r2) * dz_prev
    end do

    call inter_var(u_, x_prev, z_prev, nxp, nzp, u)
    call inter_var(z_, x_prev, z_prev, nxp, nzp, zn)
    call inter_var(p_, x_prev, z_prev, nxp, nzp, pn)    !interpolate
    call inter_var(b_, x_prev, z_prev, nxp, nzp, bn)
    call inter_var(j_, x_prev, z_prev, nxp, nzp, jn)

    return
  end subroutine inter

  subroutine inter_var(in_var, x_prev, z_prev, nxp, nzp, out_var)
    !Bilinearly interpolate the actual fields onto a new grid
    !See Numerical Recipes in Fortran77 Chap. 3.6 p.116
    use parameters
    implicit none
    
    integer (i1), intent(in)  :: nxp, nzp
    real    (r2), intent(in)  :: in_var(0:nxp,0:nzp), x_prev(0:nxp), &
                                 z_prev(0:nzp)
    real    (r2), intent(out) :: out_var(0:nx,0:nz)
    integer (i1)              :: j, k, j2, k2
    real    (r2)              :: int1, int2

    do k = 0, nz   !new 'z' index
      k2 = INT(nzp*k/nz,i1)   !old 'z' index
      do j = 0, nx   !new 'x' index
        j2 = INT(nxp*j/nx,i1)   !old 'x' index
        int1 = (x(j) - x_prev(j2)) / (x_prev(j2+1) - x_prev(j2))  !interpolating
        int2 = (z(k) - z_prev(k2)) / (z_prev(k2+1) - z_prev(k2))  !constants
        out_var(j,k) = (1.0_r2 - int1) * (1.0_r2 - int2) * in_var(j2,k2) + &
                       int1 * (1.0_r2 - int2) * in_var(j2+1,k2) + &
                       int1 * int2 * in_var(j2+1,k2+1) + &   !bilinear
                       (1.0_r2 - int1) * int2 * in_var(j2,k2+1) !interpolation
      end do
    end do

    return
  end subroutine inter_var

  subroutine u_BCS(u, t)
    !Boundary conditions for total azimuthal velocity (including CCF)
    use parameters
    implicit none

    real    (r2), intent(in)  :: t
    real    (r2), intent(out) :: u(0:nx,0:nz)
    integer (i1)              :: k

    !u(0,:) = Re1 + Re1_mod * cos(om1 * t) + &
    !         eps1 * Re1 * (1.0_r2 / eta - 1.0_r2) * cos(freq1 * z(:))
    !u(nx,:) = Re2 + Re2_mod * cos(om2 * t) - &
    !         eps2 * cos(freq2 * z(:))

    u(0,:) = Re_1(t)
    u(nx,:) = Re_2(t)
    
    if (abs(tau - 1.0_r2) < epsilon(tau)) then
      if (rot_ends) then
        u(:,0) = Re_1(t) * s(:) / eta
        u(:,nz) = 0.0_r2 !Re_1(t) * s(:) / eta
      else
        u(:,0) = 0.0_r2
        u(:,nz) = 0.0_r2
      end if
    end if

    return
  end subroutine u_BCS

  subroutine z_BCS(zn, pn, t)
    !Boundary conditions for azimuthal vorticity
    use parameters
    implicit none

    real (r2), intent(in)  :: t, pn(0:nx,0:nz)
    real (r2), intent(out) :: zn(0:nx,0:nz)

    zn(0,:) = -(8.0_r2 * pn(1,:) - pn(2,:)) / (2.0_r2 * s(0) * dx2)
    zn(nx,:) = -(8.0_r2 * pn(nx1,:) - pn(nx-2,:)) / (2.0_r2 * s(nx) * dx2)

    if (abs(tau - 1.0_r2) < epsilon(tau)) then
      zn(:,0) = -(8.0_r2 * pn(:,1) - pn(:,2)) / &
                 (2.0_r2 * (s(:)) * dz2)
      zn(:,nz) = -(8.0_r2 * pn(:,nz1) - pn(:,nz-2)) / &
                  (2.0_r2 * (s(:)) * dz2)
    else
      zn(:,0) = (-tau / (s(:) * (1.0_r2 - tau))) * &
                (0.5_r2 * (-pn(:,2) + 4.0_r2 * pn(:,1)) / delz)
      zn(:,nz) = (tau / (s(:) * (1.0_r2 - tau))) * &
                (0.5_r2 * (pn(:,nz-2) - 4.0_r2 * pn(:,nz1)) / delz)
    end if

    return
  end subroutine z_BCS

  subroutine p_BCS(p)
    !Boundary conditions for stream-function, psi
    use parameters
    implicit none

    real (r2), intent(out) :: p(0:nx,0:nz)

    p(0,:) = 0.0_r2
    p(nx,:) = 0.0_r2

    p(:,0) = 0.0_r2
    p(:,nz) = 0.0_r2

    return
  end subroutine p_BCS

  subroutine b_BCS(bn)
    !Boundary conditions for azimuthal magnetic field
    use parameters
    implicit none

    real (r2), intent(out) :: bn(0:nx,0:nz)

    if (abs(tau - 0.0_r2) < epsilon(tau)) then
      bn(:,0) = 0.0_r2
      bn(:,nz) = 0.0_r2
    end if

    return
  end subroutine b_BCS

  subroutine j_BCS(jn)
    !Boundary conditions for azimuthal current
    use parameters
    implicit none

    real (r2), intent(out) :: jn(0:nx,0:nz)

    jn(0,:) = 0.0_r2
    jn(nx,:) = 0.0_r2

    if (abs(tau - 1.0_r2) < epsilon(tau)) then
      jn(:,0) = 0.0_r2
      jn(:,nz) = 0.0_r2
    end if

    return
  end subroutine j_BCS

  function Re_1(t)
    !Time-dependent Reynolds number of the inner cylinder
    use parameters
    implicit none
                      
    real (r2), intent(in) :: t
    real (r2)             :: Re_1 

    Re_1 = Re1 + Re1_mod * cos(om1 * t)
                                                          
    return                    
  end function Re_1

  function Re_2(t)
    !Time-dependent Reynolds number of the outer cylinder
    use parameters
    implicit none

    real (r2), intent(in) :: t
    real (r2)             :: Re_2

    Re_2 = Re2 + Re2_mod * cos(om2 * t)

    return
  end function Re_2

  function f1(index)
    use parameters
    implicit none

    integer (i1), intent(in) :: index
    real    (r2)             :: f1

    f1 = eps1 * cos(freq1*z(index))

    return
  end function f1   

  function f2(index)
    use parameters
    implicit none

    integer (i1), intent(in) :: index
    real    (r2)             :: f2

    f2 = eps2 * cos(freq2*z(index)-pi)

    return
  end function f2 

end module ic_bc
