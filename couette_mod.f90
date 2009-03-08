program couette_mod
  !Main program file.  See README for a description of the code and how to run
  !it.
  use parameters
  use stream
  use matrices
  use io
  use ic_bc
  use variables
  use linear
  use nonlinear
  use solve
  use magnetic
  use current
  implicit none
  
  integer :: j, k, p=0, p_start=0, alloc_err=0
  integer, dimension(nx1*nz1) :: p_pivot
  integer, allocatable :: b_pivot(:), j_pivot(:)
  real :: t=0.0
  real, allocatable :: b_mat(:,:), j_mat(:,:)
  real, dimension(2*nx1+nx1+1, nx1*nz1) :: p_mat
  type (mat_comp) :: Ux, Zx
  type (uz_mat_comp) :: Uz
  type (zz_mat_comp) :: Zz
  logical :: state_exist
 
  if (abs(tau - 0.0) < epsilon(tau)) then
    write(6, '(A7, f4.2, A20)') 'tau = ', tau, '- Infinite cylinder'
  else if (abs(tau - 1.0) < epsilon(tau)) then
    write(6, '(A7, f4.2, A22)') 'tau = ', tau, '- Finite aspect ratio'
  else if ((abs(tau - 0.0) > epsilon(tau)) .and. &
           (abs(tau - 1.0) > epsilon(tau))) then
    write(6, '(A7, f4.2, A15)') 'tau = ', tau, '- Variable tau'
  else
    write(6, '(A7, f4.2)') 'tau = ', tau
    stop 'ERROR: Invalid value of tau'
  end if   

  print*, 'Setting up finite-difference mesh...'

  call get_xzs()   !get finite-difference mesh

  print*, 'Setting up ICS...'

  inquire(file='RUNNING', exist=state_exist)
  if (state_exist) then
    print*, 'File RUNNING exists in the current directory.  Assuming &
            &a restart is required...'
    call state_restart(ut%new, zt%new, psi%new, bt%new, jt%new, p_start)
  else
    call ICS(ut%new, zt%new, psi%new, bt%new, jt%new, p_start)

    if (.not. restart) then
      inquire(file='end_state.dat', exist=state_exist)
      !exit if not doing restart but end_state.dat exists
      if (state_exist) stop 'ERROR: restart=.false. but end_state.dat exists.'
      print*, 'Setting up BCS...'
      call u_BCS(ut%new, 0.0)
      call p_BCS(psi%new)
      call z_BCS(zt%new, psi%new, 0.0)   !get boundary conditions
      call b_BCS(bt%new)
      call j_BCS(jt%new)
    end if
  end if

  call particle_setup()

  print*, 'Opening runtime files...'

  call open_files()   !open files for writing

  print*, 'Allocating matrix dimensions...'

  if (abs(tau - 0.0) < epsilon(tau)) then
    allocate(b_mat(2*nxp1+nxp1+1,0:nxp1*nz1-1), stat=alloc_err)
    if (alloc_err /= 0) stop 'ERROR: b_mat allocation error'
    allocate(j_mat(2*nx1+nx1+1,nx1*nzp1), stat=alloc_err)
    if (alloc_err /= 0) stop 'ERROR: j_mat allocation error'
    allocate(b_pivot(nxp1*nz1))
    allocate(j_pivot(nx1*nzp1))
  else if (abs(tau - 1.0) < epsilon(tau)) then   !or infinite cylinders
    allocate(b_mat(2*nxp1+nxp1+1,0:nxp1*nzp1-1), stat=alloc_err)
    if (alloc_err /= 0) stop 'ERROR: b_mat allocation error'
    allocate(j_mat(2*nx1+nx1+1,nx1*nz1), stat=alloc_err)
    if (alloc_err /= 0) stop 'ERROR: j_mat allocation error'
    allocate(b_pivot(nxp1*nzp1))
    allocate(j_pivot(nx1*nz1))
  else if ((abs(tau - 0.0) > epsilon(tau)) .and. &
           (abs(tau - 1.0) > epsilon(tau))) then
    allocate(b_mat(2*nxp1+nxp1+1,0:nxp1*nzp1-1), stat=alloc_err)
    if (alloc_err /= 0) stop 'ERROR: b_mat allocation error'
    allocate(j_mat(2*nx1+nx1+1,nx1*nzp1), stat=alloc_err)
    if (alloc_err /= 0) stop 'ERROR: j_mat allocation error'
    allocate(b_pivot(nxp1*nzp1))
    allocate(j_pivot(nx1*nzp1))
  else
    stop 'ERROR: Allocation error due to value of tau'
  end if

  print*, 'Setting up matrices...'   !matrices for tridiagonal systems
  call matrix_setup(Ux, Uz, Zx, Zz)  !azimuthal velocity, vorticity

  call psi_mat_setup(p_mat, p_pivot)

  if (abs(tau - 0.0) < epsilon(tau)) then
    call b_mat_setup(b_mat, b_pivot)
    call j_mat_setup(j_mat, j_pivot)
  else if (abs(tau - 1.0) < epsilon(tau)) then     !left-hand side matrices
    call fin_b_mat_setup(b_mat, b_pivot)              !for stream-function,
    call fin_j_mat_setup(j_mat, j_pivot)              !current and magnetic
  else if ((abs(tau - 0.0) > epsilon(tau)) .and. & !Poisson equations
           (abs(tau - 1.0) > epsilon(tau))) then
    call fin_b_mat_setup(b_mat, b_pivot)
    call j_mat_setup(j_mat, j_pivot)
  else
    stop 'ERROR: Matrix setup subroutine call error due to value of tau'
  end if

  call copy_var(ut%old, ut%new)
  call copy_var(ut%old2, ut%new)
  call copy_var(ut%inter, ut%new)
  call copy_var(zt%old, zt%new)
  call copy_var(zt%old2, zt%new)
  call copy_var(zt%inter, zt%new)   !initialise all variables equal to
  call copy_var(psi%old, psi%new)   !'new' from ICS and BCS
  call copy_var(psi%old2, psi%new)
  call copy_var(bt%old, bt%new)
  call copy_var(bt%old, bt%new)
  call copy_var(jt%old, jt%new)
  call copy_var(jt%old, jt%new)

  print*, 'Entering time loop...'

  do p = p_start, Ntot        !start main time loop
    call terminate(p, t)      !loop terminates if file 'RUNNING'
                              !does not exist in run directory
    if (end_proc == 1) exit   !end master process

                             !save cross-sections, surfaces if
    call save_run(p, t)      !file 'save' exists in run directory

    t = p * dt   !increment time

    if (mod(p, save_rate) == 0) call write_data(p, p_start, t)

    !print*, vr(nx/2,nz/2), vrold(nx/2,nz/2)
                                                  !do necessary writing
    call copy_var(vrold, vr)                      !of data to files
    call copy_var(vzold, vz)   !update variables

    if (xsect_save) then   !save cross-sections
      if (mod(p, save_rate_2) == 0) then
        call save_xsect(vr, vz, psi%old, ut%new, zt%new, &
                        bt%old, jt%old, t, p)
        !call save_surface(psi%old, ut%new, zt%new, vr, vz, &
        !                  bt%old, jt%old, p, t)
      end if
    end if

    if (save3d) then   !save 3D surface
      if (mod(p, save_rate_2) == 0) then
        call save_3d(vr, ut%new, vz, psi%old, p)
        call save_vapor_3d(vr, ut%new, vz, psi%old, p)
      end if
    end if

    call copy_var(ut%old, ut%new)
    call copy_var(zt%old, zt%new)   !update variables

    !get linear part for u in x-direction
    call get_rhs_ux(ut%old, ut%new)
    !get nonlinear part for u in x-direction
    call get_nlin_ux(ut%old, ut%old2, psi%old, psi%old2, &
                     bt%old, bt%old2, ut%nlin_new)
    !get linear part for Z in x-direction
    call get_rhs_Zx(zt%old, zt%new)
    !get nonlinear part for Z in x-direction
    call get_nlin_Zx(t, ut%old, ut%old2, psi%old, psi%old2, zt%old, &
                     zt%old2, jt%old, jt%old2, zt%nlin_new)
    !solve for u in x-direction
    call solve_ux(ut%old, ut%new, ut%nlin_new, t, Ux)
    !solve for Z in x-direction
    call solve_Zx(zt%old, zt%new, zt%nlin_new, psi%old, t, Zx)

    call copy_var(ut%old2, ut%inter)
    call copy_var(zt%old2, zt%inter)   !update variables

    !solve for u in z-direction
    call solve_uz(ut%old, ut%new, t, Uz)
    !solve for Z in z-direction
    call solve_Zz(zt%old, psi%old, zt%new, t, Zz)

    call copy_var(ut%inter, ut%new)
    call copy_var(zt%inter, zt%new)
    call copy_var(psi%old2, psi%old)   !update variables
    call copy_var(bt%old2, bt%old)
    call copy_var(jt%old2, jt%old)

    call p_BCS(psi%new)
    call b_BCS(bt%new)   !update boundary conditions
    call j_BCS(jt%new)
   
    !solve Poisson equations depending on tau
    call p_poisson(zt%new, psi%new, p_mat, p_pivot)
    if (abs(tau - 0.0) < epsilon(tau)) then
      call b_poisson(ut%new, bt%new, b_mat, b_pivot)
      call j_poisson(psi%new, jt%new, j_mat, j_pivot)
    else if (abs(tau - 1.0) < epsilon(tau)) then
      call fin_b_poisson(ut%new, bt%new, b_mat, b_pivot)
      call fin_j_poisson(psi%new, jt%new, j_mat, j_pivot)
    else if ((abs(tau - 0.0) > epsilon(tau)) .and. &
             (abs(tau - 1.0) > epsilon(tau))) then
      call fin_b_poisson(ut%new, bt%new, b_mat, b_pivot)
      call j_poisson(psi%new, jt%new, j_mat, j_pivot)
    else
      stop 'ERROR: Poisson solver subroutine call error due to value of tau'
    end if

    call copy_var(psi%old, psi%new)
    call copy_var(bt%old, bt%new)   !update variables
    call copy_var(jt%old, jt%new)

    if (mod(p, periodic_save) == 0) then
      call end_state(ut%new, zt%new, psi%old, bt%old, jt%old, p, 0)
    end if

    if (p == Ntot) then
      call end_state(ut%new, zt%new, psi%old, bt%old, jt%old, p, 1)
      call save_xsect(vr, vz, psi%old, ut%new, zt%new, &  !save final state
                      bt%old, jt%old, t, p)
      call save_surface(psi%old, ut%new, zt%new, &
                        vr, vz, bt%old, jt%old, p, t)
    end if
  end do   !end time loop

  print*, 'Deallocating allocated arrays...'

  if (allocated(b_mat)) deallocate(b_mat, stat=alloc_err) !deallocate arrays
  if (alloc_err /= 0) stop 'ERROR: b_mat deallocation error'
  if (allocated(j_mat)) deallocate(j_mat, stat=alloc_err)
  if (alloc_err /= 0) stop 'ERROR: j_mat deallocation error'

  print*, 'Closing runtime files...'
  call close_files()   !close runtime files

  print*, 'DONE!'

end program couette_mod
