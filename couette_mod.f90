PROGRAM couette_mod
  !Main program file.  See README for a description of the code and how to run
  !it.
  USE parameters
  USE stream
  USE matrices
  USE io
  USE ic_bc
  USE variables
  USE linear
  USE nonlinear
  USE solve
  USE magnetic
  USE current
  IMPLICIT NONE
  
  INTEGER (i1)                      :: j, k, p=0, p_start=0, alloc_err=0, &
                                       desc_p(7), desc_b(7), desc_j(7)
  INTEGER (i1),         EXTERNAL    :: NUMROC
  REAL    (r2)                      :: p_fill(laf), b_fill(b_laf), &
                                       j_fill(laf), wtime(10), t=0.0_r2
  REAL    (r2),         ALLOCATABLE :: p_mat(:,:), b_mat(:,:), j_mat(:,:)
  REAL    (r2),         EXTERNAL    :: SLINQUIRE
  TYPE    (MAT_COMP)                :: Ux, Zx
  TYPE    (UZ_MAT_COMP)             :: Uz
  TYPE    (ZZ_MAT_COMP)             :: Zz
  LOGICAL                           :: state_exist
 
  IF (mycol == 0) THEN
    PRINT*
    PRINT*, 'Initialising ScaLAPACK process grid...'
  END IF

  CALL SLBOOT()   !initialise ScaLAPACK timer
  CALL SL_INIT(ictxt, nprow, npcol)   !initialise process grid
  CALL BLACS_GRIDINFO(ictxt, nprow, npcol, myrow, mycol)  !get process position
                                                          !in process grid
  IF ((myrow < 0) .OR. (mycol < 0)) STOP   !if process not in grid then exit

  IF (mycol == 0) THEN
    IF (ABS(tau - 0.0_r2) < EPSILON(tau)) THEN
      WRITE(6, '(A7, f4.2, A20)') 'tau = ', tau, '- Infinite cylinder'
    ELSE IF (ABS(tau - 1.0_r2) < EPSILON(tau)) THEN
      WRITE(6, '(A7, f4.2, A22)') 'tau = ', tau, '- Finite aspect ratio'
    ELSE IF ((ABS(tau - 0.0_r2) > EPSILON(tau)) .AND. &
             (ABS(tau - 1.0_r2) > EPSILON(tau))) THEN
      WRITE(6, '(A7, f4.2, A15)') 'tau = ', tau, '- Variable tau'
    ELSE
      WRITE(6, '(A7, f4.2)'), 'tau = ', tau
      STOP 'ERROR: Invalid value of tau'
    END IF   

    PRINT*, 'Opening runtime files...'

    CALL open_files()   !open files for writing

    PRINT*, 'Setting up ICS...'
  END IF

  !CALL BLACS_BARRIER(ictxt, 'A')  !wait until all processes get here

  IF (mycol == 0) PRINT*, 'Setting up finite-difference mesh...'

  CALL get_xzs()   !get finite-difference mesh

  IF (mycol == 0) THEN                                           !initial
    CALL ICS(ut%new, zt%new, psi%new, bt%new, jt%new, p_start)   !conditions

    IF (.NOT. restart) THEN
      INQUIRE(FILE='end_state.dat', EXIST=state_exist)
      !exit if not doing restart but end_state.dat exists
      IF (state_exist) STOP 'ERROR: restart=.FALSE. but end_state.dat exists.'
      PRINT*, 'Setting up BCS...'
      CALL u_BCS(ut%new, 0.0_r2)
      CALL p_BCS(psi%new)
      CALL z_BCS(zt%new, psi%new, 0.0_r2)   !get boundary conditions
      CALL b_BCS(bt%new)
      CALL j_BCS(jt%new)
    END IF
    PRINT*, 'Allocating matrix dimensions...'
  END IF

  !CALL BLACS_BARRIER(ictxt, 'A')  !wait until all processes get here

  CALL SLTIMER(1)
  IF (ABS(tau - 0.0_r2) < EPSILON(tau)) THEN
    b_M = NUMROC(2*nxp1+1, 2*nxp1+1, myrow, 0, nprow)
    b_N = NUMROC(nxp1*nz1, nb, mycol, 0, npcol)
    j_M = NUMROC(2*nx1+1, 2*nx1+1, myrow, 0, nprow) !allocate
    j_N = NUMROC(nx1*nzp1, nb, mycol, 0, npcol)     !matrix dimensions local
    ALLOCATE(b_mat(b_M,b_N), STAT=alloc_err)        !to each process
    IF (alloc_err /= 0) STOP 'ERROR: b_mat allocation error'
    ALLOCATE(j_mat(j_M,j_N), STAT=alloc_err)        !depending on finite
    IF (alloc_err /= 0) STOP 'ERROR: j_mat allocation error'
  ELSE IF (ABS(tau - 1.0_r2) < EPSILON(tau)) THEN   !or infinite cylinders
    b_M = NUMROC(2*nxp1+1, 2*nxp1+1, myrow, 0, nprow)
    b_N = NUMROC(nxp1*nzp1, nb, mycol, 0, npcol)
    j_M = NUMROC(2*nx1+1, 2*nx1+1, myrow, 0, nprow)
    j_N = NUMROC(nx1*nz1, nb, mycol, 0, npcol)
    ALLOCATE(b_mat(b_M,b_N), STAT=alloc_err) 
    IF (alloc_err /= 0) STOP 'ERROR: b_mat allocation error'
    ALLOCATE(j_mat(j_M,j_N), STAT=alloc_err)
    IF (alloc_err /= 0) STOP 'ERROR: j_mat allocation error'
  ELSE IF ((ABS(tau - 0.0_r2) > EPSILON(tau)) .AND. &
           (ABS(tau - 1.0_r2) > EPSILON(tau))) THEN
    b_M = NUMROC(2*nxp1+1, 2*nxp1+1, myrow, 0, nprow)
    b_N = NUMROC(nxp1*nzp1, nb, mycol, 0, npcol)
    j_M = NUMROC(2*nx1+1, 2*nx1+1, myrow, 0, nprow)
    j_N = NUMROC(nx1*nzp1, nb, mycol, 0, npcol)
    ALLOCATE(b_mat(b_M,b_N), STAT=alloc_err)
    IF (alloc_err /= 0) STOP 'ERROR: b_mat allocation error'
    ALLOCATE(j_mat(j_M,j_N), STAT=alloc_err)
    IF (alloc_err /= 0) STOP 'ERROR: j_mat allocation error'
  ELSE
    STOP 'ERROR: Allocation error due to value of tau'
  END IF

  p_M = NUMROC(2*nx1+1, 2*nx1+1, myrow, 0, nprow)
  p_N = NUMROC(nx1*nz1, nb, mycol, 0, npcol)
  ALLOCATE(p_mat(p_M,p_N), STAT=alloc_err)
    IF (alloc_err /= 0) STOP 'ERROR: p_mat allocation error'
  CALL SLTIMER(1)

  IF (mycol == 0) THEN
    PRINT*, 'Setting up matrices...'   !matrices for tridiagonal systems
    CALL matrix_setup(Ux, Uz, Zx, Zz)  !azimuthal velocity, vorticity
  END IF

  CALL SLTIMER(2)
  CALL psi_mat_setup(p_mat, desc_p, p_fill)

  IF (ABS(tau - 0.0_r2) < EPSILON(tau)) THEN
    CALL b_mat_setup(b_mat, desc_b, b_fill)
    CALL j_mat_setup(j_mat, desc_j, j_fill)
  ELSE IF (ABS(tau - 1.0_r2) < EPSILON(tau)) THEN     !left-hand side matrices
    CALL fin_b_mat_setup(b_mat, desc_b, b_fill)       !for stream-function,
    CALL fin_j_mat_setup(j_mat, desc_j, j_fill)       !current and magnetic
  ELSE IF ((ABS(tau - 0.0_r2) > EPSILON(tau)) .AND. & !Poisson equations
           (ABS(tau - 1.0_r2) > EPSILON(tau))) THEN
    CALL fin_b_mat_setup(b_mat, desc_b, b_fill)
    CALL j_mat_setup(j_mat, desc_j, j_fill)
  ELSE
    STOP 'ERROR: Matrix setup subroutine call error due to value of tau'
  END IF
  CALL SLTIMER(2)

  IF (mycol == 0) THEN
    CALL copy_var(ut%old, ut%new)
    CALL copy_var(ut%old2, ut%new)
    CALL copy_var(ut%inter, ut%new)
    CALL copy_var(zt%old, zt%new)
    CALL copy_var(zt%old2, zt%new)
    CALL copy_var(zt%inter, zt%new)   !initialise all variables equal to
    CALL copy_var(psi%old, psi%new)   !'new' from ICS and BCS
    CALL copy_var(psi%old2, psi%new)
    CALL copy_var(bt%old, bt%new)
    CALL copy_var(bt%old, bt%new)
    CALL copy_var(jt%old, jt%new)
    CALL copy_var(jt%old, jt%new)

    PRINT*, 'Entering time loop...'
  END IF

  DO p = p_start, Ntot        !start main time loop
  ! CALL BLACS_BARRIER(ictxt, 'A')  !wait until all processes get here
    CALL terminate(p, t)     !loop terminates if file 'RUNNING'
                              !does not exist in run directory
    IF (mycol == 0) THEN
      IF (end_proc == 1) EXIT   !end master process
    END IF

    IF (mycol > 0) THEN
      IF (end_proc == 1) THEN
        WRITE(6, '(A16, i2)'), 'Ending process ', mycol !end other processes
        EXIT
      END IF
    END IF

    IF (mycol == 0) THEN     !save cross-sections, surfaces if
      CALL save_run(p, t)      !file 'SAVE' exists in run directory
    END IF

    t = p * dt   !increment time

    IF (mycol == 0) THEN
      IF (MOD(p, save_rate) == 0) CALL write_data(p, p_start, t)
                                                    !do necessary writing
      CALL copy_var(vrold, vr)                      !of data to files
      CALL copy_var(vzold, vz)   !update variables

      IF (xsect_save) THEN   !save cross-sections
        IF (MOD(p, save_rate_2) == 0) THEN
          CALL save_xsect(vr, vz, psi%old, ut%new, zt%new, &
                          bt%old, jt%old, t, p)
          !CALL save_surface(psi%old, ut%new, zt%new, vr, vz, &
          !                  bt%old, jt%old, p, t)
        END IF
      END IF

      IF (save3d) THEN   !save 3D surface
        IF (MOD(p, save_rate_2) == 0) THEN
          CALL save_3d(vr, ut%new, vz, psi%old, p)
        END IF
      END IF

      CALL copy_var(ut%old, ut%new)
      CALL copy_var(zt%old, zt%new)   !update variables

      !get linear part for u in x-direction
      CALL get_rhs_ux(ut%old, ut%new)
      !get nonlinear part for u in x-direction
      CALL get_nlin_ux(ut%old, ut%old2, psi%old, psi%old2, &
                       bt%old, bt%old2, ut%nlin_new)
      !get linear part for Z in x-direction
      CALL get_rhs_Zx(zt%old, zt%new)
      !get nonlinear part for Z in x-direction
      CALL get_nlin_Zx(t, ut%old, ut%old2, psi%old, psi%old2, zt%old, &
                       zt%old2, jt%old, jt%old2, zt%nlin_new)
      !solve for u in x-direction
      CALL solve_ux(ut%old, ut%new, ut%nlin_new, t, Ux)
      !solve for Z in x-direction
      CALL solve_Zx(zt%old, zt%new, zt%nlin_new, psi%old, t, Zx)

      CALL copy_var(ut%old2, ut%inter)
      CALL copy_var(zt%old2, zt%inter)   !update variables

      !solve for u in z-direction
      CALL solve_uz(ut%old, ut%new, t, Uz)
      !solve for Z in z-direction
      CALL solve_Zz(zt%old, psi%old, zt%new, t, Zz)

      CALL copy_var(ut%inter, ut%new)
      CALL copy_var(zt%inter, zt%new)
      CALL copy_var(psi%old2, psi%old)   !update variables
      CALL copy_var(bt%old2, bt%old)
      CALL copy_var(jt%old2, jt%old)

      CALL p_BCS(psi%new)
      CALL b_BCS(bt%new)   !update boundary conditions
      CALL j_BCS(jt%new)
    END IF
   
    CALL BLACS_BARRIER(ictxt, 'A')  !wait until all processes get here

    CALL SLTIMER(4)
    !broadcast u and Z to other processes for use as RHS in Poisson equations
    IF (npcol > 1) THEN
      !PRINT*,mycol,'start'
      CALL DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, zt%new, nxp1, 0, 0)
      CALL DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, ut%new, nxp1, 0, 0)
      !PRINT*,mycol,'end'
    END IF
    CALL SLTIMER(4)

    CALL BLACS_BARRIER(ictxt, 'A')  !wait until all processes get here

    CALL SLTIMER(5)
    !solve Poisson equations depending on tau
    CALL p_poisson(zt%new, psi%new, p_mat, desc_p, p_fill)
    IF (ABS(tau - 0.0_r2) < EPSILON(tau)) THEN
      CALL b_poisson(ut%new, bt%new, b_mat, desc_b, b_fill)
      CALL j_poisson(psi%new, jt%new, j_mat, desc_j, j_fill)
    ELSE IF (ABS(tau - 1.0_r2) < EPSILON(tau)) THEN
      CALL fin_b_poisson(ut%new, bt%new, b_mat, desc_b, b_fill)
      CALL fin_j_poisson(psi%new, jt%new, j_mat, desc_j, j_fill)
    ELSE IF ((ABS(tau - 0.0_r2) > EPSILON(tau)) .AND. &
             (ABS(tau - 1.0_r2) > EPSILON(tau))) THEN
      CALL fin_b_poisson(ut%new, bt%new, b_mat, desc_b, b_fill)
      CALL j_poisson(psi%new, jt%new, j_mat, desc_j, j_fill)
    ELSE
      STOP 'ERROR: Poisson solver subroutine call error due to value of tau'
    END IF
    CALL SLTIMER(5)

    IF (mycol == 0) THEN
      CALL copy_var(psi%old, psi%new)
      CALL copy_var(bt%old, bt%new)   !update variables
      CALL copy_var(jt%old, jt%new)
      IF (p == Ntot) THEN
        CALL end_state(ut%new, zt%new, psi%old, bt%old, jt%old, p)
        CALL save_xsect(vr, vz, psi%old, ut%new, zt%new, &  !save final state
                        bt%old, jt%old, t, p)
        CALL save_surface(psi%old, ut%new, zt%new, &
                          vr, vz, bt%old, jt%old, p, t)
      END IF
    END IF
  END DO   !end time loop

  IF (mycol == 0) PRINT*, 'Deallocating allocated arrays...'

  IF (ALLOCATED(p_mat)) DEALLOCATE(p_mat, STAT=alloc_err)
  IF (alloc_err /= 0) STOP 'ERROR: p_mat deallocation error'
  IF (ALLOCATED(b_mat)) DEALLOCATE(b_mat, STAT=alloc_err) !deallocate arrays
  IF (alloc_err /= 0) STOP 'ERROR: b_mat deallocation error'
  IF (ALLOCATED(j_mat)) DEALLOCATE(j_mat, STAT=alloc_err)
  IF (alloc_err /= 0) STOP 'ERROR: j_mat deallocation error'

  IF (mycol == 0) THEN
    PRINT*, 'Closing runtime files...'
    CALL close_files()   !close runtime files
  END IF

  IF (mycol == 0) THEN
    PRINT*, '*****************************************************************'
    PRINT*, 'ScaLAPACK timer information'
    WRITE(6,*) 'proc setup end_proc broadcast solve p_sum j_sum b_sum j_copy'
  END IF
  WRITE(6,'(i2,8f9.2)') mycol, SLINQUIRE('W', 2), &
                        SLINQUIRE('W', 3), SLINQUIRE('W', 4), &
                        SLINQUIRE('W', 5), SLINQUIRE('W', 6), &
                        SLINQUIRE('W', 7), SLINQUIRE('W', 8), &
                        SLINQUIRE('W', 9)
  IF (mycol == 0) THEN
    PRINT*, '*****************************************************************'
  END IF

  IF (mycol == 0) PRINT*, 'Releasing ScaLAPACK process grid...'

  CALL BLACS_BARRIER(ictxt, 'A')  !wait until all processes get here
  CALL BLACS_GRIDEXIT(ictxt)   !release current process grid
  CALL BLACS_EXIT(0)   !release all process grids

  PRINT*, 'DONE!'

END PROGRAM couette_mod
