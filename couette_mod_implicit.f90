PROGRAM couette_mod
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

type (mat_comp) :: Ux, Zx
type (uz_mat_comp) :: Uz
type (zz_mat_comp) :: Zz
double precision :: p_fill(laf), b_fill(b_laf), j_fill(laf), &
                    wtime(10), t = 0d0
double precision, allocatable :: p_mat(:,:), b_mat(:,:), j_mat(:,:)
double precision, external :: SLINQUIRE
integer, external :: NUMROC
integer :: j, k, p = 0, p_start = 0, desc_p(7), desc_b(7), desc_j(7)
logical :: state_exist

call SLBOOT()   !initialise ScaLAPACK timer
call SL_INIT(ictxt, nprow, npcol)   !initialise process grid
call BLACS_GRIDINFO(ictxt, nprow, npcol, myrow, mycol)  !get process position
                                                        !in process grid
if ((myrow < 0) .or. (mycol < 0)) stop   !if process not in grid then exit

if (mycol == 0) then
   print*
   if (tau == 0) then
      write(6, '(A7, f4.2, A21)') 'tau = ', tau, '- Infinite cylinder'
   else if (tau == 1) then
      write(6, '(A7, f4.2, A22)') 'tau = ', tau, '- Finite aspect ratio'
   else
      write(6, '(A7, f4.2)') 'tau = ', tau
   end if   

   call open_files()   !open files for writing

   print*, 'Setting up ICS...'
end if

!call BLACS_BARRIER(ictxt, 'A')  !wait until all processes get here

call get_xzs()   !get finite-difference mesh

if (mycol == 0) then                                            !initial
   call ICS(ut%new, zt%new, psi%new, bt%new, jt%new, p_start)   !conditions

   if (.not. restart) then
      inquire(file='end_state.dat', exist=state_exist)
      !exit if not doing restart but end_state.dat exists
      if (state_exist) STOP 'restart=.false. but end_state.dat exists.'
      print*, 'Setting up BCS...'
      call u_BCS(ut%new, 0d0)
      call p_BCS(psi%new)
      call z_BCS(zt%new, psi%new, 0d0)   !get boundary conditions
      call b_BCS(bt%new)
      call j_BCS(jt%new)
   end if
   print*, 'Allocating matrix dimensions'
end if

!call BLACS_BARRIER(ictxt, 'A')  !wait until all processes get here

call SLTIMER(1)
if (tau == 0d0) then
   b_M = NUMROC(2*nxp1+1, 2*nxp1+1, myrow, 0, nprow)
   b_N = NUMROC(nxp1*nz1, nb, mycol, 0, npcol)
   j_M = NUMROC(2*nx1+1, 2*nx1+1, myrow, 0, nprow)      !allocate the matrix
   j_N = NUMROC(nx1*nzp1, nb, mycol, 0, npcol)          !dimensions local
   allocate(b_mat(b_M,b_N))                             !to each process
   allocate(j_mat(j_M,j_N))                             !depending on finite
else if (tau == 1d0) then                               !or infinite
   b_M = NUMROC(2*nxp1+1, 2*nxp1+1, myrow, 0, nprow)    !cylinders
   b_N = NUMROC(nxp1*nzp1, nb, mycol, 0, npcol)
   j_M = NUMROC(2*nx1+1, 2*nx1+1, myrow, 0, nprow)
   j_N = NUMROC(nx1*nz1, nb, mycol, 0, npcol)
   allocate(b_mat(b_M,b_N)) 
   allocate(j_mat(j_M,j_N))
else if ((tau /= 0d0) .and. (tau /= 1d0)) then
   b_M = NUMROC(2*nxp1+1, 2*nxp1+1, myrow, 0, nprow)
   b_N = NUMROC(nxp1*nzp1, nb, mycol, 0, npcol)
   j_M = NUMROC(2*nx1+1, 2*nx1+1, myrow, 0, nprow)
   j_N = NUMROC(nx1*nzp1, nb, mycol, 0, npcol)
   allocate(b_mat(b_M,b_N))
   allocate(j_mat(j_M,j_N))
end if

p_M = NUMROC(2*nx1+1, 2*nx1+1, myrow, 0, nprow)
p_N = NUMROC(nx1*nz1, nb, mycol, 0, npcol)
allocate(p_mat(p_M,p_N))
call SLTIMER(1)

if (mycol == 0) then
   print*, 'Setting up matrices...'   !matrices for tridiagonal systems
   call matrix_setup(Ux, Uz, Zx, Zz)  !azimuthal velocity, vorticity
end if

call SLTIMER(2)
call psi_mat_setup(p_mat, desc_p, p_fill)

if (tau == 0d0) then
   call b_mat_setup(b_mat, desc_b, b_fill)
   call j_mat_setup(j_mat, desc_j, j_fill)
else if (tau == 1d0) then                         !left-hand side matrices
   call fin_b_mat_setup(b_mat, desc_b, b_fill)    !for stream-function,
   call fin_j_mat_setup(j_mat, desc_j, j_fill)    !current and magnetic
else if ((tau /= 0d0) .and. (tau /= 1d0)) then    !Poisson equations
   call fin_b_mat_setup(b_mat, desc_b, b_fill)
   call j_mat_setup(j_mat, desc_j, j_fill)
end if
call SLTIMER(2)

if (mycol == 0) then
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

   print*, 'Entering time loop'
end if

do p = p_start, Ntot        !start main time loop
!   call BLACS_BARRIER(ictxt, 'A')  !wait until all processes get here
   call terminate(p, t)     !loop terminates if file 'RUNNING'
                            !does not exist in run directory
   if (mycol == 0) then
      if (end_proc == 1) exit   !end master process
   end if

   if (mycol > 0) then
      if (end_proc == 1) then
         print*, 'Ending process ', mycol   !end additional processes
         exit
      end if
   end if

   if (mycol == 0) then     !save cross-sections, surfaces if
      call save_run(p, t)   !file 'SAVE' exists in run directory
   end if

   t = p * dt   !increment time

   if (mycol == 0) then
      if (mod(p, save_rate) == 0) call write_data(p, p_start, t)
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
   end if
   
   call BLACS_BARRIER(ictxt, 'A')  !wait until all processes get here

   call SLTIMER(4)
!broadcast u and Z to other processes for use as RHS in Poisson equations
   if (npcol > 1) then
!print*,mycol,'start'
      call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, zt%new, nxp1, 0, 0)
      call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, ut%new, nxp1, 0, 0)
!print*,mycol,'end'
   end if
   call SLTIMER(4)

   call BLACS_BARRIER(ictxt, 'A')  !wait until all processes get here

   call SLTIMER(5)
!solve Poisson equations depending on tau
   call p_poisson(zt%new, psi%new, p_mat, desc_p, p_fill)
   if (tau == 0d0) then
      call b_poisson(ut%new, bt%new, b_mat, desc_b, b_fill)
      call j_poisson(psi%new, jt%new, j_mat, desc_j, j_fill)
   else if (tau == 1d0) then
      call fin_b_poisson(ut%new, bt%new, b_mat, desc_b, b_fill)
      call fin_j_poisson(psi%new, jt%new, j_mat, desc_j, j_fill)
   else if ((tau /= 0d0) .and. (tau /= 1d0)) then
      call fin_b_poisson(ut%new, bt%new, b_mat, desc_b, b_fill)
      call j_poisson(psi%new, jt%new, j_mat, desc_j, j_fill)
   end if
   call SLTIMER(5)

   if (mycol == 0) then
      call copy_var(psi%old, psi%new)
      call copy_var(bt%old, bt%new)   !update variables
      call copy_var(jt%old, jt%new)
      if (p == Ntot) then
         call end_state(ut%new, zt%new, psi%old, bt%old, jt%old, p)
         call save_xsect(vr, vz, psi%old, ut%new, zt%new, &  !save final state
                         bt%old, jt%old, t, p)
         call save_surface(psi%old, ut%new, zt%new, &
                           vr, vz, bt%old, jt%old, p, t)
      end if
   end if
end do   !end time loop

deallocate(p_mat, b_mat, j_mat)   !deallocate allocated arrays

if (mycol == 0) then
   call close_files()   !close runtime files
end if

if (mycol == 0) then
   write(6,*) 'proc setup end_proc broadcast solve p_sum j_sum b_sum j_copy'
end if
write(6,'(i2,8f9.2)') mycol, SLINQUIRE('W', 2), &
                      SLINQUIRE('W', 3), SLINQUIRE('W', 4), &
                      SLINQUIRE('W', 5), SLINQUIRE('W', 6), &
                      SLINQUIRE('W', 7), SLINQUIRE('W', 8), &
                      SLINQUIRE('W', 9)

call BLACS_BARRIER(ictxt, 'A')  !wait until all processes get here
call BLACS_GRIDEXIT(ictxt)   !release current process grid
call BLACS_EXIT(0)   !release all process grids

END PROGRAM couette_mod
