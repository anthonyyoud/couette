PROGRAM couette_mod
use parameters
use pressure
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
double precision :: t = 0d0, p_fill(laf), b_fill(b_laf), j_fill(laf), &
                    wtime(10)
double precision, allocatable :: p_mat(:,:), b_mat(:,:), j_mat(:,:)
double precision, external :: SLINQUIRE
integer, external :: NUMROC
integer :: j, k, p = 0, p_start = 0, &
           desc_p(7), desc_b(7), desc_j(7)
logical :: state_exist

call SLBOOT()
call SL_INIT(ictxt, nprow, npcol)
call BLACS_GRIDINFO(ictxt, nprow, npcol, myrow, mycol)
if ((myrow < 0) .or. (mycol < 0)) stop

if (mycol == 0) then
   print*
   if (tau == 0) then
      write(6, '(A7, f4.2, A21)') 'tau = ', tau, '- Infinite cylinder'
   else if (tau == 1) then
      write(6, '(A7, f4.2, A22)') 'tau = ', tau, '- Finite aspect ratio'
   else
      write(6, '(A7, f4.2)') 'tau = ', tau
   end if   

   call open_files()

   print*, 'Setting up ICS...'
end if

call get_xzs()

if (mycol == 0) then
   call ICS(ut%new, zt%new, psi%new, bt%new, jt%new, p_start)

   if (.not. restart) then
      inquire(file='end_state.dat', exist=state_exist)
      if (state_exist) then
         STOP 'restart=.false. but end_state.dat exists.'
      end if
      print*, 'Setting up BCS...'
      call u_BCS(ut%new, 0d0)
      call p_BCS(psi%new)
      call z_BCS(zt%new, psi%new, 0d0)
      call b_BCS(bt%new)
      call j_BCS(jt%new)
   end if
   print*, 'Allocating matrix dimensions'
end if

call SLTIMER(1)
if (tau == 0d0) then
   b_M = NUMROC(2*nxp1+1, 2*nxp1+1, myrow, 0, nprow)
   b_N = NUMROC(nxp1*nz1, nb, mycol, 0, npcol)
   j_M = NUMROC(2*nx1+1, 2*nx1+1, myrow, 0, nprow)
   j_N = NUMROC(nx1*nzp1, nb, mycol, 0, npcol)
   allocate(b_mat(b_M,b_N))
   allocate(j_mat(j_M,j_N))
else if (tau == 1d0) then 
   b_M = NUMROC(2*nxp1+1, 2*nxp1+1, myrow, 0, nprow)
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
   print*, 'Setting up matrices...'
   call matrix_setup(Ux, Uz, Zx, Zz)
end if

call SLTIMER(2)
call psi_mat_setup(p_mat, desc_p, p_fill)

if (tau == 0d0) then
   call b_mat_setup(b_mat, desc_b, b_fill)
   call j_mat_setup(j_mat, desc_j, j_fill)
else if (tau == 1d0) then
   call fin_b_mat_setup(b_mat, desc_b, b_fill)
   call fin_j_mat_setup(j_mat, desc_j, j_fill)
else if ((tau /= 0d0) .and. (tau /= 1d0)) then
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
   call copy_var(zt%inter, zt%new)
   call copy_var(psi%old, psi%new)
   call copy_var(psi%old2, psi%new)
   call copy_var(bt%old, bt%new)
   call copy_var(bt%old, bt%new)
   call copy_var(jt%old, jt%new)
   call copy_var(jt%old, jt%new)

   print*, 'Entering time loop'
end if

do p = p_start, Ntot
   call terminate(p, t)

   if (mycol == 0) then
      if (end_proc == 1) exit
   end if

   if (mycol > 0) then
      if (end_proc == 1) then
         print*, 'Ending process ', mycol
         exit
      end if
   end if

   if (mycol == 0) then
      call save_run(p, t)
   end if

   t = p * dt

   if (mycol == 0) then
      if (mod(p, save_rate) == 0) call write_data(p, p_start, t)

      call copy_var(vrold, vr)
      call copy_var(vzold, vz)

      if (xsect_save) then
         if (mod(p, save_rate_2) == 0) then
            call save_xsect(vr, vz, psi%old, t, p)
            !call save_surface(psi%old, ut%new, zt%new, vr, vz, &
            !                  bt%old, jt%old, p, t)
         end if
      end if

      call copy_var(ut%old, ut%new)
      call copy_var(zt%old, zt%new)

      call get_rhs_ux(ut%old, ut%new)
      call get_nlin_ux(ut%old, ut%old2, psi%old, psi%old2, &
                       bt%old, bt%old2, ut%nlin_new)
      call get_rhs_Zx(zt%old, zt%new)
      call get_nlin_Zx(t, ut%old, ut%old2, psi%old, psi%old2, zt%old, &
                    zt%old2, jt%old, jt%old2, zt%nlin_new)
      call solve_ux(ut%old, ut%new, ut%nlin_new, t, Ux)
      call solve_Zx(zt%old, zt%new, zt%nlin_new, psi%old, t, Zx)

      call copy_var(ut%old2, ut%inter)
      call copy_var(zt%old2, zt%inter)

      call solve_uz(ut%old, ut%new, t, Uz)
      call solve_Zz(zt%old, psi%old, zt%new, t, Zz)

      call copy_var(ut%inter, ut%new)
      call copy_var(zt%inter, zt%new)
      call copy_var(psi%old2, psi%old)
      call copy_var(bt%old2, bt%old)
      call copy_var(jt%old2, jt%old)

      call p_BCS(psi%new)
      call b_BCS(bt%new)
      call j_BCS(jt%new)
   end if

   call SLTIMER(4)
   call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, zt%new, nxp1, 0, 0)
   call DGEBR2D(ictxt, 'A', ' ', nxp1, nzp1, ut%new, nxp1, 0, 0)
   call SLTIMER(4)

   call SLTIMER(5)
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

   if (diag) then
      call calc_rhs_u(ut%new, ut%old2, psi%old, t, p)
      call calc_rhs_Z(zt%new, zt%old2, ut%old2, psi%old, t, p)
      call calc_rhs_psi(psi%new, zt%new, t, p)
   end if

   if (mycol == 0) then
      call copy_var(psi%old, psi%new)
      call copy_var(bt%old, bt%new)
      call copy_var(jt%old, jt%new)
      if (p == Ntot) then
         call end_state(ut%new, zt%new, psi%old, bt%old, jt%old, p)
         call save_xsect(vr, vz, psi%old, t, p)
         call save_surface(psi%old, ut%new, zt%new, &
                           vr, vz, bt%old, jt%old, p, t)
      end if
   end if
end do

if (mycol == 0) then
   call close_files()
end if

if (mycol == 0) then
   write(6,*) 'proc setup end_proc matrices solve p_sum j_sum j_copy b_sum'
end if
write(6,'(i3,8f9.4)') mycol, SLINQUIRE('W', 2), &
                      SLINQUIRE('W', 3), SLINQUIRE('W', 4), &
                      SLINQUIRE('W', 5), SLINQUIRE('W', 6), &
                      SLINQUIRE('W', 7), SLINQUIRE('W', 9), &
                      SLINQUIRE('W', 8)

call BLACS_BARRIER(ictxt, 'A')
call BLACS_GRIDEXIT(ictxt)
call BLACS_EXIT(0)

END PROGRAM couette_mod
