MODULE parameters
implicit none
save

!****************************************************************************
!Data types by kind
!****************************************************************************
integer, parameter :: i1 = selected_int_kind(9)       !4-byte integer
!integer, parameter :: r2 = selected_real_kind(6,37)   !4-byte real
integer, parameter :: r2 = selected_real_kind(15,307) !8-byte real
!****************************************************************************

!****************************************************************************
!Process rows and columns for ScaLAPACK
!****************************************************************************
integer (i1), parameter :: nprow       = 1
integer (i1), parameter :: npcol       = 1
!integer (i1), parameter :: nb          = 1681
!****************************************************************************

!****************************************************************************
!Parameters to set
!****************************************************************************
real (r2),    parameter :: pi          = 3.14159265358979_r2
real (r2),    parameter :: alpha       = 0.0_r2 !3.13_r2
real (r2),    parameter :: gamma       = 0.46_r2 !(2.0_r2 * pi) / alpha
real (r2),    parameter :: eta         = 0.5_r2
real (r2),    parameter :: Q           = 0.0_r2
real (r2),    parameter :: Re1         = 0.0_r2
real (r2),    parameter :: Re2         = 0.0_r2 !-1.0_r2*(1.0_r2/eta)*Re1
real (r2),    parameter :: Re1_mod     = 700.0_r2
real (r2),    parameter :: Re2_mod     = 0.0_r2
real (r2),    parameter :: om1         = 4.0_r2
real (r2),    parameter :: om2         = 0.0_r2
real (r2),    parameter :: dt          = 0.0001_r2
real (r2),    parameter :: seed        = 1e-1_r2
real (r2),    parameter :: end_time    = 3.94_r2
real (r2),    parameter :: tau_init    = 1.0_r2
real (r2),    parameter :: tau_step    = 1.0_r2
real (r2),    parameter :: tau_end     = 1.0_r2
integer (i1), parameter :: nx          = 80
integer (i1), parameter :: nt          = 20
integer (i1), parameter :: nz          = 40
integer (i1), parameter :: save_rate   = 10
integer (i1), parameter :: save_rate_2 = 78
logical,      parameter :: xsect_save  = .true.
logical,      parameter :: save3d      = .false.
logical,      parameter :: iso_hel     = .false.
logical,      parameter :: restart     = .true.
logical,      parameter :: auto_tau    = .false.
logical,      parameter :: save_part   = .false.
!****************************************************************************

!****************************************************************************
!Rarely used parameters
!****************************************************************************
integer (i1), parameter :: nb = (nx+1)*(nz+1)

real (r2),    parameter :: eps1 = 0.0_r2 !0.30529_r2
real (r2),    parameter :: freq1 = 0.0_r2 !2.0_r2 * pi / 2.8_r2
real (r2),    parameter :: eps2 = 0.0_r2
real (r2),    parameter :: freq2 = 0.0_r2
real (r2),    parameter :: x_par_pos = 0.5_r2
real (r2),    parameter :: z_par_pos = 0.3_r2
!****************************************************************************

!****************************************************************************
!NOTHING BELOW HERE SHOULD BE CHANGED
!****************************************************************************
integer (i1), parameter :: Ntot = end_time / dt
real (r2),    parameter :: delx = 1.0_r2 / nx
real (r2),    parameter :: delt = 2.0_r2*pi / nt
real (r2),    parameter :: delz = gamma / nz
real (r2),    parameter :: dx2 = delx ** 2
real (r2),    parameter :: dz2 = delz ** 2
real (r2),    parameter :: rx = dt / delx
real (r2),    parameter :: rz = dt / delz
real (r2),    parameter :: rxx = dt / dx2
real (r2),    parameter :: rzz = dt / dz2
integer (i1), parameter :: nx1 = nx-1
integer (i1), parameter :: nz1 = nz-1
integer (i1), parameter :: nxp1 = nx+1
integer (i1), parameter :: nzp1 = nz+1
integer (i1), parameter :: xlb = 1
integer (i1), parameter :: zlb = 0
real (r2) :: tau = tau_init
real (r2) :: x_pos = x_par_pos * nx
real (r2) :: z_pos = z_par_pos * nz

integer (i1), parameter :: laf   = nb*(nx1+nx1)+6*nx1*nx1
integer (i1), parameter :: b_laf = nb*(nxp1+nxp1)+6*nxp1*nxp1
integer (i1), parameter :: lwork_fac   = nx1*nx1
integer (i1), parameter :: lwork_b_fac = nxp1*nxp1
integer (i1), parameter :: lwork_sol   = nx1
integer (i1), parameter :: lwork_b_sol = nxp1
integer (i1) :: ictxt
integer (i1) :: myrow
integer (i1) :: mycol
integer (i1) :: b_M
integer (i1) :: b_N
integer (i1) :: j_M
integer (i1) :: j_N
integer (i1) :: p_M
integer (i1) :: p_N
integer (i1) :: end_proc = 0
!****************************************************************************

END MODULE parameters
