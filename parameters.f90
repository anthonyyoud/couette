MODULE parameters
implicit none
save

!****************************************************************************
!Data types by kind
!****************************************************************************
integer, parameter :: i1 = selected_int_kind(9)       !4-byte integer
integer, parameter :: r1 = selected_real_kind(6,37)   !4-byte real
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
real (r2),    parameter :: pi          = 3.14159265358979d0
real (r2),    parameter :: alpha       = 0d0 !3.13d0
real (r2),    parameter :: gamma       = 0.46d0 !(2d0 * pi) / alpha
real (r2),    parameter :: eta         = 0.5d0
real (r2),    parameter :: Q           = 0d0
real (r2),    parameter :: Re1         = 0d0
real (r2),    parameter :: Re2         = 0d0 !-1d0*(1d0/eta)*Re1
real (r2),    parameter :: Re1_mod     = 700d0
real (r2),    parameter :: Re2_mod     = 0d0
real (r2),    parameter :: om1         = 4d0
real (r2),    parameter :: om2         = 0d0
real (r2),    parameter :: dt          = 0.0001d0
real (r2),    parameter :: seed        = 1d-1
real (r2),    parameter :: end_time    = 3.154d0
real (r2),    parameter :: tau_init    = 1d0
real (r2),    parameter :: tau_step    = 1d0
real (r2),    parameter :: tau_end     = 1d0
integer (i1), parameter :: nx          = 80
integer (i1), parameter :: nt          = 20
integer (i1), parameter :: nz          = 40
integer (i1), parameter :: save_rate   = 10
integer (i1), parameter :: save_rate_2 = 39
logical,      parameter :: xsect_save  = .false.
logical,      parameter :: save3d      = .false.
logical,      parameter :: iso_hel     = .false.
logical,      parameter :: restart     = .false.
logical,      parameter :: auto_tau    = .false.
logical,      parameter :: save_part   = .false.
!****************************************************************************

!****************************************************************************
!Rarely used parameters
!****************************************************************************
integer (i1), parameter :: nb = (nx+1)*(nz+1)

real (r2),    parameter :: eps1 = 0d0 !0.30529d0
real (r2),    parameter :: freq1 = 0d0 !2d0 * pi / 2.8d0
real (r2),    parameter :: eps2 = 0d0
real (r2),    parameter :: freq2 = 0d0
real (r2),    parameter :: x_par_pos = 0.5d0
real (r2),    parameter :: z_par_pos = 0.3d0
!****************************************************************************

!****************************************************************************
!NOTHING BELOW HERE SHOULD BE CHANGED
!****************************************************************************
integer (i1), parameter :: Ntot = end_time / dt
real (r2),    parameter :: delx = 1d0 / nx
real (r2),    parameter :: delt = 2d0*pi / nt
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
