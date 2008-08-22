module parameters
  !Parameters to set.
  implicit none
  save

  !****************************************************************************
  !Data types by kind
  !****************************************************************************
  integer,      parameter :: i1 = selected_int_kind(9)       !4-byte integer
  !integer (i1), parameter :: r2 = selected_real_kind(6,37)   !4-byte real
  integer (i1), parameter :: r2 = selected_real_kind(15,307) !8-byte real
  !****************************************************************************

  !****************************************************************************
  !Parameters to set
  !****************************************************************************
  real    (r2), parameter :: pi          = 3.14159265358979_r2
  real    (r2), parameter :: alpha       = 0.0_r2 !3.73_r2
  real    (r2), parameter :: gamma       = 3.0_r2 !(2.0_r2 * pi) / alpha
  real    (r2), parameter :: eta         = 0.615_r2
  real    (r2), parameter :: Q           = 0.0_r2
  real    (r2)            :: Re1         = 500.0_r2
  real    (r2), parameter :: Re_incr     = 1.0_r2
  real    (r2), parameter :: growth_tol  = 1E-8_r2
  real    (r2), parameter :: Re2         = 0.0_r2 !-1.0_r2*(1.0_r2/eta)*Re1
  real    (r2), parameter :: Re1_mod     = 0.0_r2
  real    (r2), parameter :: Re2_mod     = 0.0_r2
  real    (r2), parameter :: om1         = 0.0_r2
  real    (r2), parameter :: om2         = 0.0_r2
  real    (r2), parameter :: dt          = 0.0001_r2
  real    (r2), parameter :: seed        = 1E-3_r2
  real    (r2), parameter :: end_time    = 20.0_r2
  real    (r2), parameter :: tau_init    = 1.0_r2
  real    (r2), parameter :: tau_step    = 1.0_r2
  real    (r2), parameter :: tau_end     = 1.0_r2
  integer (i1), parameter :: nx          = 40
  integer (i1), parameter :: nt          = 20
  integer (i1), parameter :: nz          = 120
  integer (i1), parameter :: save_rate   = 10
  integer (i1), parameter :: save_rate_2 = 104
  integer (i1), parameter :: periodic_save = 10
  logical,      parameter :: rot_ends    = .false.
  logical,      parameter :: xsect_save  = .false.
  logical,      parameter :: save3d      = .false.
  logical,      parameter :: iso_hel     = .false.
  logical,      parameter :: restart     = .true.
  logical,      parameter :: auto_tau    = .false.
  logical,      parameter :: auto_Re     = .false.
  logical,      parameter :: dec_Re      = .false.
  logical,      parameter :: hyst_Re     = .false.
  logical,      parameter :: divergence  = .false.
  !****************************************************************************

  !****************************************************************************
  !Rarely used parameters
  !****************************************************************************
  real    (r2), parameter :: eps1      = 0.0_r2 !0.30529_r2
  real    (r2), parameter :: freq1     = 0.0_r2 !2.0_r2 * pi / 2.8_r2
  real    (r2), parameter :: eps2      = 0.0_r2
  real    (r2), parameter :: freq2     = 0.0_r2
  integer (i1), parameter :: num_pars  = 50_i1
  logical,      parameter :: save_part = .false.
  !****************************************************************************

  !****************************************************************************
  !NOTHING BELOW HERE SHOULD BE CHANGED
  !****************************************************************************
  integer (i1), parameter :: Ntot = end_time / dt
  real    (r2), parameter :: delx = 1.0_r2 / nx
  real    (r2), parameter :: delt = 2.0_r2*pi / nt
  real    (r2), parameter :: delz = gamma / nz
  real    (r2), parameter :: dx2 = delx ** 2
  real    (r2), parameter :: dz2 = delz ** 2
  real    (r2), parameter :: rx = dt / delx
  real    (r2), parameter :: rz = dt / delz
  real    (r2), parameter :: rxx = dt / dx2
  real    (r2), parameter :: rzz = dt / dz2
  real    (r2), parameter :: one_eta = 1.0_r2 - eta
  integer (i1), parameter :: nx1 = nx-1
  integer (i1), parameter :: nz1 = nz-1
  integer (i1), parameter :: nxp1 = nx+1
  integer (i1), parameter :: nzp1 = nz+1
  integer (i1), parameter :: xlb = 1
  integer (i1), parameter :: zlb = 0
  real    (r2)            :: tau = tau_init
  logical (r2)            :: init_Re = .true.
  logical (r2)            :: zero_Re = .false.
  logical (r2)            :: gm_set = .false.
  logical (r2)            :: gp_set = .false.
 
  integer (i1) :: end_proc = 0
  integer (i1) :: saturated = 0
  integer (i1), parameter :: SPr = selected_real_kind(6,37) !4-byte real
  integer (i1), parameter :: DPr = selected_real_kind(15,307) !8-byte real
  !****************************************************************************

end module parameters
