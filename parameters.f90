module parameters
  !Parameters to set.
  implicit none
  save

  !****************************************************************************
  !Parameters to set
  !****************************************************************************
  real, parameter :: pi          = 3.14159265358979
  real, parameter :: alpha       = 3.16 !3.73
  real, parameter :: gamma       = (2.0 * pi) / alpha
  real, parameter :: eta         = 0.5
  real, parameter :: Q           = 0.0
  real            :: Re1         = 150.0
  real, parameter :: Re_incr     = 1.0
  real, parameter :: growth_tol  = 1E-8
  real, parameter :: Re2         = 0.0 !-1.0*(1.0/eta)*Re1
  real, parameter :: Re1_mod     = 0.0
  real, parameter :: Re2_mod     = 0.0
  real, parameter :: om1         = 0.0
  real, parameter :: om2         = 0.0
  real, parameter :: dt          = 0.0001
  real, parameter :: seed        = 1E-3
  real, parameter :: end_time    = 20.0
  real, parameter :: tau_init    = 0.0
  real, parameter :: tau_step    = 1.0
  real, parameter :: tau_end     = 1.0
  integer, parameter :: nx_init     = 40
  integer, parameter :: nt_init     = 40
  integer, parameter :: nz_init     = 80
  integer, parameter :: save_rate   = 10
  integer, parameter :: save_rate_2 = 1000
  integer, parameter :: periodic_save = 10
  logical, parameter :: rot_ends    = .false.
  logical, parameter :: xsect_save  = .false.
  logical, parameter :: save3d      = .true.
  logical, parameter :: iso_hel     = .false.
  logical, parameter :: restart     = .false.
  logical, parameter :: auto_tau    = .false.
  logical, parameter :: auto_Re     = .false.
  logical, parameter :: dec_Re      = .false.
  logical, parameter :: hyst_Re     = .false.
  logical, parameter :: divergence  = .false.
  !****************************************************************************

  !****************************************************************************
  !Rarely used parameters
  !****************************************************************************
  real, parameter :: eps1      = 0.0 !0.30529
  real, parameter :: freq1     = 0.0 !2.0 * pi / 2.8
  real, parameter :: eps2      = 0.0
  real, parameter :: freq2     = 0.0
  integer, parameter :: num_pars  = 50
  logical, parameter :: save_part = .false.
  !****************************************************************************

  !****************************************************************************
  !NOTHING BELOW HERE SHOULD BE CHANGED
  !****************************************************************************
  integer, parameter :: nx     = nx_init-1
  integer, parameter :: nt     = nt_init-1
  integer, parameter :: nz     = nz_init-1
  integer, parameter :: Ntot = end_time / dt
  real, parameter :: delx = 1.0 / nx
  real, parameter :: delt = 2.0*pi / nt
  real, parameter :: delz = gamma / nz
  real, parameter :: dx2 = delx ** 2
  real, parameter :: dz2 = delz ** 2
  real, parameter :: rx = dt / delx
  real, parameter :: rz = dt / delz
  real, parameter :: rxx = dt / dx2
  real, parameter :: rzz = dt / dz2
  real, parameter :: one_eta = 1.0 - eta
  integer, parameter :: nx1 = nx-1
  integer, parameter :: nz1 = nz-1
  integer, parameter :: nxp1 = nx+1
  integer, parameter :: ntp1 = nt+1
  integer, parameter :: nzp1 = nz+1
  integer, parameter :: xlb = 1
  integer, parameter :: zlb = 0
  real :: tau = tau_init
  logical :: init_Re = .true.
  logical :: zero_Re = .false.
  logical :: gm_set = .false.
  logical :: gp_set = .false.
 
  integer :: end_proc = 0
  integer :: saturated = 0
  !****************************************************************************

end module parameters
