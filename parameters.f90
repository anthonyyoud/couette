MODULE parameters
  !Parameters to set.
  IMPLICIT NONE
  SAVE

  !****************************************************************************
  !Data types by kind
  !****************************************************************************
  INTEGER,      PARAMETER :: i1 = SELECTED_INT_KIND(9)       !4-byte integer
  !INTEGER (i1), PARAMETER :: r2 = SELECTED_REAL_KIND(6,37)   !4-byte real
  INTEGER (i1), PARAMETER :: r2 = SELECTED_REAL_KIND(15,307) !8-byte real
  !****************************************************************************

  !****************************************************************************
  !Process rows and columns for ScaLAPACK
  !****************************************************************************
  INTEGER (i1), PARAMETER :: nprow       = 1
  INTEGER (i1), PARAMETER :: npcol       = 1
  !INTEGER (i1), PARAMETER :: nb          = 10151
  !****************************************************************************

  !****************************************************************************
  !Parameters to set
  !****************************************************************************
  REAL    (r2), PARAMETER :: pi          = 3.14159265358979_r2
  REAL    (r2), PARAMETER :: alpha       = 0.0_r2 !3.1631_r2
  REAL    (r2), PARAMETER :: gamma       = 12.30_r2 !(2.0_r2 * pi) / alpha
  REAL    (r2), PARAMETER :: eta         = 0.75_r2
  REAL    (r2), PARAMETER :: Q           = 0.0_r2
  REAL    (r2)            :: Re1         = 0.0_r2
  REAL    (r2), PARAMETER :: Re_incr     = 1.0_r2
  REAL    (r2), PARAMETER :: growth_tol  = 1E-8_r2
  REAL    (r2), PARAMETER :: Re2         = 0.0_r2 !-1.0_r2*(1.0_r2/eta)*Re1
  REAL    (r2), PARAMETER :: Re1_mod     = 100.0_r2
  REAL    (r2), PARAMETER :: Re2_mod     = 0.0_r2
  REAL    (r2), PARAMETER :: om1         = 3.0_r2
  REAL    (r2), PARAMETER :: om2         = 0.0_r2
  REAL    (r2), PARAMETER :: dt          = 0.0001_r2
  REAL    (r2), PARAMETER :: seed        = 1E-1_r2
  REAL    (r2), PARAMETER :: end_time    = 1000.0_r2
  REAL    (r2), PARAMETER :: tau_init    = 1.0_r2
  REAL    (r2), PARAMETER :: tau_step    = 1.0_r2
  REAL    (r2), PARAMETER :: tau_end     = 1.0_r2
  INTEGER (i1), PARAMETER :: nx          = 20
  INTEGER (i1), PARAMETER :: nt          = 20
  INTEGER (i1), PARAMETER :: nz          = 246
  INTEGER (i1), PARAMETER :: save_rate   = 10
  INTEGER (i1), PARAMETER :: save_rate_2 = 104
  LOGICAL,      PARAMETER :: xsect_save  = .FALSE.
  LOGICAL,      PARAMETER :: save3d      = .FALSE.
  LOGICAL,      PARAMETER :: iso_hel     = .FALSE.
  LOGICAL,      PARAMETER :: restart     = .FALSE.
  LOGICAL,      PARAMETER :: auto_tau    = .FALSE.
  LOGICAL,      PARAMETER :: auto_Re     = .FALSE.
  LOGICAL,      PARAMETER :: dec_Re      = .FALSE.
  LOGICAL,      PARAMETER :: hyst_Re     = .FALSE.
  !****************************************************************************

  !****************************************************************************
  !Rarely used parameters
  !****************************************************************************
  INTEGER (i1), PARAMETER :: nb        = (nx+1)*(nz+1)

  REAL    (r2), PARAMETER :: eps1      = 0.0_r2 !0.30529_r2
  REAL    (r2), PARAMETER :: freq1     = 0.0_r2 !2.0_r2 * pi / 2.8_r2
  REAL    (r2), PARAMETER :: eps2      = 0.0_r2
  REAL    (r2), PARAMETER :: freq2     = 0.0_r2
  INTEGER (i1), PARAMETER :: num_pars  = 50_i1
  LOGICAL,      PARAMETER :: save_part = .FALSE.
  !****************************************************************************

  !****************************************************************************
  !NOTHING BELOW HERE SHOULD BE CHANGED
  !****************************************************************************
  INTEGER (i1), PARAMETER :: Ntot = end_time / dt
  REAL    (r2), PARAMETER :: delx = 1.0_r2 / nx
  REAL    (r2), PARAMETER :: delt = 2.0_r2*pi / nt
  REAL    (r2), PARAMETER :: delz = gamma / nz
  REAL    (r2), PARAMETER :: dx2 = delx ** 2
  REAL    (r2), PARAMETER :: dz2 = delz ** 2
  REAL    (r2), PARAMETER :: rx = dt / delx
  REAL    (r2), PARAMETER :: rz = dt / delz
  REAL    (r2), PARAMETER :: rxx = dt / dx2
  REAL    (r2), PARAMETER :: rzz = dt / dz2
  REAL    (r2), PARAMETER :: one_eta = 1.0_r2 - eta
  INTEGER (i1), PARAMETER :: nx1 = nx-1
  INTEGER (i1), PARAMETER :: nz1 = nz-1
  INTEGER (i1), PARAMETER :: nxp1 = nx+1
  INTEGER (i1), PARAMETER :: nzp1 = nz+1
  INTEGER (i1), PARAMETER :: xlb = 1
  INTEGER (i1), PARAMETER :: zlb = 0
  REAL    (r2)            :: tau = tau_init
  LOGICAL (r2)            :: init_Re = .TRUE.
  LOGICAL (r2)            :: zero_Re = .FALSE.
  LOGICAL (r2)            :: gm_set = .FALSE.
  LOGICAL (r2)            :: gp_set = .FALSE.
 
  INTEGER (i1), PARAMETER :: laf   = nb*(nx1+nx1)+6*nx1*nx1
  INTEGER (i1), PARAMETER :: b_laf = nb*(nxp1+nxp1)+6*nxp1*nxp1
  INTEGER (i1), PARAMETER :: lwork_fac   = nx1*nx1
  INTEGER (i1), PARAMETER :: lwork_b_fac = nxp1*nxp1
  INTEGER (i1), PARAMETER :: lwork_sol   = nx1
  INTEGER (i1), PARAMETER :: lwork_b_sol = nxp1
  INTEGER (i1) :: ictxt
  INTEGER (i1) :: myrow
  INTEGER (i1) :: mycol
  INTEGER (i1) :: b_M
  INTEGER (i1) :: b_N
  INTEGER (i1) :: j_M
  INTEGER (i1) :: j_N
  INTEGER (i1) :: p_M
  INTEGER (i1) :: p_N
  INTEGER (i1) :: end_proc = 0
  INTEGER (i1) :: saturated = 0
  INTEGER (i1), PARAMETER :: SPr = SELECTED_REAL_KIND(6,37) !4-byte real
  INTEGER (i1), PARAMETER :: DPr = SELECTED_REAL_KIND(15,307) !8-byte real
  !****************************************************************************

END MODULE parameters
