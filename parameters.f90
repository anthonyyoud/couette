MODULE parameters
implicit none
save

double precision, parameter :: pi = 3.14159265358979d0
double precision, parameter :: alpha = 3.16d0
double precision, parameter :: gamma = (2d0 * pi) / alpha
double precision, parameter :: eta = 0.5d0
double precision, parameter :: Re1 = 72.5d0
double precision, parameter :: Re2 = 0d0
double precision, parameter :: Re1_mod = 0d0
double precision, parameter :: Re2_mod = 0d0
double precision, parameter :: om1 = 0d0
double precision, parameter :: om2 = 0d0
double precision, parameter :: dt = 0.0001d0
double precision, parameter :: seed = 1d-1
double precision, parameter :: end_time = 10d0
double precision, parameter :: tau = 0d0
integer, parameter :: Ntot = end_time / dt
integer, parameter :: save_rate = 100
integer, parameter :: save_rate_2 = 1000
integer, parameter :: nx = 20
integer, parameter :: nz = 40
logical, parameter :: diag = .false.
logical, parameter :: xsect_save = .false.
logical, parameter :: restart = .false.

double precision, parameter :: delx = 1d0 / (nx+0)
double precision, parameter :: delz = gamma / (nz+0)
double precision, parameter :: dx2 = delx ** 2
double precision, parameter :: dz2 = delz ** 2
double precision, parameter :: rx = dt / delx
double precision, parameter :: rz = dt / delz
double precision, parameter :: rxx = dt / dx2
double precision, parameter :: rzz = dt / dz2
integer, parameter :: nx1 = nx-1
integer, parameter :: nz1 = nz-1
integer, parameter :: xlb = 1
integer, parameter :: zlb = 0

!contains

!SUBROUTINE get_param()
!implicit none

!rx = dt / delx
!rz = dt / delz
!rxx = dt / dx2
!rzz = dt / dz2

!END SUBROUTINE get_param

END MODULE parameters
