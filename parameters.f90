MODULE parameters
implicit none
save

double precision, parameter :: pi 	   = 3.14159265358979d0
double precision, parameter :: alpha 	   = 0d0 !3.13d0 / 4d0
double precision, parameter :: gamma 	   = 40d0 !(2d0 * pi) / alpha
double precision, parameter :: eta 	   = 0.75d0
double precision, parameter :: Re1 	   = 100d0
double precision, parameter :: Re2 	   = 0d0
double precision, parameter :: Re1_mod 	   = 0d0
double precision, parameter :: Re2_mod 	   = 0d0
double precision, parameter :: om1 	   = 0d0
double precision, parameter :: om2 	   = 0d0
double precision, parameter :: dt 	   = 0.001d0
double precision, parameter :: seed 	   = 1d-1
double precision, parameter :: end_time    = 1000d0
double precision, parameter :: tau_init	   = 1d0
double precision, parameter :: tau_step    = 1d0
double precision, parameter :: tau_end     = 1d0
integer, 	  parameter :: nx 	   = 20
integer, 	  parameter :: nz 	   = 800
integer, 	  parameter :: save_rate   = 10
integer, 	  parameter :: save_rate_2 = 35
logical, 	  parameter :: diag 	   = .false.
logical, 	  parameter :: xsect_save  = .false.
logical, 	  parameter :: restart 	   = .false.
logical,	  parameter :: auto_tau    = .false.
logical,	  parameter :: save_part   = .false.
double precision, parameter :: eps1	   = 50d0
double precision, parameter :: freq1	   = 3.13d0
double precision, parameter :: eps2	   = 0d0
double precision, parameter :: freq2	   = 0d0
double precision, parameter :: x_par_pos   = 0.5d0
double precision, parameter :: z_par_pos   = 0.3d0

integer, 	  parameter :: Ntot = end_time / dt
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
double precision :: tau = tau_init
double precision :: x_pos = x_par_pos * nx
double precision :: z_pos = z_par_pos * nz

!contains

!SUBROUTINE get_param()
!implicit none

!rx = dt / delx
!rz = dt / delz
!rxx = dt / dx2
!rzz = dt / dz2

!END SUBROUTINE get_param

END MODULE parameters
