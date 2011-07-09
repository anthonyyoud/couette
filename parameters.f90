!Copyright 2011 Anthony Youd/Newcastle University
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

module parameters
  !Parameters to set.
  implicit none
  save

  !****************************************************************************
  !Parameters to set
  !****************************************************************************
  double precision, parameter :: pi          = 3.14159265358979323846d0
  double precision, parameter :: alpha       = 3.16d0 !3.73d0
  double precision, parameter :: gamma       = (2d0 * pi) / alpha
  double precision, parameter :: eta         = 0.75d0
  double precision, parameter :: Q           = 0d0
  double precision            :: Re1         = 380d0
  double precision, parameter :: Re_incr     = 1d0
  double precision, parameter :: growth_tol  = 1d-8
  double precision, parameter :: Re2         = 0d0 !-1d0*(1d0/eta)*Re1
  double precision, parameter :: Re1_mod     = 0d0
  double precision, parameter :: Re2_mod     = 0d0
  double precision, parameter :: om1         = 0d0
  double precision, parameter :: om2         = 0d0
  double precision, parameter :: dt          = 1d-4
  double precision, parameter :: seed        = 1d-1
  double precision, parameter :: end_time    = 50d0
  double precision, parameter :: tau_init    = 0d0
  double precision, parameter :: tau_step    = 1d0
  double precision, parameter :: tau_end     = 1d0
  integer, parameter :: nx_init     = 80
  integer, parameter :: nt_init     = 80
  integer, parameter :: nz_init     = 160
  integer, parameter :: save_rate   = 10
  integer, parameter :: save_rate_2 = 18
  integer, parameter :: periodic_save = 10
  logical, parameter :: rot_ends    = .false.
  logical, parameter :: xsect_save  = .false.
  logical, parameter :: save3d      = .false.
  logical, parameter :: savevapor3d = .false.
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
  double precision, parameter :: eps1      = 0d0 !0d30529
  double precision, parameter :: freq1     = 0d0 !2d0 * pi / 2.8d0
  double precision, parameter :: eps2      = 0d0
  double precision, parameter :: freq2     = 0d0
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
  double precision, parameter :: delx = 1d0 / nx
  double precision, parameter :: delt = 2d0*pi / nt
  double precision, parameter :: delz = gamma / nz
  double precision, parameter :: dx2 = delx ** 2
  double precision, parameter :: dz2 = delz ** 2
  double precision, parameter :: rx = dt / delx
  double precision, parameter :: rz = dt / delz
  double precision, parameter :: rxx = dt / dx2
  double precision, parameter :: rzz = dt / dz2
  double precision, parameter :: one_eta = 1d0 - eta
  integer, parameter :: nx1 = nx-1
  integer, parameter :: nz1 = nz-1
  integer, parameter :: nxp1 = nx+1
  integer, parameter :: ntp1 = nt+1
  integer, parameter :: nzp1 = nz+1
  integer, parameter :: xlb = 1
  integer, parameter :: zlb = 0
  double precision :: tau = tau_init
  logical :: init_Re = .true.
  logical :: zero_Re = .false.
  logical :: gm_set = .false.
  logical :: gp_set = .false.
 
  integer :: end_proc = 0
  integer :: saturated = 0
  !****************************************************************************

end module parameters
