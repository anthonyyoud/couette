Description
===========
Axisymmetric, hydromagnetic, modulated Couette flow between infinite or finite
cylinders in the small Prandtl number limit.  Time-stepping is via 2nd order
accurate implicit Crank-Nicolson for the linear terms and 2nd order accurate
explicit Adams-Bashforth for the non-linear terms.  For the diffusive equations
the code uses operator factorisation to allow a tridiagonal system to be
solved; the linear algebra package LAPACK is used to solve the Poisson
equations associated with the stream function, current and magnetic field.  The
spatial discretisation is via 2nd order accurate centred finite differences.
It includes a 'homotopy' parameter, `tau`, to continuously deform the
boundaries from the infinite cylinder case (`tau = 0`) to the finite cylinder
case (`tau = 1`).  A basic particle path subroutine to track the trajectory of
a particle is also included.  Spatial modulation of the inner and/or outer
Reynolds numbers in the axial direction is possible to mimic a wavy cylinder
boundary.

Requirements
============
* LAPACK - the linear algebra package.
* BLAS - the basic linear algebra subprograms.
* A Fortran 90/95 compiler.

Files
=====
* __couette_mod.f90__ - Main program file.
* __current.f90__ - Routines to solve the Poisson equation for the azimuthal current.
* __derivs.f90__ - Routines to calculate derivatives.
* __ic_bc.f90__ - Routines to set up initial and boundary conditions.
* __io.f90__ - Routines to do with input/output.
* __linear.f90__ - Routines to set up the linear parts of the RHS of the azimuthal velocity and vorticity fields.
* __magnetic.f90__ - Routines to solve the Poisson equation for the azimuthal magnetic field.
* __Makefile__ - Makefile to compile the code (see below).
* __matrices.f90__ - Routines to set up matrices involved in time-stepping and solving Poisson equations.
* __nonlinear.f90__ - As linear.f90 but for the nonlinear terms.
* __parameters.f90__ - Parameters to set (see below).
* __README.md__ - This file.
* __setup__ - Setup script to compile and set up ready to run in separate directory.
* __solve.f90__ - Routines to solve for the azimuthal velocity and vorticity fields, as well as the Thomas algorithm.
* __stream.f90__ - Routines to solve Poisson equation for stream function.
* __variables.f90__ - Derived types and routines to do with variables.

Main parameters
===============
* __alpha__ - Wavenumber (`2pi/wavelength`) in the infinite case.  Set equal zero if using finite cylinders.
* __gamma__ - Aspect ratio (`height/gap`) for the finite case.  Set equal to `2pi/alpha` if using infinite cylinders.
* __eta__ - Radius ratio (`R1/R2`).
* __Q__ - Chandrasekhar number.  Measure of the imposed magnetic field.
* __Re1,2__ - Reynolds number in `Re1,2(t)=Re1,2+Re1,2mod*cos(om1,2*t)`.
* __Re1,2_mod__ - Modulation amplitude as above.
* __om1,2__ - Frequency of modulation as above.
* __Re_incr__ - How much `Re1` should be incremented or decremented when searching for critical values.
* __growth_tol__ - How close two successive growth rates should be before `Re1` is altered in searching for critical values.
* __dt__ - Time step.  In general `10^-4` is a good choice.  Once the Reynolds numbers are large (>700, say) and/or the spatial resolution is increased significantly (>80 points in z) then `dt = 10^-5` is a better choice for stability.
* __seed__ - Initial seed for initial conditions.  Set small `O(10^-10)` for calculating linear growth rates.  Set `O(10^-3)` for non-linear saturation.  In practice this can be set to zero, since the boundary conditions of the azimuthal velocity at the cylinder wall(s) can start the flow.
* __end_time__ - Final viscous diffusion time.
* __tau_init__ - Initial value of homotopy parameter, tau.  `0<=tau<=1`.  `tau = 0` => infinite cylinders, `tau = 1` => finite cylinders.
* __tau_step__ - For steady case, how much tau should be increased after saturation at each tau.
* __tau_end__ - For steady case, the final value of tau desired.
* __nx_init__ - Number of radial grid points.
* __nt_init__ - Number of azimuthal grid points for use when a 3D OpenDX isosurface is desired.  This is not actually used in any computation in the code other than for the isosurface plots.
* __nz_init__ - Number of axial grid points.  For infinite `2*nx` is sufficient.  For finite `gamma*nx` is usually required for full resolution.
* __save_rate__ - After how many time steps should velocities be saved?
* __save_rate_2__ - After how many time steps should cross-sections be saved?
* __xsect_save__ - Should cross-sections (surfaces) of fields be saved?
* __save3d__ - Should a 3D isosurface be saved (OpenDX)?
* __iso_hel__ - If `save3d = .true.` should the isosurface be of the helicity?  If `iso_hel = .false.` then the stream function is saved.
* __restart__ - Should we restart from a previous run?  If `.true.`, file `end_state.dat` should be copied to restart directory.
* __auto_tau__ - For steady case, should tau be automatically stepped after saturation at each tau.  Set in conjunction with `tau_step` and `tau_end`.
* __auto_Re__ - For steady case, should `Re1` be stepped automatically?
* __dec_Re__ - Can the critical value for `Re1` only be found by a quasi-static decrease of `Re1`?
* __hyst_Re__ - Is the critical value of `Re1` in a hysteresis region?  (Specifically for 1- and 2-cell flows at very small aspect ratio).

Rarely used parameters
======================
* __eps1,2__ - Amplitude of oscillation of `Re_1,2(t,z)` in axial direction.
* __freq1,2__ - Frequency of oscillation of `Re_1,2(t,z)` in axial direction.
* __x\_par\_pos__ - Initial radial position of a particle in the fluid as a fraction of gap width.
* __z\_par\_pos__ - Initial axial position of a particle in the fluid as a fraction of `gamma` (finite) or `alpha` (infinite).
* __save_part__ - Should a particle path be saved?

Makefile
========
Use the makefile to compile the code and handle module dependencies.

* OBJECT - The compiled executable name.
* OBJS - The object files that should be linked.
* FC - The Fortran compiler.
* FFLAGS - Compiler flags.
* LDFLAGS - Any extra flags required by the linker.

To run
======
Set `parameters.f90` then run

    source setup <directory>

which compiles code and copies `parameters.f90` and moves `OBJECT` to
`<directory>`.  `./couette_mod` runs code.

If restarting from a previous run then be sure to set `restart = .true.` in
`parameters.f90`, and have the file `end_state.dat` from the preceding run in
the run/submit directory.

Errors are output if either:

1. `end_state.dat` exists but `restart = .false.` or
2. `restart = .true.` but `end_state.dat` does not exist.

If, at any time during a run, the file `SAVE` is present in the run directory
(e.g. by giving the command `touch SAVE`), the program will note this and
output data files for contours and surfaces at that time.  Once the data files
are output `SAVE` is removed.  This can be done at any time and as many times
as desired during the run.

Any run should be cleanly stopped by removing the empty file `RUNNING` which is
present in the run directory once the job has started.  This avoids having to
use `Ctrl-C` to interrupt the program.  Using this procedure will ensure that
all run-time files are cleaned up and unneeded directories removed.

Output files
============
* __end_state.dat__	- Saves time index, p, time-step, dt, and fields u, Z, psi, current and magnetic field for restart.
* __energy.dat__ - Saves the kinetic energy due to CCF and the total kinetic energy including CCF at each time.
* __max_psi.dat__ - Saves maximum value of psi (and vr, vz) over whole field.  Mainly for use with IDL plots.
* __particle.dat__ - Saves the radial and axial position of a particle at each time.
* __time_tau.dat__ - If `auto_tau = .true.` saves time at which each step in tau occured.
* __torque.dat__ - Saves the torques on the inner (G1) and outer cylinders (G2), the torque due to CCF (Gc), and the ratio G1/Gc.  For steady states G1=G2.
* __u_growth.dat__ - Saves time, radial velocity (outflow, inflow), growth rate, axial velocity, stream function, azimuthal velocity, vorticity, azimuthal current and magnetic field, and Reynolds number.  Velocities are saved at the points defined in `io.f90` - subroutine `save_growth`.
* __p1234567.dat__ - Stream function field saved at the time defined by `(time-step)*(number following 'p')` in filename.  (Only if `xsect_save == .true.`).
* __u1234567.dat__ - As above but for azimuthal velocity field.
* __z1234567.dat__ - As above but for vorticity field.
* __vr1234567.dat__ - As above but for radial velocity field.
* __vz1234567.dat__ - As above but for axial velocity field.
* __j1234567.dat__ - As above but for azimuthal current field.
* __b1234567.dat__ - As above but for azimuthal magnetic field.
* __xsect1234567.dat__ - Cross-sections of all fields except vorticity saved at the time defined as above.  (Only if `xsect_save == .true.`).
