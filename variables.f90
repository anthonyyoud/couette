module variables
  !Set up types and other routines to manipulate variables.
  use parameters
  implicit none

  private
  public :: copy_var, vr_vz, integrate_r, integrate_z, &
            particle_setup

  real, parameter, private   :: c1=3.0/8.0, & !Constants for
                                c2=7.0/6.0, & !numerical
                                c3=23.0/24.0  !integration

  type, public :: var
    !variables
    real :: new(0:nx,0:nz)
    real :: old(0:nx,0:nz)
    real :: old2(0:nx,0:nz)
    real :: inter(0:nx,0:nz)
    real :: nlin_new(0:nx,0:nz)
    real :: nlin_old(0:nx,0:nz)
  end type var

  type, public :: deriv
    !derivatives
    real :: x(0:nx,0:nz)
    real :: xx(0:nx,0:nz)
    real :: z(0:nx,0:nz)
    real :: zz(0:nx,0:nz)
    real :: zx(0:nx,0:nz)
    real :: zxx(0:nx,0:nz)
    real :: zzz(0:nx,0:nz)
  end type deriv

  type, public :: mat_comp
    !matrix components
    real :: lo(2:nx1)
    real :: di(nx1)
    real :: up(nx-2)
  end type mat_comp

  type, public :: uz_mat_comp
    !matrix components for u in the z-direction
    real :: lo(nz)
    real :: di(0:nz)
    real :: up(0:nz1)
  end type uz_mat_comp

  type, public :: zz_mat_comp
    !matrix components for Z in the z-direction
    real :: lo(2:nz1)
    real :: di(nz1)
    real :: up(nz-2)
  end type zz_mat_comp

  type (var), public, save :: ut, zt, psi, bt, jt
  real, public, save :: vr(0:nx,0:nz), vz(0:nx,0:nz), &
                        vrold(0:nx,0:nz) = 0.0, &
                        vzold(0:nx,0:nz) = 0.0, &
                        xold(num_pars), zold(num_pars)

  contains

  subroutine copy_var(var_out, var_in)
    !Copy a variable
    use parameters
    implicit none

    real, intent(in) :: var_in(0:nx,0:nz)
    real, intent(out) :: var_out(0:nx,0:nz)

    var_out = var_in

    return
  end subroutine copy_var

  subroutine vr_vz(p, vr, vz)
    !Calculate radial and axial velocity components from the stream function
    use parameters
    use ic_bc, only : x, z, s
    use derivs, only : deriv_x, deriv_z
    implicit none

    real, intent(in) :: p(0:nx,0:nz)
    real, intent(out) :: vr(0:nx,0:nz), vz(0:nx,0:nz)
    integer :: j, k
    real :: vrx(0:nx,0:nz), vzz(0:nx,0:nz), div

    do k = 1, nz1
      vr(:,k) = (-1.0 / (2.0 * s(:) * delz)) * (p(:,k+1) - p(:,k-1))
    end do

    if (abs(tau - 1.0) > epsilon(tau)) then
      vr(:,0) = (-1.0 / (2.0 * s(:) * delz)) * &
                (-3.0 * p(:,0) + 4.0 * p(:,1) - p(:,2))
      vr(:,nz) = (-1.0 / (2.0 * s(:) * delz)) * &
                 (3.0 * p(:,nz) - 4.0 * p(:,nz1) + p(:,nz-2))
    else
      vr(:,0) = 0.0
      vr(:,nz) = 0.0
    end if

    do j = 1, nx1
      vz(j,:) = (1.0 / (2.0 * s(j) * delx)) * (p(j+1,:) - p(j-1,:))
    end do

    vz(0,:) = 0.0 !(1.0 / (2.0 * s(0) * delx)) * &
                  !(-3.0 * p(0,k) + 4.0 * p(1,k) - p(2,k))
    vz(nx,:) = 0.0 !(1.0 / (2.0 * s(nx) * delx)) * &
                  !(3.0 * p(nx,k) - 4.0 * p(nx1,k) + p(nx-2,k))


    if (divergence) then
      !do k = 0, nz
      !  do j = 0, nx
      !    vr(j,k) = z(k)**3
      !    vz(j,k) = x(j)
      !  end do
      !end do

      call deriv_x(vr, vrx)
      call deriv_z(vz, vzz)
      
      div = 0.0
      
      do k = 0, nz
        do j = 1, nx1
          div = div + one_eta * vrx(j,k) / s(j) + 0.5 * vrx(j,k) / delx + &
                       0.5 * vzz(j,k) / delz
        end do
      end do
    
      open (97, status = 'unknown', position = 'append', &
                file = 'divergence.dat')
      write (97, '(e17.9)') div
      close (97)
    end if

    return
  end subroutine vr_vz

  subroutine integrate_r(in_var, r_int)
    !Integrate a (2D) variable in r
    use parameters
    use ic_bc, only : s
    implicit none
                                                                                
    real, intent(in) :: in_var(0:nx,0:nz)
    real, intent(out) :: r_int(0:nz)
    integer :: j, k
    real :: var(0:nx,0:nz)
                                                                                
    do j = 0, nx
      var(j,:) = in_var(j,:) * s(j) / (1.0 - eta)
    end do

    do k = 0, nz
      r_int(k) = (c1 * var(0,k) + &
                  c2 * var(1,k) + &
                  c3 * var(2,k) + &
                  sum(var(3:nx-3,k)) + &
                  c3 * var(nx-2,k) + &
                  c2 * var(nx-1,k) + &
                  c1 * var(nx,k)) * delx
    end do
                                                                                
    return
  end subroutine integrate_r

  subroutine integrate_z(in_var, z_int)
    !Integrate a (1D) variable in z
    use parameters
    implicit none
                                                                                
    real, intent(in) :: in_var(0:nz)
    real, intent(out) :: z_int
                                                                                
    z_int = (c1 * in_var(0) + &
             c2 * in_var(1) + &
             c3 * in_var(2) + &
             sum(in_var(3:nz-3)) + &
             c3 * in_var(nz-2) + &
             c2 * in_var(nz-1) + &
             c1 * in_var(nz)) * delz
                                                                                
    return
  end subroutine integrate_z

  subroutine particle_setup()
    !Set up particle positions
    use parameters
    implicit none

    integer :: j, k
    real :: x_par_pos, z_par_pos
    
    call random_seed()

    do j = 1, num_pars
      call random_number(x_par_pos)
      xold(j) = x_par_pos * nx
    end do

    call random_seed()
    
    do k = 1, num_pars
      call random_number(z_par_pos)
      zold(k) = z_par_pos * nz
    end do

    return
  end subroutine particle_setup
  
end module variables
