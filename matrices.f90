module matrices
  !Routines to set up the matrices for use in the solution of the equations.
  implicit none

  private
  public :: matrix_setup, psi_mat_setup, b_mat_setup, fin_b_mat_setup, &
            j_mat_setup, fin_j_mat_setup

  contains

  subroutine matrix_setup(ux, uz, zx, zz)
    !Set up the upper, lower and diagonal parts of the tridiagonal matrices
    !for the solution of the azimuthal velocity and vorticity equations
    use parameters
    use variables
    use ic_bc, only : s
    implicit none

    type (mat_comp),    intent(out) :: ux, zx
    type (uz_mat_comp), intent(out) :: uz
    type (zz_mat_comp), intent(out) :: zz

    ux%di(:) = 1.0 + rxx + one_eta**2 * dt * 0.5 / s(1:nx1)**2
    ux%lo(:) = -0.5 * rxx + one_eta * rx * 0.25 / s(2:nx1)
    ux%up(:) = -0.5 * rxx - one_eta * rx * 0.25 / s(1:nx-2)

    uz%di(:) = 1.0 + rzz
    uz%lo(:) = -0.5 * rzz
    uz%up(:) = -0.5 * rzz

    if (abs(tau - 1.0) > epsilon(tau)) then
      uz%di(0) = 1.0 + rzz + (rz * tau / (1.0 - tau))  !extra entries -
      uz%di(nz) = 1.0 + rzz + (rz * tau / (1.0 - tau)) !Neumann BCS ends
      uz%lo(nz) = -rzz
      uz%up(0) = -rzz
    end if

    zx%di(:) = 1.0 + rxx + one_eta**2 * dt * 0.5 / s(1:nx1)**2
    zx%lo(:) = -0.5 * rxx + one_eta * rx * 0.25 / s(2:nx1)
    zx%up(:) = -0.5 * rxx - one_eta * rx * 0.25 / s(1:nx-2)

    zz%di(:) = 1.0 + rzz
    zz%lo(:) = -0.5 * rzz
    zz%up(:) = -0.5 * rzz

    return
  end subroutine matrix_setup

  subroutine psi_mat_setup(p_mat, IPIV)
    !Setup of LHS matrix in solution of stream function Poisson equation
    use parameters
    use ic_bc, only : s
    implicit none

    real, intent(out) :: p_mat(2*nx1+nx1+1,nx1*nz1)
    integer, intent(out) :: IPIV(nx1*nz1)
    integer :: j, k, info
    real :: alp(0:nx), gam(0:nx), beta, delta

    alp = dz2 + 0.5 * delx * dz2 * one_eta / s   !coefficients
    gam = dz2 - 0.5 * delx * dz2 * one_eta / s   !in matrix
    beta = -2.0 * (dz2 + dx2)
    delta = dx2

    do j = 1, nx1*nz1
      p_mat(2*nx1+1,j) = beta   !diagonal
    end do

    do j = 1, nx1*nz1-1
      p_mat(2*nx1,j+1) = gam(mod(j, nx1))   !upper-diagonal
    end do

    do j = nx1, nx1*nz1-nx1, nx1
      p_mat(2*nx1,j+1) = 0.0   !upper-diagonal, BCS
    end do

    do j = 2, nx1*nz1
       p_mat(2*nx1+2,j-1) = alp(mod(j-1, nx1) + 1) !lower-diagonal
    end do
    
    do j = nx, nx1*nz1-nx1+1, nx1
       p_mat(2*nx1+2,j-1) = 0.0   !lower-diagonal, BCS
    end do
    
    do j = 1, nx1*nz1-nx1
       p_mat(2*nx1+1-nx1,j+nx1) = delta  !upper band
    end do
    
    do j = nx, nx1*nz1
       p_mat(2*nx1+1+nx1,j-nx1) = delta  !lower band
    end do

    !LU factorisation of the matrix for use in D/SGBTRS to solve
    call SGBTRF(nx1*nz1, nx1*nz1, nx1, nx1, p_mat, 2*nx1+nx1+1, IPIV, info)
    if (info /= 0) then
      print*, 'ERROR: SPr factorisation error psi_SGBTRF, INFO=', info
      stop
    end if

    return
  end subroutine psi_mat_setup

  subroutine b_mat_setup(b_mat, IPIV)
    !Setup of LHS matrix in solution of magnetic Poisson equation.
    !Algorithm as for stream-function above.
    use parameters
    use ic_bc, only : s
    implicit none

    real, intent(out) :: b_mat(2*nxp1+nxp1+1,0:nxp1*nz1-1)
    integer, intent(out) :: IPIV(nxp1*nz1)
    integer :: j, k, info
    real :: alp(0:nx), beta(0:nx), gam(0:nx), delta

    alp(:) = dz2 - 0.5 * delx * dz2 * one_eta / s(:)
    beta(:) = -2.0 * (dz2 + dx2) - dx2 * dz2 * one_eta**2 / s(:)**2
    gam(:) = dz2 + 0.5 * delx * dz2 * one_eta / s(:)
    delta = dx2

    do j = 0, nxp1*nz1-1
      b_mat(2*nx+3,j) = beta(mod(j, nxp1))  !diagonal
    end do
    
    do j = 0, nxp1*nz1-2
      b_mat(2*nx+2,j+1) = gam(mod(j, nxp1))  !upper-diagonal
    end do
    
    do j = 1, nxp1*nz1-1
      b_mat(2*nx+4,j-1) = alp(mod(j, nxp1))  !lower-diagonal
    end do
    
    do j = 0, nxp1*nz1-nxp1-1
      b_mat(nx+2,j+nxp1) = delta  !lower band
    end do
    
    do j = nxp1, nxp1*nz1-1
      b_mat(3*nx+4,j-nxp1) = delta  !upper band
    end do
    
    do j = nx, nxp1*nz1-nxp1-1, nxp1
      b_mat(2*nx+2,j+1) = 0.0  
    end do
    
    do j = nxp1, nxp1*nz1-nxp1, nxp1
      b_mat(2*nx+4,j-1) = 0.0
    end do
    
    do j = 0, nxp1*nz1-nxp1, nxp1
      b_mat(2*nx+3,j) = (2.0 * alp(0) * delx * (1.0 - eta) / &
                         s(0)) + beta(0)
    end do
    
    do j = nx, nxp1*nz1-1, nxp1
      b_mat(2*nx+3,j) = (-2.0 * gam(nx) * delx * (1.0 - eta) / &
                         s(nx)) + beta(nx)
    end do
    
    do j = 0, nxp1*nz1-nxp1, nxp1
      b_mat(2*nx+2,j+1) = alp(0) + gam(0)
    end do
    
    do j = nx, nxp1*nz1-1, nxp1
      b_mat(2*nx+4,j-1) = alp(nx) + gam(nx)
    end do

    call SGBTRF(nxp1*nz1, nxp1*nz1, nxp1, nxp1, b_mat, 2*nxp1+nxp1+1, &
                IPIV, info)
    if (info /= 0) then
      print*, 'ERROR: SPr factorisation error mag_inf_SGBTRF, INFO=', info
      stop
    end if

    return
  end subroutine b_mat_setup

  subroutine fin_b_mat_setup(b_mat, IPIV)
    use parameters
    use ic_bc, only : s
    implicit none

    real, intent(out) :: b_mat(2*nxp1+nxp1+1,0:nxp1*nzp1-1)
    integer, intent(out) :: IPIV(nxp1*nzp1)
    integer :: j, k, info
    real :: alp(0:nx), beta(0:nx), gam(0:nx), delta

    alp(:) = dz2 - 0.5 * delx * dz2 * one_eta / s(:)
    beta(:) = -2.0 * (dz2 + dx2) - dx2 * dz2 * one_eta**2 / s(:)**2
    gam(:) = dz2 + 0.5 * delx * dz2 * one_eta / s(:)
    delta = dx2

    !diagonal
    do j = 0, nxp1*nzp1-1
      b_mat(2*nx+3,j) = beta(mod(j, nxp1))
    end do
    
    !upper diagonal
    do j = 0, nxp1*nzp1-2
      b_mat(2*nx+2,j+1) = gam(mod(j, nxp1))
    end do
    
    !lower diagonal
    do j = 1, nxp1*nzp1-1
      b_mat(2*nx+4,j-1) = alp(mod(j, nxp1))
    end do
    
    !upper diagonal branch
    do j = 0, nxp1*nzp1-nxp1-1
      b_mat(nx+2,j+nxp1) = delta
    end do
    
    !lower diagonal branch
    do j = nxp1, nxp1*nzp1-1
      b_mat(3*nx+4,j-nxp1) = delta
    end do
    
    !upper diagonal, j=nx, so j=0 not present
    do j = nx, nxp1*nzp1-nxp1-1, nxp1
      b_mat(2*nx+2,j+1) = 0.0
    end do
    
    !lower diagonal, j=0, so j=nx not present
    do j = nxp1, nxp1*nzp1-nxp1, nxp1
      b_mat(2*nx+4,j-1) = 0.0
    end do
    
    !diagonal, j=0
    do j = 0, nxp1*nzp1-nxp1, nxp1
      b_mat(2*nx+3,j) = (2.0 * alp(0) * delx * (1.0 - eta) / &
                         s(0)) + beta(0)
    end do
    
    !diagonal, k=0
    do j = 0, nx
      b_mat(2*nx+3,j) = beta(mod(j, nxp1)) - &
                        2.0 * delta * delz * (1.0 - tau) / tau
    end do
    
    !diagonal, j=nx
    do j = nx, nxp1*nzp1-1, nxp1
      b_mat(2*nx+3,j) = (-2.0 * gam(nx) * delx * (1.0 - eta) / s(nx)) + &
                         beta(nx)
    end do
    
    !diagonal, k=nz
    do j = nxp1*nzp1-nxp1, nxp1*nzp1-1
      b_mat(2*nx+3,j) = beta(mod(j, nxp1)) - &
                        2.0 * delta * delz * (1.0 - tau) / tau
    end do
    
    !diagonal, j=0, k=0
    b_mat(2*nx+3,0) = beta(0) - &
                      2.0 * delta * delz * (1.0 - tau) / tau + &
                     (2.0 * alp(0) * delx * (1.0 - eta) / s(0))
    
    !diagonal, j=0, k=nz
    b_mat(2*nx+3,nxp1*nzp1-nxp1) = beta(0) - &
                                   2.0 * delta * delz * (1.0 - tau) / &
                                   tau + &
                                  (2.0 * alp(0) * delx * (1.0 - eta) / &
                                  s(0))
    
    !diagonal, j=nx, k=0
    b_mat(2*nx+3,nx) = beta(nx) - &
                      2.0 * delta * delz * (1.0 - tau) / tau - &
                     (2.0 * gam(nx) * delx * (1.0 - eta) / s(nx))
    
    !diagonal, j=nx, k=nz
    b_mat(2*nx+3,nxp1*nzp1-1) = beta(nx) - &
                                2.0 * delta * delz * (1.0 - tau) / &
                                tau - &
                               (2.0 * gam(nx) * delx * (1.0 - eta) / &
                                s(nx))
    
    !upper diagonal, j=0
    do j = 0, nxp1*nzp1-nxp1, nxp1
      b_mat(2*nx+2,j+1) = alp(0) + gam(0)
    end do
    
    !lower diagonal, j=nx
    do j = nx, nxp1*nzp1-1, nxp1
      b_mat(2*nx+4,j-1) = alp(nx) + gam(nx)
    end do
    
    !upper diagonal branch, k=0
    do j = 0, nx
      b_mat(nx+2,j+nxp1) = 2.0 * delta
    end do
    
    !lower diagonal branch, k=nz
    do j = nxp1*nzp1-nxp1, nxp1*nzp1-1
      b_mat(3*nx+4,j-nxp1) = 2.0 * delta
    end do

    call SGBTRF(nxp1*nzp1, nxp1*nzp1, nxp1, nxp1, b_mat, 2*nxp1+nxp1+1, &
                IPIV, info)
    if (info /= 0) then
      print*, 'ERROR: SPr factorisation error mag_fin_SGBTRF, INFO=', info
      stop
    end if

    return
  end subroutine fin_b_mat_setup

  subroutine j_mat_setup(j_mat, IPIV)
    use parameters
    use ic_bc, only : s
    implicit none

    real, intent(out) :: j_mat(2*nx1+nx1+1,nx1*nzp1)
    integer, intent(out) :: IPIV(nx1*nzp1)
    integer :: j, k, info
    real :: alp(0:nx), beta(0:nx), gam(0:nx), delta

    alp(:) = dz2 - 0.5 * delx * dz2 * one_eta / s(:)
    beta(:) = -2.0 * (dz2 + dx2) - dx2 * dz2 * one_eta**2 / s(:)**2
    gam(:) = dz2 + 0.5 * delx * dz2 * one_eta / s(:)
    delta = dx2

    do j = 1, nx1*nzp1
      j_mat(2*nx-1,j) = beta(mod(j-1, nx1)+1)
    end do
    
    do j = 1, nx1
      j_mat(2*nx-1,j) = beta(mod(j-1, nx1)+1) - &
                        2.0 * delta * delz * tau / (1.0 - tau)
    end do
    
    do j = nx1*nzp1-nx1+1, nx1*nzp1
      j_mat(2*nx-1,j) = beta(mod(j-1, nx1)+1) - &
                        2.0 * delta * delz * tau / (1.0 - tau)
    end do
    
    do j = 1, nx1*nzp1-1
      j_mat(2*nx-2,j+1) = gam(mod(j, nx1))
    end do
    
    do j = nx1, nx1*nzp1-nx1, nx1
      j_mat(2*nx-2,j+1) = 0.0
    end do
    
    do j = 2, nx1*nzp1
      j_mat(2*nx,j-1) = alp(mod(j-1, nx1) + 1)
    end do
    
    do j = nx, nx1*nzp1-nx1+1, nx1
      j_mat(2*nx,j-1) = 0.0
    end do
    
    do j = 1, nx1*nzp1-nx1
      j_mat(nx,j+nx1) = delta
    end do
    
    do j = nx, nx1*nzp1
      j_mat(3*nx-2,j-nx1) = delta
    end do
    
    do j = 1, nx1
      j_mat(nx,j+nx1) = 2.0 * delta
    end do
    
    do j = nx1*nzp1-nx1+1, nx1*nzp1
      j_mat(3*nx-2,j-nx1) = 2.0 * delta
    end do

    call SGBTRF(nx1*nzp1, nx1*nzp1, nx1, nx1, j_mat, 2*nx1+nx1+1, IPIV, &
                info)
    if (info /= 0) then
      print*, 'ERROR: SPr factorisation error cur_inf_SGBTRF, INFO=', info
      stop
    end if

    return
  end subroutine j_mat_setup

  subroutine fin_j_mat_setup(j_mat, IPIV)
    use parameters
    use ic_bc, only : s
    implicit none

    real, intent(out) :: j_mat(2*nx1+nx1+1,nx1*nz1)
    integer, intent(out) :: IPIV(nx1*nz1)
    integer :: j, k, info
    real :: alp(0:nx), beta(0:nx), gam(0:nx), delta

    alp(:) = dz2 - 0.5 * delx * dz2 * one_eta / s(:)
    beta(:) = -2.0 * (dz2 + dx2) - dx2 * dz2 * one_eta**2 / s(:)**2
    gam(:) = dz2 + 0.5 * delx * dz2 * one_eta / s(:)
    delta = dx2

    do j = 1, nx1*nz1
      j_mat(2*nx-1,j) = beta(mod(j-1, nx1)+1)
    end do
    
    do j = 1, nx1*nz1-1
      j_mat(2*nx-2,j+1) = gam(mod(j, nx1))
    end do
    
    do j = nx1, nx1*nz1-nx1, nx1
      j_mat(2*nx-2,j+1) = 0.0
    end do
    
    do j = 2, nx1*nz1
      j_mat(2*nx,j-1) = alp(mod(j-1, nx1) + 1)
    end do
    
    do j = nx, nx1*nz1-nx1+1, nx1
      j_mat(2*nx,j-1) = 0.0
    end do
    
    do j = 1, nx1*nz1-nx1
      j_mat(nx,j+nx1) = delta
    end do
    
    do j = nx, nx1*nz1
      j_mat(3*nx-2,j-nx1) = delta
    end do

    call SGBTRF(nx1*nz1, nx1*nz1, nx1, nx1, j_mat, 2*nx1+nx1+1, IPIV, info)
    if (info /= 0) then
      print*, 'ERROR: SPr factorisation error cur_fin_SGBTRF, INFO=', info
      stop
    end if

    return
  end subroutine fin_j_mat_setup

end module matrices
