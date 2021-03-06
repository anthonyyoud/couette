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

module derivs
  !Routines which calculate finite-difference derivatives.  These only
  !calculate the numerator of the derivatives; the denominator dx, dz etc. is
  !included in other routines whenever the derivative is used.
  implicit none

  private
  public :: deriv_x, deriv_z, deriv_xx, deriv_zz

  contains

  subroutine deriv_x(f, fx)
    !First x derivative
    use parameters
    implicit none

    double precision, intent(in)  :: f(0:nx,0:nz)
    double precision, intent(out) :: fx(0:nx,0:nz)
    integer :: j, k

    do j = 1, nx1
      fx(j,:) = f(j+1,:) - f(j-1,:)
    end do

    return
  end subroutine deriv_x

  subroutine deriv_xx(f, fxx)
    !Second x derivative
    use parameters
    implicit none

    double precision, intent(in)  :: f(0:nx,0:nz)
    double precision, intent(out) :: fxx(0:nx,0:nz)
    integer :: j, k

    do j = 1, nx1
      fxx(j,:) = f(j+1,:) - 2d0 * f(j,:) + f(j-1,:)
    end do

    return
  end subroutine deriv_xx

  subroutine deriv_z(f, fz)
    !First z derivative
    use parameters
    implicit none

    double precision, intent(in)  :: f(0:nx,0:nz)
    double precision, intent(out) :: fz(0:nx,0:nz)
    integer :: j, k

    do k = 1, nz1
      fz(1:nx1,k) = f(1:nx1,k+1) - f(1:nx1,k-1)
      fz(1:nx1,0) = -f(1:nx1,2) + 4d0 * f(1:nx1,1) - 3d0 * f(1:nx1,0)
      fz(1:nx1,nz) = f(1:nx1,nz-2) - 4d0 * f(1:nx1,nz1) + 3d0 * &
                     f(1:nx1,nz)   !forward/backward difference at boundaries
    end do

    return
  end subroutine deriv_z

  subroutine deriv_zz(f, fzz)
    !Second z derivative
    use parameters
    implicit none

    double precision, intent(in)  :: f(0:nx,0:nz)
    double precision, intent(out) :: fzz(0:nx,0:nz)
    integer :: j, k

    do k = 1, nz1
      fzz(1:nx1,k) = f(1:nx1,k+1) - 2d0 * f(1:nx1,k) + f(1:nx1,k-1)
      fzz(1:nx1,0) = -f(1:nx1,3) + 4d0 * f(1:nx1,2) - &
                     5d0 * f(1:nx1,1) + 2d0 * f(1:nx1,0) !forward/backward
      fzz(1:nx1,nz) = -f(1:nx1,nz-3) + 4d0 * f(1:nx1,nz-2) - & !difference at
                      5d0 * f(1:nx1,nz1) + 2d0 * f(1:nx1,nz)  !boundaries
    end do

    return
  end subroutine deriv_zz

end module derivs
