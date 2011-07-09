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

module stream
  !Routines to do with solving the Poisson equation for the stream function.
  !Algorithms for the solution of the stream function Poisson equation are as
  !for the current in current.f90
  implicit none

  private
  public :: p_poisson

  contains

  subroutine p_poisson(Z_mat, psi, p_mat, IPIV)
    !Solve Poisson equation for the stream function, psi for all tau
    use parameters
    use ic_bc, only : p_BCS, s
    implicit none

    double precision, intent(in) :: Z_mat(0:nx,0:nz), &
                                    p_mat(2*nx1+nx1+1,nx1*nz1)
    double precision, intent(out) :: psi(0:nx,0:nz)
    integer, intent(in) :: IPIV(nx1*nz1)
    integer :: j, k, info
    double precision :: zvec(nx1*nz1)

    do k = 1, nz1
      do j = 1, nx1
        zvec(nx1*(k-1)+j) = -s(j) * dx2 * dz2 * Z_mat(j,k)
      end do
    end do

    call DGBTRS('N', nx1*nz1, nx1, nx1, 1, p_mat, 2*nx1+nx1+1, &
                 IPIV, zvec, nx1*nz1, info)
    if (info /= 0) then
      print*, 'ERROR: Solve error psi_DGBTRS, INFO=', info
      stop
    end if

    do k = 1, nz1
      do j = 1, nx1
        psi(j,k) = zvec(nx1*(k-1)+j)
      end do
    end do

    call p_BCS(psi)

    return
  end subroutine p_poisson

end module stream
