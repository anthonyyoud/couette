MODULE derivs
IMPLICIT NONE

PRIVATE
PUBLIC :: deriv_x, deriv_z, deriv_xx, deriv_zz

CONTAINS

SUBROUTINE deriv_x(f, fx)
!First x derivative
USE parameters
IMPLICIT NONE

REAL (r2), INTENT(IN)  :: f(0:nx,0:nz)
REAL (r2), INTENT(OUT) :: fx(0:nx,0:nz)
INTEGER (i1)           :: j, k

DO j = 1, nx1
   fx(j,:) = f(j+1,:) - f(j-1,:)
END DO

RETURN
END SUBROUTINE deriv_x

SUBROUTINE deriv_xx(f, fxx)
!Second x derivative
USE parameters
IMPLICIT NONE

REAL (r2), INTENT(IN)  :: f(0:nx,0:nz)
REAL (r2), INTENT(OUT) :: fxx(0:nx,0:nz)
INTEGER (i1)           :: j, k

DO j = 1, nx1
   fxx(j,:) = f(j+1,:) - 2.0_r2 * f(j,:) + f(j-1,:)
END DO

RETURN
END SUBROUTINE deriv_xx

SUBROUTINE deriv_z(f, fz)
!First z derivative
USE parameters
IMPLICIT NONE

REAL (r2), INTENT(IN)  :: f(0:nx,0:nz)
REAL (r2), INTENT(OUT) :: fz(0:nx,0:nz)
INTEGER (i1)           :: j, k

DO k = 1, nz1
   fz(1:nx1,k) = f(1:nx1,k+1) - f(1:nx1,k-1)
   fz(1:nx1,0) = -f(1:nx1,2) + 4.0_r2 * f(1:nx1,1) - 3.0_r2 * f(1:nx1,0)
   fz(1:nx1,nz) = f(1:nx1,nz-2) - 4.0_r2 * f(1:nx1,nz1) + 3.0_r2 * &
                   f(1:nx1,nz)   !forward/backward difference at boundaries
END DO

RETURN
END SUBROUTINE deriv_z

SUBROUTINE deriv_zz(f, fzz)
!Second z derivative
USE parameters
IMPLICIT NONE

REAL (r2), INTENT(IN)  :: f(0:nx,0:nz)
REAL (r2), INTENT(OUT) :: fzz(0:nx,0:nz)
INTEGER (i1)           :: j, k

DO k = 1, nz1
   fzz(1:nx1,k) = f(1:nx1,k+1) - 2.0_r2 * f(1:nx1,k) + f(1:nx1,k-1)
   fzz(1:nx1,0) = -f(1:nx1,3) + 4.0_r2 * f(1:nx1,2) - &
                    5.0_r2 * f(1:nx1,1) + 2.0_r2 * f(1:nx1,0) !forward/backward
   fzz(1:nx1,nz) = -f(1:nx1,nz-3) + 4.0_r2 * f(1:nx1,nz-2) - &  !difference at
                     5.0_r2 * f(1:nx1,nz1) + 2.0_r2 * f(1:nx1,nz)  !boundaries
END DO

RETURN
END SUBROUTINE deriv_zz

END MODULE derivs
