MODULE ccf
!use parameters
implicit none

contains

SUBROUTINE get_A(t, A, A_)
use parameters
implicit none
double precision, intent(out) :: A, A_
double precision, intent(in) :: t
double precision :: t_

t_ = t - dt

A = (1d0 / (1d0 - eta**2)) * ((Re2 + Re2_mod * cos(om2 * t)) - &
    eta * (Re1 + Re1_mod * cos(om1 * t)))

A_ = (1d0 / (1d0 - eta**2)) * ((Re2 + Re2_mod * cos(om2 * t_)) - &
    eta * (Re1 + Re1_mod * cos(om1 * t_)))

!A = -171.4285714d0
!A_ = -171.4285714d0

return
END SUBROUTINE get_A

SUBROUTINE get_B(t, B, B_)
use parameters
implicit none
double precision, intent(out) :: B, B_
double precision, intent(in) :: t
double precision :: t_

t_ = t - dt

B = (eta / (1d0 - eta**2)) * ((Re1 + Re1_mod * cos(om1 * t)) - &
    eta * (Re2 + Re2_mod * cos(om2 * t)))

B_ = (eta / (1d0 - eta**2)) * ((Re1 + Re1_mod * cos(om1 * t_)) - &
    eta * (Re2 + Re2_mod * cos(om2 * t_)))

!B = 171.4285714d0
!B_ = 171.4285714d0

return
END SUBROUTINE get_B

SUBROUTINE get_vc(Avel, A_vel, Bvel, B_vel, vc, vc_, s)
use parameters
implicit none
double precision, intent(out) :: vc(0:nx), vc_(0:nx)
double precision, intent(in) :: Avel, A_vel, Bvel, B_vel, s(0:nx)

vc(:) = Avel * s(:) + Bvel / s(:)
vc_(:) = A_vel * s(:) + B_vel / s(:)

!print*,A(t)
!print*,B(t)
!print*,vc
return
END SUBROUTINE get_vc

SUBROUTINE get_F(t, s, F)
use parameters
implicit none
double precision, intent(in) :: t, s(0:nx)
double precision, intent(out) :: F(0:nx)

F(:) = (1d0 / (1d0 - eta**2)) * &
       ((-Re1_mod * (eta * s(:) - eta / s(:)) * &
       (cos(om1 * (t + dt)) - cos(om1 * t))) + &
       (Re2_mod * (s(:) - (eta**2 / s(:))) * &
       (cos(om2 * (t + dt)) - cos(om2 * t))))

return
END SUBROUTINE get_F

!FUNCTION F(t, s_el)
!use parameters
!implicit none
!double precision :: F
!double precision, intent(in) :: t, s_el !s(0:nx), t

!F = 0d0 !(1d0 / (1d0 - eta**2)) * &
!    !((-Re1_mod * (eta * s_el - eta / s_el) * &
!    !(cos(om1 * (t + dt)) - cos(om1 * t))) + &
!    !(Re2_mod * (s_el - (eta**2 / s_el)) * &
!    !(cos(om2 * (t + dt)) - cos(om2 * t))))

!return
!END FUNCTION F

END MODULE ccf
