! Module holds the derivatives of the stress invariants

module mod_stress_invar_deriv
   use stdlib_kinds, only: dp
   use mod_voigt_utils  , only : calc_voigt_square

   implicit none

contains

   pure function calc_dp_by_dsig() result(dp_by_dsig)
      real(dp) :: dp_by_dsig(6)

      ! Local variables
      integer :: i

      ! Zero all the values
      dp_by_dsig(:) = 0.0_dp

      ! Add the 1/3 to the first three elements
      do i = 1, 3
         dp_by_dsig(i) = 1.0_dp / 3.0_dp
      end do

   end function calc_dp_by_dsig

   pure function calc_dq_by_dsig(dev, q) result(dq_by_dsig)
      ! Calc dq/dSigma
      real(dp), intent(in) :: dev(6), q
      real(dp) :: dq_by_dsig(6)

      ! Local variables
      real(dp) :: dJ2_dsigma(6)

      ! Calc dJ2/dsigma
      dJ2_dSigma = calc_dJ2_by_dsig(dev)

      ! Calc dq/dSigma
      dq_by_dsig = 3.0_dp/ (2.0_dp * q) * dJ2_dsigma

   end function calc_dq_by_dsig

   pure function calc_dJ2_by_dsig(dev) result(dJ2_by_dsig)
      real(dp), intent(in) :: dev(6)
      real(dp) :: dJ2_by_dsig(6)

      dJ2_by_dsig = dev

      ! Double the shear terms
      dJ2_by_dsig(4:6) = 2.0_dp * dJ2_by_dsig(4:6)

   end function calc_dJ2_by_dsig

   pure function calc_dJ3_by_dsig(dev) result(dJ3_by_dsig)
      real(dp), intent(in) :: dev(6)
      real(dp) :: dJ3_by_dsig(6)

      ! Local variables
      real(kind=  dp) :: II(6), dev2(6), TrS2

      !Fill S.S
      dev2 = calc_voigt_square( dev )

      !Compute dJ3_by_dsig
      TrS2 = dev2(1) + dev2(2) + dev2(3)

      II=0.0d0!Identity tensor
      II(1)=1.0
      II(2)=1.0
      II(3)=1.0

      ! This is equaivalent to s^{2} - 2/3 J_{2} \matr{1}
      ! J_{2}(\matr{s}) = 2 * I_{1}(\matr{ s^{2} })
      ! See Appendix B. Invariant Notes Moore, Jonathan Thesis for more details
      dJ3_by_dsig = dev2 - ( TrS2*II / 3.0d0 )

      ! Need to double the shear terms because voigt notation is being used and therefore shear terms are linked together
      dJ3_by_dsig(4:6) = 2.0_dp * dJ3_by_dsig(4:6)

   end function calc_dJ3_by_dsig

   pure function calc_dJ3_by_dsig_full(stress) result(dJ3_by_dsig)
      ! This function calculates dJ3/dSigma assuming that the input voigt vector follows incremental driver convention
      ! As such the output follows incremental driver convention as well

      real(dp), intent(in) :: stress(6)
      real(dp) :: dJ3_by_dsig(6)

      ! Local variables
      real(kind=  dp) :: t(6)
      real(dp),parameter :: ONE_NINTH = 1.0_dp/9.0_dp, &
         ONE_THIRD = 1.0_dp/3.0_dp, &
         TWO       = 2.0_dp, &
         FOUR      = 4.0_dp

      ! Store the stress in a shorter name to make it easier to type
      t = stress

      ! Calc dJ3/dSigma_{11}
      dJ3_by_dsig(1) = ONE_NINTH * (TWO * t(1)**2 - t(2)**2 - t(3)**2 - &
         TWO * t(1) * t(2) - 2 * t(1) * t(3) + FOUR * t(2) * t(3)) + &
         ONE_THIRD * (t(4)**2 + t(5)**2 - TWO * t(6)**2)

      ! Calc dJ3/dSigma_{22}
      dJ3_by_dsig(2) = ONE_NINTH * (-t(1)**2 + TWO * t(2)**2 - t(3)**2 - &
         TWO * t(1) * t(2) + FOUR * t(1) * t(3) - TWO * t(2) * t(3)) + &
         ONE_THIRD * (t(4)**2 - TWO * t(5)**2 + t(6)**2)

      ! Calc dJ3/dSigma_{33}
      dJ3_by_dsig(3) = ONE_NINTH * (-t(1)**2 -t(2)**2 + TWO * t(3)**2 + &
         FOUR * t(1) * t(2) - TWO * t(1) * t(3) - TWO * t(2) * t(3)) + &
         ONE_THIRD * (-TWO * t(4)**2 + t(5)**2 + t(6)**2)

      ! Calc dJ3/dSigma_{12}
      dJ3_by_dsig(4) = ONE_THIRD * (TWO * t(1) * t(4) + TWO * t(4) * t(2) - FOUR * t(4) * t(3)) + TWO * t(5) * t(6)

      ! Calc dJ3/dSigma_{13}
      dJ3_by_dsig(5) = ONE_THIRD * (TWO * t(1) * t(5) - FOUR * t(5) * t(2) + TWO * t(5) * t(3)) + TWO * t(4) * t(6)

      ! Calc dJ3/dSigma_{23}
      dJ3_by_dsig(6) = ONE_THIRD * (-FOUR * t(1) * t(6) + TWO * t(2) * t(6) + TWO * t(6) * t(3)) + TWO * t(4) * t(5)

   end function calc_dJ3_by_dsig_full

   pure function calc_dtheta_by_dsig(dJ3_by_dsig, dev, J3, J2, theta) result(dtheta_by_dsig)
      real(dp), intent(in) :: dJ3_by_dsig(6), dev(6)
      real(dp), intent(in) :: J3, J2, theta
      real(dp) :: dtheta_by_dsig(6)

      ! Local variables
      real(dp) :: cos_term, outside_term, inside_term_1(6), dJ2_dSigma(6)
      real(dp), parameter :: tolerance = 1e-12_dp
      real(dp), parameter :: THREE = 3.0_dp, &
         TWO   = 2.0_dp, &
         ZERO  = 0.0_dp
      ! Calc cos(3 \theta)
      cos_term = cos(3 * theta)

      ! TODO: Turn the tolerance back on once I've checked that they match
      ! ! If cos term is zero( Trx compression or tension) set to tiny value
      ! if ( abs(cos_term) <= tolerance) then
      !    cos_term = tolerance
      ! end if

      ! Calc the fraction before the parenthesis
      outside_term = sqrt(THREE) / (TWO * cos_term * J2**1.5_dp)

      ! Calc the first term inside the parenthesis
      dJ2_dSigma = calc_dJ2_by_dsig(dev)

      inside_term_1 = THREE * J3 / (2 * J2) * dJ2_dSigma

      ! Calc the full term
      dtheta_by_dsig = outside_term * (inside_term_1 - dJ3_by_dsig)

   end function calc_dtheta_by_dsig

end module mod_stress_invar_deriv
