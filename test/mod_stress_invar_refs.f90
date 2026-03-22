! Reference implementations of stress invariant derivatives.
! These use the full Voigt-expanded form and exist solely for
! cross-checking the production implementations in tests.
! Do NOT use in production code — they are tightly coupled to the
! internal Voigt ordering [11,22,33,12,13,23].

module mod_stress_invar_refs
   use iso_fortran_env, only: dp => real64
   implicit none
   private
   public :: calc_dJ3_by_dsig_full

contains

   ! dJ3/dSigma computed directly from the full stress vector
   ! (Voigt order [11,22,33,12,13,23]) without first extracting
   ! the deviatoric tensor. Used as an independent cross-check in tests.
   pure function calc_dJ3_by_dsig_full(stress) result(dJ3_by_dsig)
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

end module mod_stress_invar_refs
