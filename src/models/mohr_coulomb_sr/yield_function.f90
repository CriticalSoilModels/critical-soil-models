! Module holds the yield function equation and the derivatives of the yield function
module mod_yield_function
   use mod_csm_kinds, only: wp
   use mod_stress_invariants, only : calc_J3_inv, calc_sig_inv, calc_J2_inv
   use mod_stress_invar_deriv, only: calc_dp_by_dsig, calc_dq_by_dsig, calc_dJ3_by_dsig, calc_dlode_angle_by_dsig
   use mod_voigt_utils   , only: calc_dev_stress
   implicit none

contains

   pure function calc_dF_by_dtheta(M_tc, p, theta) result(dF_by_dtheta)
      real(wp), intent(in) :: M_tc, p, theta
      real(wp) :: dF_by_dtheta

      real(wp), parameter :: pi = 2.0_wp * acos(0.0_wp), &
         THREE_HALVES = 1.5_wp, &
         ONE_QUARTER  = 0.25_wp, &
         TWO_TENTHS   = 0.2_wp

      dF_by_dtheta = 0.45_wp * p * M_tc * ( ( cos(THREE_HALVES * theta + ONE_QUARTER * pi) )**TWO_TENTHS) &
         * sin( THREE_HALVES * theta + ONE_QUARTER * pi)

   end function calc_dF_by_dtheta

   pure subroutine calc_dF_by_dsig(M_tc, eta_y, sig, n_vec)
      !! Returns dF/dsig = dF/dp * dp/dsig + dF/dq * dq/dsig + dF/dtheta * dtheta/dsig
      !! n_vec is a 6-vector in Voigt notation
      implicit none
      real(wp), intent(in)  :: M_tc, eta_y, sig(6)
      real(wp), intent(out) :: n_vec(6)

      real(wp) :: p, q, lode_angle, J2, J3
      real(wp) :: dJ3_by_dsig(6), dF_by_dtheta, dp_by_dsig(6), dq_by_dsig(6)
      real(wp) :: dev(6), dlode_angle_by_dsig(6)

      call calc_sig_inv(sig, p, q, lode_angle)

      dF_by_dtheta = calc_dF_by_dtheta(M_tc, p, lode_angle)

      dp_by_dsig = calc_dp_by_dsig()

      dev = calc_dev_stress(sig, p)

      dq_by_dsig = calc_dq_by_dsig(dev, q)

      J2 = calc_J2_inv(dev)
      J3 = calc_J3_inv(dev)

      dJ3_by_dsig = calc_dJ3_by_dsig(dev)

      dlode_angle_by_dsig = calc_dlode_angle_by_dsig(dJ3_by_dsig, dev, J3, J2, lode_angle)

      n_vec = (eta_y * dp_by_dsig) + dq_by_dsig + (dF_by_dtheta * dlode_angle_by_dsig)

   end subroutine calc_dF_by_dsig


   pure subroutine calc_yield_function(q, p, eta_y, F)
      !! Returns the value of the yield function F = q + eta_y * p
      !! Sign convention: compression is negative in UMAT
      implicit none
      real(wp), intent(in)  :: q, p, eta_y
      real(wp), intent(out) :: F

      F = q + eta_y * p

   end subroutine calc_yield_function

end module mod_yield_function
