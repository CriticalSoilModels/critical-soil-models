! Module holds the yield function equation and the derivatives of the yield function
module mod_yield_function
   use stdlib_kinds, only: dp
   use mod_stress_invariants, only : calc_J3, calc_sig_invariants, calc_J2
   use mod_stress_invar_deriv, only: calc_dp_by_dsig, calc_dq_by_dsig, calc_dJ3_by_dsig, calc_dlode_angle_by_dsig
   use mod_voigt_utils   , only: calc_dev_stress
   implicit none

contains

   pure function calc_dF_by_dtheta(M_tc, p, theta) result(dF_by_dtheta)
      real(dp), intent(in) :: M_tc, p, theta
      real(dp) :: dF_by_dtheta

      real(dp), parameter :: pi = 2.0_dp * acos(0.0_dp), &
         THREE_HALVES = 1.5_dp, &
         ONE_QUARTER  = 0.25_dp, &
         TWO_TENTHS   = 0.2_dp

      dF_by_dtheta = 0.45_dp * p * M_tc * ( ( cos(THREE_HALVES * theta + ONE_QUARTER * pi) )**TWO_TENTHS) &
         * sin( THREE_HALVES * theta + ONE_QUARTER * pi)

   end function calc_dF_by_dtheta

   pure subroutine calc_dF_by_dsig(M_tc, eta_y, sig, n_vec)
      !! Returns dF/dsig = dF/dp * dp/dsig + dF/dq * dq/dsig + dF/dtheta * dtheta/dsig
      !! n_vec is a 6-vector in Voigt notation
      implicit none
      real(dp), intent(in)  :: M_tc, eta_y, sig(6)
      real(dp), intent(out) :: n_vec(6)

      real(dp) :: p, q, lode_angle, J2, J3
      real(dp) :: dJ3_by_dsig(6), dF_by_dtheta, dp_by_dsig(6), dq_by_dsig(6)
      real(dp) :: dev(6), dlode_angle_by_dsig(6)

      call calc_sig_invariants(sig, p, q, lode_angle)

      dF_by_dtheta = calc_dF_by_dtheta(M_tc, p, lode_angle)

      dp_by_dsig = calc_dp_by_dsig()

      dev = calc_dev_stress(sig, p)

      dq_by_dsig = calc_dq_by_dsig(dev, q)

      J2 = calc_J2(dev)
      J3 = calc_J3(dev)

      dJ3_by_dsig = calc_dJ3_by_dsig(dev)

      dlode_angle_by_dsig = calc_dlode_angle_by_dsig(dJ3_by_dsig, dev, J3, J2, lode_angle)

      n_vec = (eta_y * dp_by_dsig) + dq_by_dsig + (dF_by_dtheta * dlode_angle_by_dsig)

   end subroutine calc_dF_by_dsig


   pure subroutine calc_yield_function(q, p, eta_y, F)
      !! Returns the value of the yield function F = q + eta_y * p
      !! Sign convention: compression is negative in UMAT
      implicit none
      real(dp), intent(in)  :: q, p, eta_y
      real(dp), intent(out) :: F

      F = q + eta_y * p

   end subroutine calc_yield_function

end module mod_yield_function
