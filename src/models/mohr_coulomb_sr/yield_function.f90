! Module holds the yield function equation and the derivatives of the yield function
module mod_yield_function
   use stdlib_kinds, only: dp
   use mod_stress_invariants, only : calc_J3, calc_stress_invariants, calc_J2
   use mod_stress_invar_deriv, only: calc_dp_by_dsig, calc_dq_by_dsig, calc_dJ3_by_dsig, calc_dlode_angle_by_dsig
   use mod_voigt_utils   , only: calc_dev_stress
   implicit none

contains

   pure function calc_dF_to_dtheta(M_tc, p, theta) result(dfdtheta)
      real(dp), intent(in) :: M_tc, p, theta
      real(dp) :: dfdtheta

      ! ! Local variables
      real(dp), parameter :: pi=2.0*acos(0.0_dp), &
         THREE_HALVES = 1.5_dp, &
         ONE_QUARTER = 0.25_dp, &
         TWO_TENTHS = 0.2_dp

      ! !Get dF/dtheta
      dfdtheta = 0.45 * p * M_tc * ( ( cos(THREE_HALVES * theta + ONE_QUARTER * PI) )**TWO_TENTHS) &
         * sin( THREE_HALVES * theta + ONE_QUARTER * PI)

   end function calc_dF_to_dtheta

   pure subroutine Get_dF_to_dSigma(M_tc, eta_y, Sig, n_vec)
      !************************************************************************
      ! Returns the derivative of the yield function with respect to the		*
      ! stress tensor 														*
      ! n=dF/dSigma =dF/dp*dp/dSigma+ dF/dq*dq/dSigma +dF/dlode_angle*dlode_angle/dSigma*
      ! n is a (1X6) vector													*
      !************************************************************************
      implicit none
      !input
      real(dp), intent(in)  :: M_tc, eta_y, Sig(6)
      real(dp), intent(out) :: n_vec(6)

      ! Local variables
      real(dp):: p, q, lode_angle, J2, J3, dJ3dsig(6), dfdtheta, &
         dpdsig(6), dqdsig(6), dev(6), dlode_angle_dSig(6)

      !Get the invariants
      call calc_stress_invariants(Sig, p, q, lode_angle)

      !Get dF/dp=eta_y and dF/dq=1
      !Get dF/dtheta
      dfdtheta= calc_dF_to_dtheta(M_tc, p, lode_angle)

      !___________________________________________________________________________
      !1) Get dp/dSig=1/3 Imat
      dpdsig = calc_dp_by_dsig()

      dev= calc_dev_stress(Sig, p)

      !2) Get dq/dsig
      dqdSig = calc_dq_by_dsig(dev, q)

      !3) Get dlode_angle/dSigma
      J2 = calc_J2(dev)

      J3 = calc_J3( dev )

      dJ3dsig = calc_dJ3_by_dsig(dev)

      !Compute dlode_angle/dsig
      dlode_angle_dSig = calc_dlode_angle_by_dsig(dJ3dsig, dev, J3, J2, lode_angle)

      !Get n_vec=dF/dSig
      n_vec = ( eta_y * dpdsig ) + dqdSig + ( dfdtheta * dlode_angle_dSig) !n_vec=dF/dSig
   end subroutine Get_dF_to_dSigma

   subroutine Get_dF_to_dSigma_3(M_tc, eta_y, Sig, n_vec)
      !************************************************************************
      ! Returns the derivative of the yield function with respect to the		*
      ! stress tensor 														*
      ! n=dF/dSigma =dF/dp*dp/dSigma+ dF/dq*dq/dSigma +dF/dlode_angle*dlode_angle/dSigma*
      ! n is a (1X6) vector													*
      !************************************************************************
      implicit none
      !input
      real(dp), intent(in)  :: M_tc, eta_y, Sig(6)
      real(dp), intent(out) :: n_vec(6)

      ! Local variables
      real(dp):: p, q, lode_angle, J2, J3, dJ3dsig(6), dfdtheta, &
         dpdsig(6), dqdsig(6), dev(6), dlode_angle_dSig(6)

      !Get the invariants
      call calc_stress_invariants(Sig, p, q, lode_angle)

      !Get dF/dp=eta_y and dF/dq=1
      !Get dF/dtheta
      dfdtheta= calc_dF_to_dtheta(M_tc, p, lode_angle)

      !___________________________________________________________________________
      !1) Get dp/dSig=1/3 Imat
      dpdsig = calc_dp_by_dsig()

      dev= calc_dev_stress(Sig, p)

      !2) Get dq/dsig
      dqdSig = calc_dq_by_dsig(dev, q)

      !3) Get dlode_angle/dSigma
      J2 = calc_J2(dev)

      J3 = calc_J3( dev )

      dJ3dsig = calc_dJ3_by_dsig(dev)

      !Compute dlode_angle/dsig
      dlode_angle_dSig = calc_dlode_angle_by_dsig(dJ3dsig, dev, J3, J2, lode_angle)

      !Get n_vec=dF/dSig
      n_vec = ( eta_y * dpdsig ) + dqdSig + ( dfdtheta * dlode_angle_dSig) !n_vec=dF/dSig

   end subroutine Get_dF_to_dSigma_3

   pure subroutine YieldFunction(q, p, eta_y, F)
      !*********************************************************************
      ! Returns the value of the yield function evaluated at q, p , eta    *
      !																	 *
      !*********************************************************************
      implicit none
      real(dp), intent(in):: q, p, eta_y
      real(dp), intent(out):: F

      F=q+eta_y*p !sign is due to compression being negative in UMAT
   end subroutine YieldFunction

end module mod_yield_function
