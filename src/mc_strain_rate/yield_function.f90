! Module holds the yield function equation and the derivatives of the yield function
module mod_yield_function
   use kind_precision_module, only: real_type => dp
   use mod_stress_invariants, only : calc_inc_driver_J3_invariant, Get_invariants, calc_J2_invariant
   use mod_stress_invar_deriv, only: calc_mean_stress_to_dSigma, calc_dq_to_dSigma, calc_dJ3_to_dSigma, calc_dtheta_to_dSigma
   use mod_voigt_functions   , only: calc_dev_stess
   implicit none

contains

   pure function calc_dF_to_dtheta(M_tc, p, theta) result(dfdtheta)
      real(kind = real_type), intent(in) :: M_tc, p, theta
      real(kind = real_type) :: dfdtheta

      ! ! Local variables
      real(kind = real_type), parameter :: pi=2.0*acos(0.0_real_type), &
         THREE_HALVES = 1.5_real_type, &
         ONE_QUARTER = 0.25_real_type, &
         TWO_TENTHS = 0.2_real_type

      ! !Get dF/dtheta
      dfdtheta = 0.45 * p * M_tc * ( ( cos(THREE_HALVES * theta + ONE_QUARTER * PI) )**TWO_TENTHS) &
         * sin( THREE_HALVES * theta + ONE_QUARTER * PI)

   end function calc_dF_to_dtheta

   pure subroutine Get_dF_to_dSigma(M_tc, eta_y, Sig, n_vec)
      !************************************************************************
      ! Returns the derivative of the yield function with respect to the		*
      ! stress tensor 														*
      ! n=dF/dSigma =dF/dp*dp/dSigma+ dF/dq*dq/dSigma +dF/dtheta*dtheta/dSigma*
      ! n is a (1X6) vector													*
      !************************************************************************
      implicit none
      !input
      real(kind = real_type), intent(in)  :: M_tc, eta_y, Sig(6)
      real(kind = real_type), intent(out) :: n_vec(6)

      ! Local variables
      real(kind = real_type):: p, q, theta, J2, J3, dJ3dsig(6), dfdtheta, &
         dpdsig(6), dqdsig(6), dev(6), dthetadSig(6)

      !Get the invariants
      call Get_invariants(Sig, p, q, theta)

      !Get dF/dp=eta_y and dF/dq=1
      !Get dF/dtheta
      dfdtheta= calc_dF_to_dtheta(M_tc, p, theta)

      !___________________________________________________________________________
      !1) Get dp/dSig=1/3 Imat
      dpdsig = calc_mean_stress_to_dSigma()

      dev= calc_dev_stess(Sig, p)

      !2) Get dq/dsig
      dqdSig = calc_dq_to_dSigma(dev, q)

      !3) Get dtheta/dSigma
      J2 = calc_J2_invariant(dev)

      J3 = calc_inc_driver_J3_invariant( dev )

      dJ3dsig = calc_dJ3_to_dSigma(dev)

      !Compute dtheta/dsig
      dthetadSig = calc_dtheta_to_dSigma(dJ3dsig, dev, J3, J2, theta)

      !Get n_vec=dF/dSig
      n_vec = ( eta_y * dpdsig ) + dqdSig + ( dfdtheta * dthetadSig) !n_vec=dF/dSig
   end subroutine Get_dF_to_dSigma

   subroutine Get_dF_to_dSigma_3(M_tc, eta_y, Sig, n_vec)
      !************************************************************************
      ! Returns the derivative of the yield function with respect to the		*
      ! stress tensor 														*
      ! n=dF/dSigma =dF/dp*dp/dSigma+ dF/dq*dq/dSigma +dF/dtheta*dtheta/dSigma*
      ! n is a (1X6) vector													*
      !************************************************************************
      implicit none
      !input
      real(kind = real_type), intent(in)  :: M_tc, eta_y, Sig(6)
      real(kind = real_type), intent(out) :: n_vec(6)

      ! Local variables
      real(kind = real_type):: p, q, theta, J2, J3, dJ3dsig(6), dfdtheta, &
         dpdsig(6), dqdsig(6), dev(6), dthetadSig(6)

      !Get the invariants
      call Get_invariants(Sig, p, q, theta)

      !Get dF/dp=eta_y and dF/dq=1
      !Get dF/dtheta
      dfdtheta= calc_dF_to_dtheta(M_tc, p, theta)

      !___________________________________________________________________________
      !1) Get dp/dSig=1/3 Imat
      dpdsig = calc_mean_stress_to_dSigma()

      dev= calc_dev_stess(Sig, p)

      !2) Get dq/dsig
      dqdSig = calc_dq_to_dSigma(dev, q)

      !3) Get dtheta/dSigma
      J2 = calc_J2_invariant(dev)

      J3 = calc_inc_driver_J3_invariant( dev )

      dJ3dsig = calc_dJ3_to_dSigma(dev)

      !Compute dtheta/dsig
      dthetadSig = calc_dtheta_to_dSigma(dJ3dsig, dev, J3, J2, theta)

      !Get n_vec=dF/dSig
      n_vec = ( eta_y * dpdsig ) + dqdSig + ( dfdtheta * dthetadSig) !n_vec=dF/dSig

   end subroutine Get_dF_to_dSigma_3

   pure subroutine YieldFunction(q, p, eta_y, F)
      !*********************************************************************
      ! Returns the value of the yield function evaluated at q, p , eta    *
      !																	 *
      !*********************************************************************
      implicit none
      real(kind = real_type), intent(in):: q, p, eta_y
      real(kind = real_type), intent(out):: F

      F=q+eta_y*p !sign is due to compression being negative in UMAT
   end subroutine YieldFunction

end module mod_yield_function
