! Module for holding the function that calcualtes the plastic potential function value and and the derivatives

module mod_plastic_potential
   use stdlib_kinds, only: dp
   use mod_stress_invariants , only: calc_stress_invariants
   use mod_stress_invar_deriv, only: calc_dp_by_dsig, calc_dq_by_dsig
   use mod_voigt_utils   , only: calc_dev_stress
   implicit none

contains
   subroutine Get_dP_to_dSigma(D, Sig, m_vec)
      !************************************************************************
      ! Returns the derivative of the plastic potential function with respect *
      ! to the stress tensor													*
      ! m=dP/dSigma =dP/dp*dp/dSigma+ dP/dq*dq/dSigma							*
      ! m is a (1X6) vector													*
      !************************************************************************
      implicit none
      !input
      real(dp), intent(in)  :: D, Sig(6)
      !output
      real(dp), intent(out) :: m_vec(6)

      !local variables
      real(dp):: p, q, lode_angle, dpdsig(6), dev(6), dqdSig(6)

      ! Get the invariants
      call calc_stress_invariants(Sig, p, q, lode_angle)

      !Get dP/dp=-D and dF/dq=1
      !___________________________________________________________________________
      !1) Get dp/dSig=1/3 Imat
      dpdsig=0.0d0

      dpdsig = calc_dp_by_dsig()

      !2) Get dq/dsig
      dev = calc_dev_stress(Sig, p)

      dqdSig = calc_dq_by_dsig(dev, q)

      !__________________________________________________________________
      !Get m_vec=dP/dSig]
      ! Dilatancy is negative there
      m_vec=(-D*dpdsig)+dqdSig !m_vec=dP/dSig

   end subroutine Get_dP_to_dSigma
end module mod_plastic_potential
