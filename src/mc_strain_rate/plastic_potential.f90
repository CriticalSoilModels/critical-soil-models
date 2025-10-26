! Module for holding the function that calcualtes the plastic potential function value and and the derivatives

module mod_plastic_potential
   use kind_precision_module , only: real_type => dp
   use mod_stress_invariants , only: Get_invariants
   use mod_stress_invar_deriv, only: calc_mean_stress_to_dSigma, calc_dq_to_dSigma
   use mod_voigt_functions   , only: calc_dev_stess
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
      real(kind = real_type), intent(in)  :: D, Sig(6)
      !output
      real(kind = real_type), intent(out) :: m_vec(6)

      !local variables
      real(kind = real_type):: p, q, theta, dpdsig(6), dev(6), dqdSig(6)
      
      ! Get the invariants
      call Get_invariants(Sig, p, q, theta)

      !Get dP/dp=-D and dF/dq=1
      !___________________________________________________________________________
      !1) Get dp/dSig=1/3 Imat
      dpdsig=0.0d0
            
      dpdsig = calc_mean_stress_to_dSigma()

      !2) Get dq/dsig
      dev = calc_dev_stess(Sig, p)

      dqdSig = calc_dq_to_dSigma(dev, q)

      !__________________________________________________________________
      !Get m_vec=dP/dSig]
      ! Dilatancy is negative there
      m_vec=(-D*dpdsig)+dqdSig !m_vec=dP/dSig

   end subroutine Get_dP_to_dSigma
end module mod_plastic_potential
