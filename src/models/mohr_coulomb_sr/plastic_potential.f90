! Module for holding the function that calculates the plastic potential derivatives

module mod_plastic_potential
   use stdlib_kinds, only: dp
   use mod_stress_invariants , only: calc_sig_invariants
   use mod_stress_invar_deriv, only: calc_dp_by_dsig, calc_dq_by_dsig
   use mod_voigt_utils   , only: calc_dev_stress
   implicit none

contains

   subroutine calc_dg_plas_by_dsig(dilation, sig, m_vec)
      !! Returns dP/dsig = dP/dp * dp/dsig + dP/dq * dq/dsig
      !! where dP/dp = -dilation and dP/dq = 1
      !! m_vec is a 6-vector in Voigt notation
      implicit none
      real(dp), intent(in)  :: dilation, sig(6)
      real(dp), intent(out) :: m_vec(6)

      real(dp) :: p, q, lode_angle, dp_by_dsig(6), dev(6), dq_by_dsig(6)

      call calc_sig_invariants(sig, p, q, lode_angle)

      dp_by_dsig = calc_dp_by_dsig()

      dev = calc_dev_stress(sig, p)

      dq_by_dsig = calc_dq_by_dsig(dev, q)

      m_vec = (-dilation * dp_by_dsig) + dq_by_dsig

   end subroutine calc_dg_plas_by_dsig

end module mod_plastic_potential
