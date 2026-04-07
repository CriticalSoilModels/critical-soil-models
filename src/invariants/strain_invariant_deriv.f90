!! Derivatives of strain invariants with respect to the strain vector.

module mod_strain_invar_deriv
   use mod_csm_kinds, only: wp
   use mod_strain_invariants, only: calc_dev_strain, calc_eps_vol_inv
   implicit none

contains

   pure function calc_deps_q_by_deps(eps_q, eps) result(deps_q_by_deps)
      !! Derivative of the deviatoric strain invariant εq with respect to strain: ∂εq/∂ε.
      !!
      !! Returns a zero vector when εq = 0 to avoid division by zero at the origin.
      real(wp), intent(in) :: eps_q     !! Deviatoric strain invariant εq
      real(wp), intent(in) :: eps(6)    !! Strain vector (Voigt)
      real(wp) :: deps_q_by_deps(6)
      real(wp) :: evol, dev(6)

      evol = calc_eps_vol_inv(eps)
      dev  = calc_dev_strain(eps, evol)

      if (eps_q > 0.0_wp) then
         deps_q_by_deps = (2.0_wp / (3.0_wp * eps_q)) * dev
      else
         deps_q_by_deps = 0.0_wp
      end if
   end function calc_deps_q_by_deps

end module mod_strain_invar_deriv
