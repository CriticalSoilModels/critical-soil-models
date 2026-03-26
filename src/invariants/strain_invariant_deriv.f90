! Module for holding the derivatives of the strain invariants

module mod_strain_invar_deriv
   use mod_csm_kinds, only: wp
   use mod_strain_invariants, only: calc_dev_strain, calc_eps_vol_inv
   implicit none

contains
   pure function calc_deps_q_by_deps(eps_q, eps) result(deps_q_by_deps)
      ! Returns dεq/dε — derivative of deviatoric strain invariant w.r.t. strain.
      real(wp), intent(in) :: eps_q, eps(6)
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
