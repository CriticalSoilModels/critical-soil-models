! Module for holding the derivatives of the strain invariants

module mod_strain_invar_deriv
   use stdlib_kinds, only: dp
   use mod_strain_invariants, only: calc_dev_strain, calc_eps_vol_invariant
   implicit none

contains
   subroutine calc_deps_q_by_deps(Epsq, Eps, deps_q_by_deps)
      !************************************************************************
      ! Returns the derivative of the deviatoric strain with respect to the   *
      ! deviatoric strain	tensor						     					*
      ! dEqdEpsq is a (1X6) vector											*
      !************************************************************************
      implicit none
      !input
      real(dp), intent(in):: Epsq, Eps(6)
      !output
      real(dp), intent(out):: deps_q_by_deps(6)
      !local variables
      real(dp):: evol, dev(6)

      evol=calc_eps_vol_invariant(Eps)

      dev=calc_dev_strain(Eps, evol)

      if (Epsq>0.0d0) then !in case of zero plastic strain
         deps_q_by_deps=(2.0/(3.0*Epsq))*dev
      else
         deps_q_by_deps=0.0d0
      endif
   end subroutine calc_deps_q_by_deps

end module mod_strain_invar_deriv
