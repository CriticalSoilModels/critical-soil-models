! Pure constitutive functions for the linear elastic model.
!
! All procedures are pure and take (params, state, ...) as explicit arguments —
! no class(), no dynamic dispatch. Safe to call from GPU kernels.
!
! POD types:
!   le_params_t  — fixed material parameters (G, nu)
!   le_state_t   — internal state (empty: linear elastic has none)

module mod_le_functions
   use mod_csm_kinds,     only: wp
   use mod_le_types,      only: le_params_t, le_state_t
   use mod_elastic_utils, only: calc_stiffness_GK, calc_K_from_G_nu
   implicit none
   private

   public :: le_params_t, le_state_t
   public :: le_yield_fn
   public :: le_flow_rule
   public :: le_plastic_potential
   public :: le_update_hardening
   public :: le_elastic_stiffness

contains

   pure function le_yield_fn(params, state, sig) result(F)
      !! Always returns -1: linear elastic never yields.
      type(le_params_t), intent(in) :: params
      type(le_state_t),  intent(in) :: state
      real(wp),          intent(in) :: sig(6)
      real(wp) :: F
      F = -1.0_wp
   end function le_yield_fn

   pure function le_flow_rule(params, state, sig) result(dF_by_dsig)
      !! Never called — yield function is always negative.
      type(le_params_t), intent(in) :: params
      type(le_state_t),  intent(in) :: state
      real(wp),          intent(in) :: sig(6)
      real(wp) :: dF_by_dsig(6)
      dF_by_dsig = 0.0_wp
   end function le_flow_rule

   pure function le_plastic_potential(params, state, sig) result(dG_by_dsig)
      !! Never called — yield function is always negative.
      type(le_params_t), intent(in) :: params
      type(le_state_t),  intent(in) :: state
      real(wp),          intent(in) :: sig(6)
      real(wp) :: dG_by_dsig(6)
      dG_by_dsig = 0.0_wp
   end function le_plastic_potential

   pure subroutine le_update_hardening(params, state, deps_p)
      !! No-op — no internal state to update.
      type(le_params_t), intent(in)    :: params
      type(le_state_t),  intent(inout) :: state
      real(wp),          intent(in)    :: deps_p(6)
   end subroutine le_update_hardening

   pure function le_elastic_stiffness(params, state) result(stiff_e)
      !! Isotropic linear elastic stiffness in Voigt notation [11,22,33,12,13,23].
      type(le_params_t), intent(in) :: params
      type(le_state_t),  intent(in) :: state
      real(wp) :: stiff_e(6,6)
      stiff_e = calc_stiffness_GK(params%G, calc_K_from_G_nu(params%G, params%nu))
   end function le_elastic_stiffness

end module mod_le_functions
