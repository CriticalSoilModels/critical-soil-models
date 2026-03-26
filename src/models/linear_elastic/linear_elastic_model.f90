! Linear elastic model — simplest concrete implementation of csm_model_t.
! No yield surface, no state variables. Always elastic.
! Used to validate the abstract base + integrator stack end-to-end.
!
! All constitutive math lives in mod_le_functions (pure, GPU-callable).
! The deferred procedure implementations here are thin wrappers only.
!
! PROPS layout:
!   1: G   shear modulus [kPa]
!   2: nu  Poisson's ratio [-]

module mod_linear_elastic_model
   use mod_csm_kinds,   only: wp
   use mod_csm_model,   only: csm_model_t
   use mod_le_functions, only: le_params_t, le_state_t, &
                                le_yield_fn, le_flow_rule, le_plastic_potential, &
                                le_update_hardening, le_elastic_stiffness
   implicit none
   private
   public :: linear_elastic_model_t, le_from_props

   type, extends(csm_model_t) :: linear_elastic_model_t
      type(le_params_t) :: params
      type(le_state_t)  :: state
   contains
      procedure :: yield_fn          => le_yield_fn_method
      procedure :: flow_rule         => le_flow_rule_method
      procedure :: plastic_potential => le_plastic_potential_method
      procedure :: update_hardening  => le_update_hardening_method
      procedure :: elastic_stiffness => le_elastic_stiffness_method
      procedure :: snapshot          => le_snapshot
      procedure :: restore           => le_restore
   end type linear_elastic_model_t

contains

   function le_from_props(props) result(model)
      real(wp), intent(in) :: props(:)
      type(linear_elastic_model_t) :: model
      model%params%G  = props(1)
      model%params%nu = props(2)
   end function le_from_props

   ! --- Thin wrappers around pure functions ---

   function le_yield_fn_method(self, sig) result(F)
      class(linear_elastic_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: F
      F = le_yield_fn(self%params, self%state, sig)
   end function le_yield_fn_method

   function le_flow_rule_method(self, sig) result(dF_by_dsig)
      class(linear_elastic_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: dF_by_dsig(6)
      dF_by_dsig = le_flow_rule(self%params, self%state, sig)
   end function le_flow_rule_method

   function le_plastic_potential_method(self, sig) result(dG_by_dsig)
      class(linear_elastic_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: dG_by_dsig(6)
      dG_by_dsig = le_plastic_potential(self%params, self%state, sig)
   end function le_plastic_potential_method

   subroutine le_update_hardening_method(self, deps_p)
      class(linear_elastic_model_t), intent(inout) :: self
      real(wp), intent(in) :: deps_p(6)
      call le_update_hardening(self%params, self%state, deps_p)
   end subroutine le_update_hardening_method

   function le_elastic_stiffness_method(self) result(stiff_e)
      class(linear_elastic_model_t), intent(in) :: self
      real(wp) :: stiff_e(6,6)
      stiff_e = le_elastic_stiffness(self%params, self%state)
   end function le_elastic_stiffness_method

   ! --- OOP-only: snapshot/restore for the generic CPU integrator ---

   subroutine le_snapshot(self, saved)
      !! No state to save — allocates an empty array for the integrator
      class(linear_elastic_model_t), intent(in)  :: self
      real(wp), allocatable,         intent(out) :: saved(:)
      allocate(saved(0))
   end subroutine le_snapshot

   subroutine le_restore(self, saved)
      !! No state to restore — no-op
      class(linear_elastic_model_t), intent(inout) :: self
      real(wp),                      intent(in)    :: saved(:)
   end subroutine le_restore

end module mod_linear_elastic_model
