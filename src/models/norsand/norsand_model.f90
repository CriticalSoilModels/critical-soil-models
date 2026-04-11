!! NorSand — concrete implementation of `csm_model_t`.
!!
!! All constitutive mathematics lives in `mod_norsand_functions` (pure,
!! GPU-callable). The deferred procedure implementations here are thin
!! wrappers only — they forward to those pure functions with the model's
!! params and state unpacked.
!!
!! ### PROPS layout
!!
!! Required (indices 1–13):
!!
!! | Index | Symbol    | Description                                   | Units |
!! |-------|-----------|-----------------------------------------------|-------|
!! | 1     | G_0       | Reference shear modulus                       | kPa   |
!! | 2     | p_ref     | Reference mean stress (positive)              | kPa   |
!! | 3     | nG        | Shear modulus exponent                        | —     |
!! | 4     | nu        | Poisson's ratio                               | —     |
!! | 5     | e_o       | Initial void ratio                            | —     |
!! | 6     | Gamma     | CSL altitude at p = 1 kPa                     | —     |
!! | 7     | lambda_c  | CSL slope (ln scale)                          | —     |
!! | 8     | R         | OCR (used to initialise p_i)                  | —     |
!! | 9     | M_tc      | Critical friction ratio (triaxial compression)| —     |
!! | 10    | N         | Nova volumetric coupling coefficient          | —     |
!! | 11    | chi_tc    | Dilatancy coefficient (triaxial compression)  | —     |
!! | 12    | H_0       | Hardening modulus intercept                   | —     |
!! | 13    | H_y       | Hardening modulus slope                       | —     |
!!
!! Optional (indices 14–15; use defaults when omitted or ≤ 0):
!!
!! | Index | Symbol    | Description                                   | Units |
!! |-------|-----------|-----------------------------------------------|-------|
!! | 14    | ftol      | Yield surface tolerance (default 1e-8)        | —     |
!! | 15    | max_iters | Max integrator substep iterations (default 500)| —    |
!!
!! ### STATEV layout (15 entries)
!!
!! | Index | Symbol        | Description                              | Units |
!! |-------|---------------|------------------------------------------|-------|
!! | 1     | G             | Current shear modulus (from pre_step)    | kPa   |
!! | 2     | K             | Current bulk modulus (from pre_step)     | kPa   |
!! | 3     | p             | Current mean effective stress            | kPa   |
!! | 4     | e             | Void ratio                               | —     |
!! | 5     | psi           | State parameter ψ = e − e_c             | —     |
!! | 6     | chi_tce       | Current dilatancy coefficient            | —     |
!! | 7     | p_i           | Image mean stress (hardening variable)   | kPa   |
!! | 8     | M_i           | Image stress ratio                       | —     |
!! | 9     | switch_yield  | Yielding flag (0.0 = false, 1.0 = true)  | —     |
!! | 10–15 | eps_p         | Accumulated plastic strain (Voigt)       | —     |

module mod_norsand_model
   use mod_csm_kinds,       only: wp
   use mod_csm_model,       only: csm_model_t
   use mod_norsand_types,   only: norsand_params_t, norsand_state_t
   use mod_norsand_functions, only: norsand_yield_fn, norsand_flow_rule, norsand_plastic_potential, &
                                    norsand_update_hardening, norsand_elastic_stiffness, &
                                    norsand_hardening_modulus, norsand_pre_step
   implicit none
   private

   public :: norsand_model_t, norsand_from_props, norsand_load_state, norsand_save_state

   type, extends(csm_model_t) :: norsand_model_t
      !! NorSand critical state model.
      !!
      !! Extends `csm_model_t`. Parameters are fixed after construction via
      !! `norsand_from_props`; state evolves during `update_hardening` and
      !! `pre_step` calls.
      type(norsand_params_t) :: params  !! Fixed model parameters
      type(norsand_state_t)  :: state   !! Evolving internal state
   contains
      procedure :: yield_fn          => norsand_yield_fn_method
      procedure :: flow_rule         => norsand_flow_rule_method
      procedure :: plastic_potential => norsand_plastic_potential_method
      procedure :: update_hardening  => norsand_update_hardening_method
      procedure :: elastic_stiffness => norsand_elastic_stiffness_method
      procedure :: snapshot          => norsand_snapshot
      procedure :: restore           => norsand_restore
      procedure :: hardening_modulus => norsand_hardening_modulus_method
      procedure :: pre_step          => norsand_pre_step_method
   end type norsand_model_t

contains

   ! ---------------------------------------------------------------------------
   ! UMAT boundary helpers — only place PROPS/STATEV index magic lives
   ! ---------------------------------------------------------------------------

   function norsand_from_props(props) result(model)
      !! Construct a `norsand_model_t` from a flat PROPS array (UMAT convention).
      real(wp), intent(in) :: props(:)  !! PROPS array (see module doc for layout)
      type(norsand_model_t) :: model

      model%params%G_0      = props(1)
      model%params%p_ref    = props(2)
      model%params%nG       = props(3)
      model%params%nu       = props(4)
      model%params%e_o      = props(5)
      model%params%Gamma    = props(6)
      model%params%lambda_c = props(7)
      model%params%R        = props(8)
      model%params%M_tc     = props(9)
      model%params%N        = props(10)
      model%params%chi_tc   = props(11)
      model%params%H_0      = props(12)
      model%params%H_y      = props(13)
   end function norsand_from_props

   subroutine norsand_load_state(model, statev)
      !! Unpack a STATEV array into the model's internal state fields.
      type(norsand_model_t), intent(inout) :: model
      real(wp),              intent(in)    :: statev(:)  !! STATEV array (see module doc for layout)

      model%state%G            = statev(1)
      model%state%K            = statev(2)
      model%state%p            = statev(3)
      model%state%e            = statev(4)
      model%state%psi          = statev(5)
      model%state%chi_tce      = statev(6)
      model%state%p_i          = statev(7)
      model%state%M_i          = statev(8)
      model%state%switch_yield = statev(9) > 0.5_wp
      model%state%eps_p        = statev(10:15)
   end subroutine norsand_load_state

   subroutine norsand_save_state(model, statev)
      !! Pack the model's internal state back into a STATEV array.
      type(norsand_model_t), intent(in)    :: model
      real(wp),              intent(inout) :: statev(:)  !! STATEV array (see module doc for layout)

      statev(1)     = model%state%G
      statev(2)     = model%state%K
      statev(3)     = model%state%p
      statev(4)     = model%state%e
      statev(5)     = model%state%psi
      statev(6)     = model%state%chi_tce
      statev(7)     = model%state%p_i
      statev(8)     = model%state%M_i
      statev(9)     = merge(1.0_wp, 0.0_wp, model%state%switch_yield)
      statev(10:15) = model%state%eps_p
   end subroutine norsand_save_state

   ! ---------------------------------------------------------------------------
   ! Thin wrappers around pure functions
   ! ---------------------------------------------------------------------------

   function norsand_yield_fn_method(self, sig) result(F)
      class(norsand_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: F
      F = norsand_yield_fn(self%params, self%state, sig)
   end function norsand_yield_fn_method

   function norsand_flow_rule_method(self, sig) result(dF_by_dsig)
      class(norsand_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: dF_by_dsig(6)
      dF_by_dsig = norsand_flow_rule(self%params, self%state, sig)
   end function norsand_flow_rule_method

   function norsand_plastic_potential_method(self, sig) result(dG_by_dsig)
      class(norsand_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: dG_by_dsig(6)
      dG_by_dsig = norsand_plastic_potential(self%params, self%state, sig)
   end function norsand_plastic_potential_method

   subroutine norsand_update_hardening_method(self, deps_p)
      class(norsand_model_t), intent(inout) :: self
      real(wp), intent(in) :: deps_p(6)
      call norsand_update_hardening(self%params, self%state, deps_p)
   end subroutine norsand_update_hardening_method

   function norsand_elastic_stiffness_method(self) result(stiff_e)
      class(norsand_model_t), intent(in) :: self
      real(wp) :: stiff_e(6,6)
      stiff_e = norsand_elastic_stiffness(self%params, self%state)
   end function norsand_elastic_stiffness_method

   function norsand_hardening_modulus_method(self, sig, dg_by_dsig) result(H)
      class(norsand_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6), dg_by_dsig(6)
      real(wp) :: H
      H = norsand_hardening_modulus(self%params, self%state, sig, dg_by_dsig)
   end function norsand_hardening_modulus_method

   subroutine norsand_pre_step_method(self, sig)
      class(norsand_model_t), intent(inout) :: self
      real(wp), intent(in) :: sig(6)
      call norsand_pre_step(self%params, self%state, sig)
   end subroutine norsand_pre_step_method

   ! ---------------------------------------------------------------------------
   ! snapshot / restore — for integrator rollback
   ! ---------------------------------------------------------------------------

   subroutine norsand_snapshot(self, saved)
      class(norsand_model_t), intent(in)    :: self
      real(wp), allocatable,  intent(out)   :: saved(:)

      allocate(saved(15))
      saved(1)     = self%state%G
      saved(2)     = self%state%K
      saved(3)     = self%state%p
      saved(4)     = self%state%e
      saved(5)     = self%state%psi
      saved(6)     = self%state%chi_tce
      saved(7)     = self%state%p_i
      saved(8)     = self%state%M_i
      saved(9)     = merge(1.0_wp, 0.0_wp, self%state%switch_yield)
      saved(10:15) = self%state%eps_p
   end subroutine norsand_snapshot

   subroutine norsand_restore(self, saved)
      class(norsand_model_t), intent(inout) :: self
      real(wp),               intent(in)    :: saved(:)

      self%state%G            = saved(1)
      self%state%K            = saved(2)
      self%state%p            = saved(3)
      self%state%e            = saved(4)
      self%state%psi          = saved(5)
      self%state%chi_tce      = saved(6)
      self%state%p_i          = saved(7)
      self%state%M_i          = saved(8)
      self%state%switch_yield = saved(9) > 0.5_wp
      self%state%eps_p        = saved(10:15)
   end subroutine norsand_restore

end module mod_norsand_model
