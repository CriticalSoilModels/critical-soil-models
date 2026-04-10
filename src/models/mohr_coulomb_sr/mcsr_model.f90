!! Mohr-Coulomb Strain Rate — concrete implementation of `csm_model_t`.
!!
!! `mcsr_model_t` extends the abstract base type and forwards all deferred
!! procedure calls to the pure functions in `mod_mcsr_functions`.
!!
!! Rate-dependent state (G, K, eta_y, dilation, M_lode) must be updated by
!! calling `mcsr_update_rate_state` in the UMAT wrapper **before** handing
!! the model to `integrate_stress`. Inside the integrator those fields are
!! treated as frozen per-step constants; only dilation (and therefore eta_y)
!! evolves further via `mcsr_update_hardening`.
!!
!! ### snapshot / restore
!!
!! The integrator requires rollback capability. The snapshot captures the 12
!! entries that can change inside a substep:
!!
!! | Slot | Field         |
!! |------|---------------|
!! | 1    | G             |
!! | 2    | K             |
!! | 3    | eta_y         |
!! | 4    | dilation      |
!! | 5    | I_coeff       |
!! | 6    | M_lode        |
!! | 7    | switch_yield  (0.0 / 1.0) |
!! | 8–13 | eps_p(6)      |
!!
!! N_i and SUM_rate are NOT captured — rate-smoothing state survives rollback.

module mod_mcsr_model
   use mod_csm_kinds,      only: wp
   use mod_csm_model,      only: csm_model_t
   use mod_mcsr_types,     only: mcsr_params_t, mcsr_state_t
   use mod_mcsr_functions, only: mcsr_yield_fn, mcsr_flow_rule, mcsr_plastic_potential, &
                                  mcsr_update_hardening, mcsr_elastic_stiffness,         &
                                  mcsr_hardening_modulus, mcsr_update_rate_state
   implicit none
   private

   public :: mcsr_model_t, mcsr_from_props, mcsr_load_state, mcsr_save_state

   type, extends(csm_model_t) :: mcsr_model_t
      !! Mohr-Coulomb Strain Rate model.
      !!
      !! Extends `csm_model_t`. Parameters are fixed after construction via
      !! `mcsr_from_props`; rate-state is updated by `mcsr_update_rate_state`
      !! once per timestep; inner-loop state evolves via `update_hardening`.
      type(mcsr_params_t) :: params    !! Fixed model parameters
      type(mcsr_state_t)  :: state     !! Evolving internal state
   contains
      procedure :: yield_fn          => mcsr_yield_fn_method
      procedure :: flow_rule         => mcsr_flow_rule_method
      procedure :: plastic_potential => mcsr_plastic_potential_method
      procedure :: update_hardening  => mcsr_update_hardening_method
      procedure :: elastic_stiffness => mcsr_elastic_stiffness_method
      procedure :: snapshot          => mcsr_snapshot
      procedure :: restore           => mcsr_restore
      procedure :: hardening_modulus => mcsr_hardening_modulus_method
   end type mcsr_model_t

contains

   ! ---------------------------------------------------------------------------
   ! UMAT boundary helpers — only place PROPS/STATEV index magic lives
   ! ---------------------------------------------------------------------------

   function mcsr_from_props(props) result(model)
      !! Construct an `mcsr_model_t` from a flat PROPS array (UMAT convention).
      !!
      !! See `umat_mcsr_legacy.f90` module doc for the full PROPS layout.
      !! `switch_original` (PROPS(15)) is validated here: 0=Olzak-Perzyna is
      !! not supported in the new architecture.
      real(wp), intent(in) :: props(:)  !! PROPS array (18 entries expected)
      type(mcsr_model_t) :: model

      model%params%G_0        = props(1)
      model%params%nu         = props(2)
      model%params%M_tc       = props(3)
      model%params%N          = props(4)
      model%params%D_min      = props(5)
      model%params%h          = props(6)
      model%params%alpha_G    = props(7)
      model%params%alpha_K    = props(8)
      model%params%alpha_D    = props(9)
      model%params%D_part     = props(10)
      model%params%G_s        = props(11)
      model%params%ref_e_rate = props(12)

      model%params%switch_smooth = (nint(props(13)) /= 0)
      model%params%N_S           = nint(props(14))

      ! switch_original = 0 (Olzak-Perzyna) is not supported in new architecture
      if (size(props) >= 15) then
         if (nint(props(15)) == 0) then
            error stop "MCSR: switch_original=0 (Olzak-Perzyna) is not supported &
                       &in the new architecture. Use switch_original=1 (Wang)."
         end if
      end if

      model%params%yield_tol  = props(16)
      model%params%max_iters  = nint(props(17))
      ! props(18) = switch_plastic_integration — ignored; always uses euler_substep

      ! Apply legacy defaults from umat_mcsr_legacy
      if (model%params%ref_e_rate == 0.0_wp) model%params%ref_e_rate = 2.5e-5_wp
      if (model%params%alpha_D == 0.0_wp .and. model%params%alpha_G > 0.0_wp) &
         model%params%alpha_D = model%params%alpha_G
      if (model%params%alpha_K == 0.0_wp .and. model%params%alpha_G > 0.0_wp) &
         model%params%alpha_K = 2.5_wp * model%params%alpha_G
   end function mcsr_from_props

   subroutine mcsr_load_state(model, statev)
      !! Unpack a STATEV array into the model's internal state.
      !!
      !! STATEV layout (14 entries):
      !!   1  = G, 2 = K, 3 = eta_y, 4 = dilation, 5 = I_coeff,
      !!   6  = switch_yield (0.0/1.0), 7–12 = eps_p(6),
      !!   13 = N_i, 14 = SUM_rate
      type(mcsr_model_t), intent(inout) :: model
      real(wp),           intent(in)    :: statev(:)

      model%state%G            = statev(1)
      model%state%K            = statev(2)
      model%state%eta_y        = statev(3)
      model%state%dilation     = statev(4)
      model%state%I_coeff      = statev(5)
      model%state%switch_yield = (statev(6) /= 0.0_wp)
      model%state%eps_p        = statev(7:12)
      model%state%N_i          = nint(statev(13))
      model%state%SUM_rate     = statev(14)
      ! M_lode is NOT stored in STATEV — it is computed each step by mcsr_update_rate_state
      model%state%M_lode       = 0.0_wp
   end subroutine mcsr_load_state

   subroutine mcsr_save_state(model, statev)
      !! Pack the model's internal state back into a STATEV array.
      type(mcsr_model_t), intent(in)    :: model
      real(wp),           intent(inout) :: statev(:)

      statev(1)    = model%state%G
      statev(2)    = model%state%K
      statev(3)    = model%state%eta_y
      statev(4)    = model%state%dilation
      statev(5)    = model%state%I_coeff
      statev(6)    = merge(1.0_wp, 0.0_wp, model%state%switch_yield)
      statev(7:12) = model%state%eps_p
      statev(13)   = real(model%state%N_i, wp)
      statev(14)   = model%state%SUM_rate
   end subroutine mcsr_save_state

   ! ---------------------------------------------------------------------------
   ! Thin wrappers around pure functions — implement deferred procedures
   ! ---------------------------------------------------------------------------

   function mcsr_yield_fn_method(self, sig) result(F)
      class(mcsr_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: F
      F = mcsr_yield_fn(self%params, self%state, sig)
   end function mcsr_yield_fn_method

   function mcsr_flow_rule_method(self, sig) result(dF_by_dsig)
      class(mcsr_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: dF_by_dsig(6)
      dF_by_dsig = mcsr_flow_rule(self%params, self%state, sig)
   end function mcsr_flow_rule_method

   function mcsr_plastic_potential_method(self, sig) result(dG_by_dsig)
      class(mcsr_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: dG_by_dsig(6)
      dG_by_dsig = mcsr_plastic_potential(self%params, self%state, sig)
   end function mcsr_plastic_potential_method

   subroutine mcsr_update_hardening_method(self, deps_p)
      class(mcsr_model_t), intent(inout) :: self
      real(wp), intent(in) :: deps_p(6)
      call mcsr_update_hardening(self%params, self%state, deps_p)
   end subroutine mcsr_update_hardening_method

   function mcsr_elastic_stiffness_method(self) result(stiff_e)
      class(mcsr_model_t), intent(in) :: self
      real(wp) :: stiff_e(6,6)
      stiff_e = mcsr_elastic_stiffness(self%params, self%state)
   end function mcsr_elastic_stiffness_method

   function mcsr_hardening_modulus_method(self, sig, dg_by_dsig) result(H)
      class(mcsr_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6), dg_by_dsig(6)
      real(wp) :: H
      call mcsr_hardening_modulus(self%params, self%state, sig, dg_by_dsig, H)
   end function mcsr_hardening_modulus_method

   ! ---------------------------------------------------------------------------
   ! snapshot / restore — for integrator rollback
   ! ---------------------------------------------------------------------------

   subroutine mcsr_snapshot(self, saved)
      class(mcsr_model_t), intent(in)    :: self
      real(wp), allocatable, intent(out) :: saved(:)

      allocate(saved(13))
      saved(1)    = self%state%G
      saved(2)    = self%state%K
      saved(3)    = self%state%eta_y
      saved(4)    = self%state%dilation
      saved(5)    = self%state%I_coeff
      saved(6)    = self%state%M_lode
      saved(7)    = merge(1.0_wp, 0.0_wp, self%state%switch_yield)
      saved(8:13) = self%state%eps_p
   end subroutine mcsr_snapshot

   subroutine mcsr_restore(self, saved)
      class(mcsr_model_t), intent(inout) :: self
      real(wp),            intent(in)    :: saved(:)

      self%state%G            = saved(1)
      self%state%K            = saved(2)
      self%state%eta_y        = saved(3)
      self%state%dilation     = saved(4)
      self%state%I_coeff      = saved(5)
      self%state%M_lode       = saved(6)
      self%state%switch_yield = (saved(7) /= 0.0_wp)
      self%state%eps_p        = saved(8:13)
   end subroutine mcsr_restore

end module mod_mcsr_model
