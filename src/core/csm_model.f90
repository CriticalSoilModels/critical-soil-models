!! Abstract base type for all constitutive soil models.
!!
!! Every concrete model (MCSS, NorSand, MC-SR, etc.) extends `csm_model_t`
!! and implements the eight deferred procedures. The integrators in
!! `mod_integrate_stress` depend only on this interface — they have no
!! knowledge of any specific model.
!!
!! ### Design notes
!! - Stress is **not** stored here; it is always passed explicitly.
!! - Fixed parameters and evolving state live in the concrete extended type.
!! - `snapshot`/`restore` let the integrator roll back state during the
!!   two-estimate substep loop without knowing what variables any given model has.

module mod_csm_model
   use mod_csm_kinds, only: wp
   implicit none

   ! ---------------------------------------------------------------------------
   ! Abstract base type
   ! ---------------------------------------------------------------------------

   type, abstract :: csm_model_t
      !! Abstract base for all CSM constitutive models.
      !!
      !! Concrete models extend this type and implement every deferred procedure.
      !! The integrator interacts exclusively through this interface.
      !!
      !! ### Stiffness procedures
      !!
      !! `elastic_stiffness` — deferred; returns D_e, used by the integrator for
      !! the elastic predictor and plastic correction. Every model must implement it.
      !!
      !! `consistent_tangent` — non-deferred; returns the stiffness tensor passed
      !! back to the FEA solver as DDSDDE. Defaults to `elastic_stiffness`.
      !! Override in models that implement the consistent elastoplastic tangent for
      !! faster FEA convergence.
      !!
      !! `pre_step` — non-deferred; called by the integrator at the start of every
      !! substep, before `elastic_stiffness` is evaluated. Default is a no-op.
      !! Override to update any state that depends on the current stress — for example,
      !! pressure-dependent elastic moduli or rate state. Receiving `sig` at substep
      !! resolution matters for large increments (e.g. implicit MPM time stepping).
   contains
      procedure(yield_fn_iface),    deferred :: yield_fn
      procedure(flow_rule_iface),   deferred :: flow_rule
      procedure(plastic_pot_iface), deferred :: plastic_potential
      procedure(harden_iface),      deferred :: update_hardening
      procedure(stiffness_iface),   deferred :: elastic_stiffness
      procedure(snapshot_iface),    deferred :: snapshot
      procedure(restore_iface),     deferred :: restore
      procedure(harden_mod_iface),  deferred :: hardening_modulus
      procedure :: consistent_tangent => default_consistent_tangent
      procedure :: pre_step          => default_pre_step
   end type csm_model_t

   ! ---------------------------------------------------------------------------
   ! Abstract interfaces — defines the contract every model must satisfy
   ! ---------------------------------------------------------------------------

   abstract interface

      function yield_fn_iface(self, sig) result(F)
         !! Evaluate the yield function F(σ).
         !!
         !! Returns a scalar: F < 0 elastic, F = 0 on the yield surface.
         import :: csm_model_t, wp
         class(csm_model_t), intent(in) :: self
         real(wp), intent(in) :: sig(6)   !! Stress vector (Voigt, internal convention)
         real(wp) :: F                    !! Yield function value
      end function yield_fn_iface

      function flow_rule_iface(self, sig) result(dF_by_dsig)
         !! Gradient of the yield surface, ∂F/∂σ.
         !!
         !! Used in the consistency condition to determine the plastic multiplier.
         import :: csm_model_t, wp
         class(csm_model_t), intent(in) :: self
         real(wp), intent(in) :: sig(6)        !! Stress vector (Voigt)
         real(wp) :: dF_by_dsig(6)             !! Yield surface gradient
      end function flow_rule_iface

      function plastic_pot_iface(self, sig) result(dg_plas_by_dsig)
         !! Plastic flow direction, ∂G/∂σ.
         !!
         !! For associated flow rules this is identical to `flow_rule`.
         !! For non-associated flow (e.g. Mohr-Coulomb) the dilation angle
         !! replaces the friction angle.
         import :: csm_model_t, wp
         class(csm_model_t), intent(in) :: self
         real(wp), intent(in) :: sig(6)        !! Stress vector (Voigt)
         real(wp) :: dg_plas_by_dsig(6)        !! Plastic flow direction
      end function plastic_pot_iface

      subroutine harden_iface(self, deps_p)
         !! Update all internal state variables given a plastic strain increment.
         !!
         !! This is the **only** place internal state changes during stress
         !! integration. The integrator calls this once per accepted substep.
         import :: csm_model_t, wp
         class(csm_model_t), intent(inout) :: self
         real(wp), intent(in) :: deps_p(6)     !! Plastic strain increment (Voigt)
      end subroutine harden_iface

      function stiffness_iface(self) result(stiff_e)
         !! Elastic stiffness matrix D_e (6×6).
         !!
         !! For models with pressure-dependent stiffness this uses the
         !! current internal state.
         import :: csm_model_t, wp
         class(csm_model_t), intent(in) :: self
         real(wp) :: stiff_e(6,6)              !! Elastic stiffness [stress units]
      end function stiffness_iface

      subroutine snapshot_iface(self, saved)
         !! Serialise all state variables into a flat real array.
         !!
         !! The integrator holds `saved` opaquely and passes it back to
         !! `restore` if a substep is rejected. Index layout is defined
         !! only inside each model's implementation.
         import :: csm_model_t, wp
         class(csm_model_t), intent(in)    :: self
         real(wp), allocatable, intent(out) :: saved(:)  !! Serialised state
      end subroutine snapshot_iface

      subroutine restore_iface(self, saved)
         !! Deserialise state from an array previously produced by `snapshot`.
         !!
         !! Inverse of `snapshot`; index layout must match exactly.
         import :: csm_model_t, wp
         class(csm_model_t), intent(inout) :: self
         real(wp),            intent(in)   :: saved(:)   !! Serialised state
      end subroutine restore_iface

      function harden_mod_iface(self, sig, dg_by_dsig) result(H)
         !! Hardening modulus H = (∂F/∂κ)(dκ/dλ).
         !!
         !! Used in the consistency condition to solve for the plastic multiplier:
         !!
         !! - H > 0 → softening (yield surface shrinks)
         !! - H < 0 → hardening (yield surface expands)
         !! - H = 0 → perfectly plastic
         !!
         !! `dg_by_dsig` is needed to evaluate dε_p/dλ when κ depends on
         !! equivalent plastic strain.
         import :: csm_model_t, wp
         class(csm_model_t), intent(in) :: self
         real(wp),           intent(in) :: sig(6)        !! Stress vector (Voigt)
         real(wp),           intent(in) :: dg_by_dsig(6) !! Plastic flow direction
         real(wp) :: H                                   !! Hardening modulus
      end function harden_mod_iface

   end interface

contains

   subroutine default_pre_step(self, sig)
      !! Default pre-step hook: no-op.
      !!
      !! Override in models whose state depends on the current stress — for
      !! example, to recompute pressure-dependent elastic moduli before each
      !! substep. Called by the integrator before `elastic_stiffness`.
      class(csm_model_t), intent(inout) :: self
      real(wp),           intent(in)    :: sig(6)
   end subroutine default_pre_step

   function default_consistent_tangent(self) result(stiff)
      !! Default consistent tangent: returns the elastic stiffness D_e.
      !!
      !! UMAT wrappers pass this to DDSDDE. For most models the elastic tangent
      !! is a safe, conservative choice — FEA convergence is slower but correct.
      !! Override this procedure in models that implement the full consistent
      !! elastoplastic tangent.
      class(csm_model_t), intent(in) :: self
      real(wp) :: stiff(6,6)
      stiff = self%elastic_stiffness()
   end function default_consistent_tangent

end module mod_csm_model
