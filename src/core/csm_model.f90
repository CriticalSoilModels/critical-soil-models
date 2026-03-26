! Abstract base type for all constitutive soil models.
! Every concrete model (MCSS, NorSand, MC-SR, etc.) extends this type
! and implements the seven deferred procedures.

module mod_csm_model
   use mod_csm_kinds, only: wp
   implicit none

   ! ---------------------------------------------------------------------------
   ! Abstract base type
   ! ---------------------------------------------------------------------------
   ! Stress is NOT stored here — it is always passed explicitly to procedures.
   ! Params (fixed) and state (evolving) live in the concrete extended type.
   !
   ! snapshot/restore exist so the integrator can roll back state during the
   ! two-estimate substep loop without any knowledge of what state variables
   ! a given model has. Index magic lives only inside each model's snapshot
   ! and restore — never in the integrator.

   type, abstract :: csm_model_t
   contains
      procedure(yield_fn_iface),    deferred :: yield_fn
      procedure(flow_rule_iface),   deferred :: flow_rule
      procedure(plastic_pot_iface), deferred :: plastic_potential
      procedure(harden_iface),      deferred :: update_hardening
      procedure(stiffness_iface),   deferred :: elastic_stiffness
      procedure(snapshot_iface),    deferred :: snapshot
      procedure(restore_iface),     deferred :: restore
   end type csm_model_t

   ! ---------------------------------------------------------------------------
   ! Abstract interfaces — defines the contract every model must satisfy
   ! ---------------------------------------------------------------------------

   abstract interface

      ! F = model%yield_fn(sig)
      ! Returns scalar yield function value. F < 0 = elastic, F = 0 = on surface.
      function yield_fn_iface(self, sig) result(F)
         import :: csm_model_t, wp
         class(csm_model_t), intent(in) :: self
         real(wp), intent(in) :: sig(6)   ! internal Voigt convention
         real(wp) :: F
      end function yield_fn_iface

      ! dF_by_dsig = model%flow_rule(sig)
      ! Gradient of the yield surface — used for the consistency condition.
      function flow_rule_iface(self, sig) result(dF_by_dsig)
         import :: csm_model_t, wp
         class(csm_model_t), intent(in) :: self
         real(wp), intent(in) :: sig(6)
         real(wp) :: dF_by_dsig(6)
      end function flow_rule_iface

      ! dg_plas_by_dsig = model%plastic_potential(sig)
      ! Flow direction. For associated flow: identical to flow_rule.
      ! For non-associated flow (e.g. MC): uses dilation angle instead of friction angle.
      function plastic_pot_iface(self, sig) result(dg_plas_by_dsig)
         import :: csm_model_t, wp
         class(csm_model_t), intent(in) :: self
         real(wp), intent(in) :: sig(6)
         real(wp) :: dg_plas_by_dsig(6)
      end function plastic_pot_iface

      ! call model%update_hardening(deps_p)
      ! Updates all internal state variables given a plastic strain increment.
      ! This is the only place internal state changes during stress integration.
      subroutine harden_iface(self, deps_p)
         import :: csm_model_t, wp
         class(csm_model_t), intent(inout) :: self
         real(wp), intent(in) :: deps_p(6)
      end subroutine harden_iface

      ! stiff_e = model%elastic_stiffness()
      ! Returns the 6x6 elastic stiffness matrix.
      ! For models with pressure-dependent stiffness this uses current state.
      function stiffness_iface(self) result(stiff_e)
         import :: csm_model_t, wp
         class(csm_model_t), intent(in) :: self
         real(wp) :: stiff_e(6,6)
      end function stiffness_iface

      ! call model%snapshot(saved)
      ! Serialises all state variables into a flat real array.
      ! The integrator holds this array opaquely — it does not interpret the contents.
      ! Index layout is defined only inside each model's implementation.
      subroutine snapshot_iface(self, saved)
         import :: csm_model_t, wp
         class(csm_model_t), intent(in)               :: self
         real(wp), allocatable, intent(out)            :: saved(:)
      end subroutine snapshot_iface

      ! call model%restore(saved)
      ! Deserialises state from a flat array previously produced by snapshot.
      ! Inverse of snapshot — index layout must match.
      subroutine restore_iface(self, saved)
         import :: csm_model_t, wp
         class(csm_model_t), intent(inout) :: self
         real(wp),            intent(in)   :: saved(:)
      end subroutine restore_iface

   end interface

end module mod_csm_model
