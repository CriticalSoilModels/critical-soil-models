! Dispatch wrapper for stress integration methods.
!
! Provides a single call site for all integrators so that model UMATs
! and test drivers do not call integrators directly.  The method is
! selected at runtime via a string, enabling comparison across methods
! without recompilation.
!
! Named constants (use these — do not write bare string literals):
!   INTEGRATION_EULER       = "euler"
!   INTEGRATION_ORTIZ_SIMO  = "ortiz_simo"

module mod_integrate_stress
   use mod_csm_kinds,     only: wp
   use mod_csm_model,     only: csm_model_t
   use mod_euler_substep, only: euler_substep
   implicit none
   private
   public :: integrate_stress
   public :: INTEGRATION_EULER, INTEGRATION_ORTIZ_SIMO

   character(len=*), parameter :: INTEGRATION_EULER      = "euler"
   character(len=*), parameter :: INTEGRATION_ORTIZ_SIMO = "ortiz_simo"

contains

   subroutine integrate_stress(model, sig, deps, ftol, stol, method)
      !! Integrate the stress state by one increment using the named method.
      !!
      !! model  — concrete csm_model_t (provides yield_fn, flow_rule, etc.)
      !! sig    — stress vector [6], updated in place
      !! deps   — strain increment [6]
      !! ftol   — yield function tolerance
      !! stol   — substep error tolerance (used by substepping methods)
      !! method — integration method name; use the INTEGRATION_* constants
      class(csm_model_t), intent(inout) :: model
      real(wp),           intent(inout) :: sig(6)
      real(wp),           intent(in)    :: deps(6)
      real(wp),           intent(in)    :: ftol, stol
      character(len=*),   intent(in)    :: method

      select case(trim(method))

      case(INTEGRATION_EULER)
         call euler_substep(model, sig, deps, ftol, stol)

      case(INTEGRATION_ORTIZ_SIMO)
         error stop "integrate_stress: ortiz_simo not yet implemented"

      case default
         error stop "integrate_stress: unknown integration method"

      end select
   end subroutine integrate_stress

end module mod_integrate_stress
