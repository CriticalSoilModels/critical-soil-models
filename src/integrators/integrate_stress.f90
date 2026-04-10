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
   use mod_euler_substep, only: euler_substep, integrator_params_t, DEFAULT_INTEGRATOR_PARAMS
   use mod_cprm,          only: cprm_step
   implicit none
   private
   public :: integrate_stress
   public :: integrator_params_t, DEFAULT_INTEGRATOR_PARAMS
   public :: INTEGRATION_EULER, INTEGRATION_ORTIZ_SIMO

   character(len=*), parameter :: INTEGRATION_EULER      = "euler"
   character(len=*), parameter :: INTEGRATION_ORTIZ_SIMO = "ortiz_simo"

contains

   subroutine integrate_stress(model, sig, deps, method, iparams)
      !! Integrate the stress state by one increment using the named method.
      !!
      !! model   — concrete csm_model_t (provides yield_fn, flow_rule, etc.)
      !! sig     — stress vector [6], updated in place
      !! deps    — strain increment [6]
      !! method  — integration method name; use the INTEGRATION_* constants
      !! iparams — tolerances and step-size limits (optional; defaults used if absent)
      class(csm_model_t),        intent(inout)        :: model
      real(wp),                  intent(inout)        :: sig(6)
      real(wp),                  intent(in)           :: deps(6)
      character(len=*),          intent(in)           :: method
      type(integrator_params_t), intent(in), optional :: iparams

      select case(trim(method))

      case(INTEGRATION_EULER)
         call euler_substep(model, sig, deps, iparams)

      case(INTEGRATION_ORTIZ_SIMO)
         call cprm_step(model, sig, deps, iparams)

      case default
         error stop "integrate_stress: unknown integration method"

      end select
   end subroutine integrate_stress

end module mod_integrate_stress
