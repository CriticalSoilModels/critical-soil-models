!! Top-level library module — single import point for programmatic consumers.
!!
!! Re-exports the complete public API of the library.  A downstream project
!! added as an fpm dependency only needs:
!!
!!   use critical_soil_models, only: wp, csm_model_t, integrate_stress, ...
!!
!! ### What is NOT re-exported here
!! - Internal pure constitutive functions (`mod_*_functions`): these are
!!   lower-level building blocks.  If you need GPU-callable pure functions
!!   import directly from the model-specific `mod_*_functions` module.
!! - Legacy UMAT wrappers: use the Abaqus UMAT entry point in `src/umat.f90`.
!! - Invariant and convention helpers (`mod_stress_invariants`, etc.):
!!   import from `src/invariants/` directly if needed.

module critical_soil_models

   ! ---------------------------------------------------------------------------
   ! Precision
   ! ---------------------------------------------------------------------------
   use mod_csm_kinds, only: wp

   ! ---------------------------------------------------------------------------
   ! Abstract base type — needed for polymorphic use of any model
   ! ---------------------------------------------------------------------------
   use mod_csm_model, only: csm_model_t

   ! ---------------------------------------------------------------------------
   ! Stress integration dispatch
   ! ---------------------------------------------------------------------------
   use mod_integrate_stress, only: integrate_stress,          &
                                   integrator_params_t,       &
                                   DEFAULT_INTEGRATOR_PARAMS, &
                                   INTEGRATION_EULER,         &
                                   INTEGRATION_CPA

   ! ---------------------------------------------------------------------------
   ! Linear Elastic
   ! ---------------------------------------------------------------------------
   use mod_le_types,            only: le_params_t, le_state_t
   use mod_linear_elastic_model, only: linear_elastic_model_t, le_from_props

   ! ---------------------------------------------------------------------------
   ! Mohr-Coulomb Strain Softening (MCSS)
   ! ---------------------------------------------------------------------------
   use mod_mcss_types,  only: mcss_params_t, mcss_state_t,     &
                               abbo_sloan_params_t, DEFAULT_AS_PARAMS
   use mod_mcss_model,  only: mcss_model_t,                    &
                               mcss_from_props,                 &
                               mcss_load_state,                 &
                               mcss_save_state

   ! ---------------------------------------------------------------------------
   ! NorSand
   ! ---------------------------------------------------------------------------
   use mod_norsand_types,  only: norsand_params_t, norsand_state_t
   use mod_norsand_model,  only: norsand_model_t,               &
                                 norsand_from_props,             &
                                 norsand_load_state,             &
                                 norsand_save_state

   ! ---------------------------------------------------------------------------
   ! Mohr-Coulomb Strain Rate (MCSR)
   ! ---------------------------------------------------------------------------
   use mod_mcsr_types,     only: mcsr_params_t, mcsr_state_t
   use mod_mcsr_model,     only: mcsr_model_t,                  &
                                 mcsr_from_props,                &
                                 mcsr_load_state,                &
                                 mcsr_save_state
   use mod_mcsr_functions, only: mcsr_update_rate_state

   implicit none

contains

   subroutine print_banner()
      implicit none
      print *, " "
      print *, "---------------------------------------"
      print *, " Critical Soil Models Module for Fortran"
      print *, "             WaveHello"
      print *, "---------------------------------------"
      print *, " "
   end subroutine print_banner

end module critical_soil_models
