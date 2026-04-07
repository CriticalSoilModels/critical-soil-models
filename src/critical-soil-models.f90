!! Top-level library module — exports model types for programmatic use.
!!
!! Re-exports the public API for consumers who embed the library directly
!! rather than linking against the Abaqus UMAT interface. UMATs are accessed
!! through the single `UMAT` entry point in `src/umat.f90`.
module critical_soil_models
   use mod_le_functions,         only: le_params_t, le_state_t, &
                                        le_yield_fn, le_flow_rule, le_plastic_potential, &
                                        le_update_hardening, le_elastic_stiffness
   use mod_linear_elastic_model, only: linear_elastic_model_t, le_from_props
   use mod_mcss_model,           only: mcss_model_t
   implicit none

contains
   subroutine print_banner()
      implicit none
      print *, " "
      print *, "---------------------------------------"
      print *, " Critical Soil Models Module for Fortran"
      print *, "             WaveHello"
      print *, "---------------------------------------"
      ! print *, " with contributions from: "
      ! add more contributors here :)
      print *, " "
   end subroutine print_banner

end module critical_soil_models
