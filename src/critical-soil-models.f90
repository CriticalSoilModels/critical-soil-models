! Top-level library module — exports model types for programmatic use.
! UMATs are accessed through the single UMAT entry point in src/umat.f90.
module critical_soil_models
   use mod_linear_elastic_model, only: linear_elastic_model_t
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
