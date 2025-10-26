module critical_soil_models
   use mod_mcss_esm, only: umat_mc_strain_softening => UMAT_MohrCoulombStrainSoftening
   use mod_lin_elastic, only: umat_linear => UMAT
   use mod_bingham, only: umat_bingham => umat
   use mod_mc_strain_rate, only: umat_mc_strain_rate
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
