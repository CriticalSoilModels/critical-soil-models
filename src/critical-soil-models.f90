module critical_soil_models
  use mod_mcss_esm, only: UMAT_MCSS => UMAT_MohrCoulombStrainSoftening
  use mod_lin_elastic, only: UMAT_Lin => UMAT

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