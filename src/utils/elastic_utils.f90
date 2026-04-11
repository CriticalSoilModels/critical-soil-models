! Utilities for isotropic linear elasticity.
! All functions are pure one-liners; no side effects.

module mod_elastic_utils
   use mod_csm_kinds, only: wp
   implicit none
   private

   public :: calc_stiffness_GK
   public :: calc_K_from_G_nu
   public :: calc_G_from_E_nu
   public :: calc_nu_from_G_K

contains

   pure function calc_stiffness_GK(G, K) result(stiff_e)
      !! Isotropic linear elastic stiffness matrix in Voigt notation [11,22,33,12,13,23]
      real(wp), intent(in) :: G  !! Shear modulus [kPa]
      real(wp), intent(in) :: K  !! Bulk modulus [kPa]
      real(wp) :: stiff_e(6,6)

      real(wp) :: lame_1, lame_2

      lame_1 = K + 4.0_wp*G/3.0_wp   ! normal diagonal
      lame_2 = K - 2.0_wp*G/3.0_wp   ! normal off-diagonal (Lamé lambda)

      stiff_e          = 0.0_wp
      stiff_e(1:3,1:3) = lame_2
      stiff_e(1,1) = lame_1;  stiff_e(2,2) = lame_1;  stiff_e(3,3) = lame_1
      stiff_e(4,4) = G;       stiff_e(5,5) = G;        stiff_e(6,6) = G

   end function calc_stiffness_GK

   pure function calc_K_from_G_nu(G, nu) result(K)
      !! K = 2G(1 + nu) / (3(1 - 2nu))
      real(wp), intent(in) :: G, nu
      real(wp) :: K
      K = 2.0_wp*G*(1.0_wp + nu) / (3.0_wp*(1.0_wp - 2.0_wp*nu))
   end function calc_K_from_G_nu

   pure function calc_G_from_E_nu(E, nu) result(G)
      !! G = E / (2(1 + nu))
      real(wp), intent(in) :: E, nu
      real(wp) :: G
      G = E / (2.0_wp*(1.0_wp + nu))
   end function calc_G_from_E_nu

   pure function calc_nu_from_G_K(G, K) result(nu)
      !! nu = (3K - 2G) / (2(3K + G))
      real(wp), intent(in) :: G, K
      real(wp) :: nu
      nu = (3.0_wp*K - 2.0_wp*G) / (2.0_wp*(3.0_wp*K + G))
   end function calc_nu_from_G_K

end module mod_elastic_utils
