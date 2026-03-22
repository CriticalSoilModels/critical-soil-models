! Module for the the functions for the strain invariants
module mod_strain_invariants
   use stdlib_kinds, only: dp
   use mod_voigt_utils  , only: calc_two_norm_tensor_strain
   implicit none

contains
   subroutine calc_eps_invariants(eps, eps_v, eps_q)
      !! Takes the strain tensor and returns volumetric and deviatoric strain invariants
      implicit none
      real(dp), dimension(6), intent(in)  :: eps
      real(dp),               intent(out) :: eps_v, eps_q

      ! Local variables
      real(dp) :: dev(6)

      eps_v = calc_eps_vol(eps)
      dev   = calc_dev_strain(eps, eps_v)
      eps_q = calc_eps_q(dev)
   end subroutine calc_eps_invariants

   pure function calc_eps_vol(strain) result(eps_vol)
      ! Calc the volumetric strain invariant
      real(dp), intent(in) :: strain(6)
      real(dp) :: eps_vol

      ! Init to zero to make sure nothing weird happens
      eps_vol = 0.0
      ! Calc the volumetric strain invariant Tr(\epsilon)
      eps_vol = strain(1) + strain(2) + strain(3)

   end function calc_eps_vol

   pure function calc_dev_strain(strain, eps_v) result(dev_strain)
      ! Calc the deviatoric strain voigt vector
      real(dp), intent(in) :: strain(6), eps_v
      real(dp) :: dev_strain(6)

      ! Local variables
      integer :: i

      ! Store a copy of the strain voigt vector
      dev_strain = strain

      ! Subtract off the mean volumetric strain
      do i = 1, 3
         dev_strain(i) = dev_strain(i) - eps_v/3.0_dp
      end do
   end function calc_dev_strain

   pure function calc_eps_q(dev_strain) result(eps_q)
      ! Calc the derivatoric strain invariant
      ! TODO: Need to check that the invariants are correct
      real(dp), intent(in) :: dev_strain(6)
      real(dp) :: eps_q

      eps_q = calc_two_norm_tensor_strain(dev_strain) * sqrt(2.0_dp / 3.0_dp)
   end function calc_eps_q
   
end module mod_strain_invariants
