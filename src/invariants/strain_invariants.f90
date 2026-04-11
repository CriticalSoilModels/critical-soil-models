! Strain invariant functions for use in constitutive models.
!
! Computation chain:
!   eps(6)
!     -> eps_v   = calc_eps_vol_inv(eps)          ! volumetric strain = tr(ε)
!     -> dev     = calc_dev_strain(eps, eps_v)    ! deviatoric strain tensor
!     -> eps_q   = calc_eps_q_inv(dev)            ! deviatoric strain invariant: sqrt(2/3 * dev:dev)
!
! All-in-one: calc_eps_inv(eps, eps_v, eps_q)

module mod_strain_invariants
   use mod_csm_kinds, only: wp
   use mod_voigt_utils, only: calc_two_norm_tensor_strain
   implicit none

contains

   pure function calc_eps_vol_inv(strain) result(eps_vol)
      !! Volumetric strain invariant: εv = tr(ε) = ε11 + ε22 + ε33
      real(wp), intent(in) :: strain(6)
      real(wp) :: eps_vol

      eps_vol = strain(1) + strain(2) + strain(3)
   end function calc_eps_vol_inv

   pure function calc_dev_strain(strain, eps_v) result(dev_strain)
      !! Deviatoric strain tensor: e = ε - (εv/3)*I
      real(wp), intent(in) :: strain(6), eps_v
      real(wp) :: dev_strain(6)

      dev_strain    = strain
      dev_strain(1) = strain(1) - eps_v/3.0_wp
      dev_strain(2) = strain(2) - eps_v/3.0_wp
      dev_strain(3) = strain(3) - eps_v/3.0_wp
   end function calc_dev_strain

   pure function calc_eps_q_inv(dev_strain) result(eps_q)
      !! Deviatoric strain invariant: εq = sqrt(2/3) * ||e||
      !! TODO: verify invariant definition against reference
      real(wp), intent(in) :: dev_strain(6)
      real(wp) :: eps_q

      eps_q = calc_two_norm_tensor_strain(dev_strain) * sqrt(2.0_wp / 3.0_wp)
   end function calc_eps_q_inv

   pure subroutine calc_eps_inv(eps, eps_v, eps_q)
      !! All-in-one: returns volumetric and deviatoric strain invariants.
      real(wp), intent(in)  :: eps(6)
      real(wp), intent(out) :: eps_v, eps_q

      real(wp) :: dev(6)

      eps_v = calc_eps_vol_inv(eps)
      dev   = calc_dev_strain(eps, eps_v)
      eps_q = calc_eps_q_inv(dev)
   end subroutine calc_eps_inv

end module mod_strain_invariants
