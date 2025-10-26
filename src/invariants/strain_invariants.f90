! Module for the the functions for the strain invariants
module mod_strain_invariants
   use kind_precision_module, only: real_type => dp
   use mod_voigt_functions  , only: TwoNormTensor_strain
   implicit none

contains
   Subroutine Get_strain_invariants(Eps, Eps_v, Eps_q)
      !*********************************************************************
      ! Takes the strain tensor and returns deviatoric and vol. strain     *
      !																	 *
      !*********************************************************************
      implicit none
      !input
      double precision, dimension(6), intent(in):: Eps
      !output
      double precision, intent(out):: Eps_v, Eps_q
      !local variables
      double precision:: dev(6)
      ! Eps_v=Eps(1)+Eps(2)+Eps(3)! vol strain
      Eps_v = calc_eps_vol_invariant(Eps)

      ! dev=Eps
      ! dev(1)=dev(1)-(Eps_v/3.0)
      ! dev(2)=dev(2)-(Eps_v/3.0)
      ! dev(3)=dev(3)-(Eps_v/3.0)!deviatoric strain tensor
      dev = calc_dev_strain(Eps, Eps_v)

      ! call TwoNormTensor_strain(dev, 6, Eps_q)
      ! Eps_q=Eps_q*sqrt(2.0/3.0) ! dev strain
      Eps_q = calc_eps_q_invariant(dev)
   end subroutine  Get_strain_invariants

   pure function calc_eps_vol_invariant(strain) result(eps_vol)
      ! Calc the volumetric strain invariant
      real(kind = real_type), intent(in) :: strain(6)
      real(kind = real_type) :: eps_vol

      ! Init to zero to make sure nothing weird happens
      eps_vol = 0.0
      ! Calc the volumetric strain invariant Tr(\epsilon)
      eps_vol = strain(1) + strain(2) + strain(3)

   end function calc_eps_vol_invariant

   pure function calc_dev_strain(strain, eps_v) result(dev_strain)
      ! Calc the deviatoric strain voigt vector
      real(kind = real_type), intent(in) :: strain(6), eps_v
      real(kind = real_type) :: dev_strain(6)

      ! Local variables
      integer :: i

      ! Store a copy of the strain voigt vector
      dev_strain = strain

      ! Subtract off the mean volumetric strain
      do i = 1, 3
         dev_strain(i) = dev_strain(i) - eps_v/3.0_real_type
      end do
   end function calc_dev_strain

   pure function calc_eps_q_invariant(dev_strain) result(eps_q)
      ! Calc the derivatoric strain invariant
      ! TODO: Need to check that the invariants are correct
      real(kind = real_type), intent(in) :: dev_strain(6)
      real(kind = real_type) :: eps_q

      call TwoNormTensor_strain(dev_strain, 6, eps_q)

      eps_q = eps_q * sqrt(2.0_real_type / 3.0_real_type)
   end function calc_eps_q_invariant
   
end module mod_strain_invariants
