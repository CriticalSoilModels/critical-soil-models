module mod_stress_invariants
   use kind_precision_module, only : dp, i32

   use mod_voigt_functions, only : TwoNormTensor, TwoNormTensor_strain, calc_dev_stess

   implicit none

contains

   pure function calc_mean_stress(stress) result(mean_stress)
      real(kind = dp), intent(in) :: stress(6)
      real(kind = dp) :: mean_stress

      ! Local variables
      integer(kind = i32) :: i

      ! Calc the mean stress
      mean_stress = 0.0_dp

      do i = 1, 3
         mean_stress = mean_stress + stress(i)
      end do

      mean_stress = mean_stress / 3.0_dp
   end function calc_mean_stress

   pure function calc_theta_s(J2, J3) result(lode_angle)
      real(kind=dp), intent(in) :: J2, J3
      real(kind=dp) :: lode_angle

      ! Local variables
      real(kind=dp) :: sin3theta

      ! Reference: D. M. Potts and L. Zdravković, Finite element analysis in geotechnical engineering. 2: Application. London: Telford, 2001.
      ! Formula for Lode angle: Potts and Zdravković page 186
      ! Information on bounding and intuition: Potts and Zdravković page 116
      ! Wikipedia has information but I think they messed up the bounds
      
      ! Trx compression: s1 >= s2 = s3    => -PI/6
      ! Shear          : s2 = (s1 + s3)/2 =>  0
      ! Trx Extension  : s1 = s2 >= s3    =>  PI/6
      
      if (J2 > 0.0_dp) then
         ! Ensure correct scaling for sin(3*theta)
         sin3theta = 0.5_dp * J3 * (3.0_dp / J2)**1.5_dp
      else
         ! Assume triaxial compression if J2 is zero or negative
         sin3theta = -1.0_dp
      endif

      ! Clamp sin3theta between -1 and 1 for numerical stability
      if (sin3theta < -1.0_dp) sin3theta = -1.0_dp
      if (sin3theta >  1.0_dp) sin3theta =  1.0_dp

      ! Lode angle calculation
      lode_angle = -asin(sin3theta) / 3.0_dp

   end function calc_theta_s

   pure function calc_inc_driver_J3_invariant(dev) result(J3)
      real(kind = dp), intent(in) :: dev(6)
      real(kind = dp) :: J3

      ! Local variables
      real(kind = dp) :: first_term, second_term, third_term, fourth_term, fifth_term, mean_stress
      real(kind = dp), parameter :: ONE   = 1.0_dp, &
         TWO   = 2.0_dp, &
         THREE = 3.0_dp
      ! Incremental driver Voigt order
      !={
      !    11 (xx),
      !    22 (yy),
      !    33 (zz),
      !    12 (xy),
      !    13 (xz),
      !    23 (yz)
      !}

      ! Calc the J3 Stress invariant
      first_term  = product( dev(1:3) )

      ! Calc -1.0 * (sigma_x - p) * tau_{yz}^{2}
      second_term = - dev(1) * dev(6)**2

      ! Calc -1.0 * (sigma_y - p) * tau_{zx}^{2}
      third_term  = - dev(2) * dev(5)**2

      ! Calc -1.0 * (sigma_z - p) * tau_{xy}^{2}
      fourth_term = - dev(3) * dev(4)**2

      ! Calc 2.0 * tau_{xy} * tau_{yz} * tau_{zx}
      fifth_term  = TWO * product( dev(4:6) )

      J3 = first_term + second_term + third_term + fourth_term + fifth_term

   end function calc_inc_driver_J3_invariant

   pure function calc_J2_invariant(dev) result(J2)
      real(kind = dp), intent(in) :: dev(6)
      real(kind = dp) :: J2

      ! Calc the sqrt of the J2 invariant
      ! TwoNormTensor calculates the sqrt{inner product of a symmetric stress tensor}
      call TwoNormTensor(dev, 6, J2)

      J2 = 0.5_dp * J2**2
   end function calc_J2_invariant

   pure function calc_q_invariant(J2) result(q)
      real(kind = dp), intent(in) :: J2
      real(kind = dp) :: q

      q = sqrt(3.0_dp * J2)

   end function calc_q_invariant

   pure function calc_Zam_J3_invariant(stress) result(J3)
      real(kind = dp), intent(in) :: stress(6)
      real(kind = dp) :: J3

      ! Local variables
      real(kind = dp) :: dev(6)

      ! store a copy of the stress
      dev = stress

      ! Remove the pressure
      dev(1:3) = dev(1:3) - sum(stress(1:3)) / 3.0_dp

      ! Calc the J3 Stress invariant
      J3 = dev(1)*dev(2)*dev(3) - dev(1)*dev(6)**2 - dev(2)*dev(4)**2 - dev(3)*dev(5)**2 + 2.0*dev(4)*dev(5)*dev(6)

   end function calc_Zam_J3_invariant

   pure subroutine Get_invariants(Sig, p, q, theta)
      !*********************************************************************
      ! Takes the stress tensor Sig and return invariants p, q, and theta  *
      !																	 *
      !*********************************************************************
      implicit none
      !input variables
      double precision, dimension(6), intent(in):: Sig
      !output variables
      double precision, intent(out)::p, q, theta
      !local variables
      double precision:: dev(6), J2, J3

      p = calc_mean_stress(Sig)

      dev = calc_dev_stess(Sig, p)

      J2 = calc_J2_invariant(dev)

      q = calc_q_invariant(J2) ! deviatoric stress

      !J3 stress invariant
      J3 = calc_inc_driver_J3_invariant(dev)

      theta = calc_theta_s(J2, J3) !Lode's angle

   end subroutine Get_invariants
   
end module mod_stress_invariants
