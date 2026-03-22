module mod_stress_invariants
   use stdlib_kinds, only: dp, i32 => int32

   use mod_voigt_utils, only : calc_two_norm_tensor, calc_two_norm_tensor_strain, calc_dev_stress

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

   pure function calc_lode_angle(J2, J3) result(lode_angle)
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

   end function calc_lode_angle

   pure function calc_J3(dev) result(J3)
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

   end function calc_J3

   pure function calc_J2(dev) result(J2)
      real(kind = dp), intent(in) :: dev(6)
      real(kind = dp) :: J2

      ! Calc the sqrt of the J2 invariant
      J2 = 0.5_dp * calc_two_norm_tensor(dev)**2
   end function calc_J2

   pure function calc_q(J2) result(q)
      real(kind = dp), intent(in) :: J2
      real(kind = dp) :: q

      q = sqrt(3.0_dp * J2)

   end function calc_q

   pure subroutine calc_sig_invariants(sig, p, q, lode_angle)
      !! Takes the stress tensor sig and returns invariants p, q, and lode_angle
      real(dp), dimension(6), intent(in) :: sig
      real(dp), intent(out) :: p, q, lode_angle

      ! Local variables
      real(dp) :: dev(6), J2, J3

      p = calc_mean_stress(sig)

      dev = calc_dev_stress(sig, p)

      J2 = calc_J2(dev)

      q = calc_q(J2) ! deviatoric stress

      !J3 stress invariant
      J3 = calc_J3(dev)

      lode_angle = calc_lode_angle(J2, J3) !Lode's angle

   end subroutine calc_sig_invariants

end module mod_stress_invariants
