! Stress invariant functions for use in constitutive models.
!
! Computation chain:
!   sig(6)
!     -> p          = calc_p_inv(sig)
!     -> dev        = calc_dev_stress(sig, p)      [mod_voigt_utils]
!     -> J2         = calc_J2_inv(dev)
!     -> J          = calc_J_inv(q)                ! = q / sqrt(3)
!     -> q          = calc_q_inv(J2)               ! = sqrt(3*J2)
!     -> J3         = calc_J3_inv(dev)
!     -> lode_angle = calc_lode_inv(J2, J3)
!
! All-in-one: calc_sig_inv(sig, p, q, lode_angle)

module mod_stress_invariants
   use mod_csm_kinds, only: wp
   use mod_voigt_utils, only: calc_two_norm_tensor, calc_dev_stress
   implicit none

contains

   pure function calc_p_inv(stress) result(p)
      !! Mean stress: p = tr(σ) / 3
      real(wp), intent(in) :: stress(6)
      real(wp) :: p

      p = (stress(1) + stress(2) + stress(3)) / 3.0_wp
   end function calc_p_inv

   pure function calc_J2_inv(dev) result(J2)
      !! Second deviatoric stress invariant: J2 = 0.5 * s:s
      real(wp), intent(in) :: dev(6)
      real(wp) :: J2

      J2 = 0.5_wp * calc_two_norm_tensor(dev)**2
   end function calc_J2_inv

   pure function calc_J_inv(q) result(J)
      !! Deviatoric stress measure J = q / sqrt(3).
      !! Used in Abbo-Sloan Mohr-Coulomb formulation where J = sqrt(J2).
      real(wp), intent(in) :: q
      real(wp) :: J

      J = q / sqrt(3.0_wp)
   end function calc_J_inv

   pure function calc_q_inv(J2) result(q)
      !! Deviatoric stress invariant: q = sqrt(3 * J2)
      real(wp), intent(in) :: J2
      real(wp) :: q

      q = sqrt(3.0_wp * J2)
   end function calc_q_inv

   pure function calc_J3_inv(dev) result(J3)
      !! Third deviatoric stress invariant: J3 = det(s)
      !! Voigt order: [11, 22, 33, 12, 13, 23]
      real(wp), intent(in) :: dev(6)
      real(wp) :: J3

      J3 = product(dev(1:3)) &
         - dev(1)*dev(6)**2  &   ! -s11 * s23^2
         - dev(2)*dev(5)**2  &   ! -s22 * s13^2
         - dev(3)*dev(4)**2  &   ! -s33 * s12^2
         + 2.0_wp * product(dev(4:6))  ! +2 * s12*s13*s23
   end function calc_J3_inv

   pure function calc_lode_inv(J2, J3) result(lode_angle)
      !! Lode angle θ ∈ [-π/6, π/6]
      !! θ = -π/6: triaxial compression (s1 >= s2 = s3)
      !! θ =    0: pure shear
      !! θ =  π/6: triaxial extension (s1 = s2 >= s3)
      !!
      !! Reference: Potts & Zdravković, FEA in Geotechnical Engineering, Vol 2, p.186
      real(wp), intent(in) :: J2, J3
      real(wp) :: lode_angle

      real(wp) :: sin3theta

      if (J2 > 0.0_wp) then
         sin3theta = 0.5_wp * J3 * (3.0_wp / J2)**1.5_wp
      else
         sin3theta = -1.0_wp
      end if

      sin3theta = max(-1.0_wp, min(1.0_wp, sin3theta))

      lode_angle = -asin(sin3theta) / 3.0_wp
   end function calc_lode_inv

   pure subroutine calc_sig_inv(sig, p, q, lode_angle)
      !! All-in-one: returns p, q, and lode_angle from stress vector sig(6).
      real(wp), intent(in)  :: sig(6)
      real(wp), intent(out) :: p, q, lode_angle

      real(wp) :: dev(6), J2, J3

      p   = calc_p_inv(sig)
      dev = calc_dev_stress(sig, p)
      J2  = calc_J2_inv(dev)
      q   = calc_q_inv(J2)
      J3  = calc_J3_inv(dev)
      lode_angle = calc_lode_inv(J2, J3)
   end subroutine calc_sig_inv

end module mod_stress_invariants
