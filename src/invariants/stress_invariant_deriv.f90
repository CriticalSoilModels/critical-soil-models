! Module holds the derivatives of the stress invariants

module mod_stress_invar_deriv
   use mod_csm_kinds, only: wp
   use mod_voigt_utils  , only : calc_voigt_square

   implicit none

contains

   pure function calc_dp_by_dsig() result(dp_by_dsig)
      !! dp/dsig = dp/ds = [1/3, 1/3, 1/3, 0, 0, 0].
      !! Mean stress p = tr(σ)/3 depends only on normal components, so its
      !! derivative is the same whether taken w.r.t. full stress σ or
      !! deviatoric stress s.
      real(wp) :: dp_by_dsig(6)

      ! Local variables
      integer :: i

      ! Zero all the values
      dp_by_dsig(:) = 0.0_wp

      ! Add the 1/3 to the first three elements
      do i = 1, 3
         dp_by_dsig(i) = 1.0_wp / 3.0_wp
      end do

   end function calc_dp_by_dsig

   pure function calc_dq_by_dsig(dev, q) result(dq_by_dsig)
      ! Calc dq/dSigma
      real(wp), intent(in) :: dev(6), q
      real(wp) :: dq_by_dsig(6)

      ! Local variables
      real(wp) :: dJ2_by_dsig(6)

      ! Calc dJ2/dsig
      dJ2_by_dsig = calc_dJ2_by_dsig(dev)

      ! Calc dq/dsig
      dq_by_dsig = 3.0_wp/ (2.0_wp * q) * dJ2_by_dsig

   end function calc_dq_by_dsig

   pure function calc_dJ2_by_dsig(dev) result(dJ2_by_dsig)
      real(wp), intent(in) :: dev(6)
      real(wp) :: dJ2_by_dsig(6)

      dJ2_by_dsig = dev

      ! Double the shear terms
      dJ2_by_dsig(4:6) = 2.0_wp * dJ2_by_dsig(4:6)

   end function calc_dJ2_by_dsig

   pure function calc_dJ_by_dsig(dev, J) result(dJ_by_dsig)
      !! dJ/dσ = dJ2/dσ / (2J), where J = sqrt(J2).
      !! Caller must guard against J ≈ 0 to avoid division by zero.
      real(wp), intent(in) :: dev(6), J
      real(wp) :: dJ_by_dsig(6)
      dJ_by_dsig = calc_dJ2_by_dsig(dev) / (2.0_wp * J)
   end function calc_dJ_by_dsig

   pure function calc_dJ3_by_dsig(dev) result(dJ3_by_dsig)
      real(wp), intent(in) :: dev(6)
      real(wp) :: dJ3_by_dsig(6)

      ! Local variables
      real(wp) :: ii(6), dev2(6), tr_s2

      !Fill S.S
      dev2 = calc_voigt_square( dev )

      !Compute dJ3_by_dsig
      tr_s2 = dev2(1) + dev2(2) + dev2(3)

      ii = 0.0_wp   ! Identity tensor
      ii(1) = 1.0_wp
      ii(2) = 1.0_wp
      ii(3) = 1.0_wp

      ! This is equaivalent to s^{2} - 2/3 J_{2} \matr{1}
      ! J_{2}(\matr{s}) = 2 * I_{1}(\matr{ s^{2} })
      ! See Appendix B. Invariant Notes Moore, Jonathan Thesis for more details
      dJ3_by_dsig = dev2 - ( tr_s2*ii / 3.0_wp )

      ! Need to double the shear terms because voigt notation is being used and therefore shear terms are linked together
      dJ3_by_dsig(4:6) = 2.0_wp * dJ3_by_dsig(4:6)

   end function calc_dJ3_by_dsig

   pure function calc_dlode_angle_by_dsig(dJ3_by_dsig, dev, J3, J2, lode_angle) result(dlode_angle_by_dsig)
      real(wp), intent(in) :: dJ3_by_dsig(6), dev(6)
      real(wp), intent(in) :: J3, J2, lode_angle
      real(wp) :: dlode_angle_by_dsig(6)

      ! Local variables
      real(wp) :: cos_term, outside_term, inside_term_1(6), dJ2_by_dsig(6)
      real(wp), parameter :: tolerance = 1e-12_wp
      real(wp), parameter :: THREE = 3.0_wp, &
         TWO   = 2.0_wp, &
         ZERO  = 0.0_wp
      ! Calc cos(3 * lode_angle)
      cos_term = cos(3.0_wp * lode_angle)

      ! TODO: Turn the tolerance back on once I've checked that they match
      ! ! If cos term is zero( Trx compression or tension) set to tiny value
      ! if ( abs(cos_term) <= tolerance) then
      !    cos_term = tolerance
      ! end if

      ! Calc the fraction before the parenthesis
      outside_term = sqrt(THREE) / (TWO * cos_term * J2**1.5_wp)

      ! Calc the first term inside the parenthesis
      dJ2_by_dsig = calc_dJ2_by_dsig(dev)

      inside_term_1 = THREE * J3 / (2 * J2) * dJ2_by_dsig

      ! Calc the full term
      dlode_angle_by_dsig = outside_term * (inside_term_1 - dJ3_by_dsig)

   end function calc_dlode_angle_by_dsig

end module mod_stress_invar_deriv
