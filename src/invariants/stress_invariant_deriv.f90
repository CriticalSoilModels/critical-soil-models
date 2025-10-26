! Module holds the derivatives of the stress invariants

module mod_stress_invar_deriv
   use kind_precision_module, only: real_type => dp
   use mod_voigt_functions  , only : square_inc_driver_voigt_vector

   implicit none

contains

   pure function calc_mean_stress_to_dSigma() result(dmean_dSigma)
      real(kind = real_type) :: dmean_dSigma(6)

      ! Local variables
      integer :: i

      ! Zero all the values
      dmean_dSigma(:) = 0.0_real_type

      ! Add the 1/3 to the first three elements
      do i = 1, 3
         dmean_dSigma(i) = 1.0_real_type / 3.0_real_type
      end do

   end function calc_mean_stress_to_dSigma

   pure function calc_dq_to_dSigma(dev, q) result(dq_dSigma)
      ! Calc dq/dSigma
      real(kind = real_type), intent(in) :: dev(6), q
      real(kind = real_type) :: dq_dSigma(6)

      ! Local variables
      real(kind = real_type) :: dJ2_dsigma(6)

      ! Calc dJ2/dsigma
      dJ2_dSigma = calc_dJ2_to_dSigma(dev)

      ! Calc dq/dSigma
      dq_dSigma = 3.0_real_type/ (2.0_real_type * q) * dJ2_dsigma

   end function calc_dq_to_dSigma

   pure function calc_dJ2_to_dSigma(dev) result(dJ2_dSigma)
      real(kind = real_type), intent(in) :: dev(6)
      real(kind = real_type) :: dJ2_dSigma(6)

      dJ2_dSigma = dev

      ! Double the shear terms
      dJ2_dSigma(4:6) = 2.0_real_type * dJ2_dSigma(4:6)

   end function calc_dJ2_to_dSigma

   pure function calc_dJ3_to_dSigma(dev) result(dJ3dSig)
      real(kind = real_type), intent(in) :: dev(6)
      real(kind = real_type) :: dJ3dSig(6)

      ! Local variables
      real(kind=  real_type) :: II(6), dev2(6), TrS2

      !Fill S.S
      dev2 = square_inc_driver_voigt_vector( dev )

      !Compute dJ3dSig
      TrS2 = dev2(1) + dev2(2) + dev2(3)

      II=0.0d0!Identity tensor
      II(1)=1.0
      II(2)=1.0
      II(3)=1.0

      ! This is equaivalent to s^{2} - 2/3 J_{2} \matr{1}
      ! J_{2}(\matr{s}) = 2 * I_{1}(\matr{ s^{2} })
      ! See Appendix B. Invariant Notes Moore, Jonathan Thesis for more details
      dJ3dsig = dev2 - ( TrS2*II / 3.0d0 )

      ! Need to double the shear terms because voigt notation is being used and therefore shear terms are linked together
      dJ3dSig(4:6) = 2.0_real_type * dJ3dSig(4:6)

   end function calc_dJ3_to_dSigma

   pure function calc_inc_driver_dJ3_to_dSigma(stress) result(dJ3_dSigma)
      ! This function calculates dJ3/dSigma assuming that the input voigt vector follows incremental driver convention
      ! As such the output follows incremental driver convention as well

      real(kind = real_type), intent(in) :: stress(6)
      real(kind = real_type) :: dJ3_dSigma(6)

      ! Local variables
      real(kind=  real_type) :: t(6)
      real(kind = real_type),parameter :: ONE_NINTH = 1.0_real_type/9.0_real_type, &
         ONE_THIRD = 1.0_real_type/3.0_real_type, &
         TWO       = 2.0_real_type, &
         FOUR      = 4.0_real_type

      ! Store the stress in a shorter name to make it easier to type
      t = stress

      ! Calc dJ3/dSigma_{11}
      dJ3_dSigma(1) = ONE_NINTH * (TWO * t(1)**2 - t(2)**2 - t(3)**2 - &
         TWO * t(1) * t(2) - 2 * t(1) * t(3) + FOUR * t(2) * t(3)) + &
         ONE_THIRD * (t(4)**2 + t(5)**2 - TWO * t(6)**2)

      ! Calc dJ3/dSigma_{22}
      dJ3_dSigma(2) = ONE_NINTH * (-t(1)**2 + TWO * t(2)**2 - t(3)**2 - &
         TWO * t(1) * t(2) + FOUR * t(1) * t(3) - TWO * t(2) * t(3)) + &
         ONE_THIRD * (t(4)**2 - TWO * t(5)**2 + t(6)**2)

      ! Calc dJ3/dSigma_{33}
      dJ3_dSigma(3) = ONE_NINTH * (-t(1)**2 -t(2)**2 + TWO * t(3)**2 + &
         FOUR * t(1) * t(2) - TWO * t(1) * t(3) - TWO * t(2) * t(3)) + &
         ONE_THIRD * (-TWO * t(4)**2 + t(5)**2 + t(6)**2)

      ! Calc dJ3/dSigma_{12}
      dJ3_dSigma(4) = ONE_THIRD * (TWO * t(1) * t(4) + TWO * t(4) * t(2) - FOUR * t(4) * t(3)) + TWO * t(5) * t(6)

      ! Calc dJ3/dSigma_{13}
      dJ3_dSigma(5) = ONE_THIRD * (TWO * t(1) * t(5) - FOUR * t(5) * t(2) + TWO * t(5) * t(3)) + TWO * t(4) * t(6)

      ! Calc dJ3/dSigma_{23}
      dJ3_dSigma(6) = ONE_THIRD * (-FOUR * t(1) * t(6) + TWO * t(2) * t(6) + TWO * t(6) * t(3)) + TWO * t(4) * t(5)

   end function calc_inc_driver_dJ3_to_dSigma

   pure function calc_dtheta_to_dSigma(dJ3_dSigma, dev, J3, J2, theta) result(dtheta_dSigma)
      real(kind = real_type), intent(in) :: dJ3_dSigma(6), dev(6)
      real(kind = real_type), intent(in) :: J3, J2, theta
      real(kind = real_type) :: dtheta_dSigma(6)

      ! Local variables
      real(kind = real_type) :: cos_term, outside_term, inside_term_1(6), dJ2_dSigma(6)
      real(kind = real_type), parameter :: tolerance = 1e-12_real_type
      real(kind = real_type), parameter :: THREE = 3.0_real_type, &
         TWO   = 2.0_real_type, &
         ZERO  = 0.0_real_type
      ! Calc cos(3 \theta)
      cos_term = cos(3 * theta)

      ! TODO: Turn the tolerance back on once I've checked that they match
      ! ! If cos term is zero( Trx compression or tension) set to tiny value
      ! if ( abs(cos_term) <= tolerance) then
      !    cos_term = tolerance
      ! end if

      ! Calc the fraction before the parenthesis
      outside_term = sqrt(THREE) / (TWO * cos_term * J2**1.5_real_type)

      ! Calc the first term inside the parenthesis
      dJ2_dSigma = calc_dJ2_to_dSigma(dev)

      inside_term_1 = THREE * J3 / (2 * J2) * dJ2_dSigma

      ! Calc the full term
      dtheta_dSigma = outside_term * (inside_term_1 - dJ3_dSigma)

   end function calc_dtheta_to_dSigma

   pure function calc_dtheta_to_dSigma_2(dJ3_dSigma, dev, J3, J2) result(dtheta_dSigma)
      real(kind = real_type), intent(in) :: dJ3_dSigma(6), dev(6), J3, J2
      real(kind = real_type) :: dtheta_dSigma(6)

      ! Local variables
      real(kind = real_type) :: outside_term_1, outside_term_2, inside(6)
      real(kind = real_type) :: dJ2_dSigma(6)
      real(kind = real_type), parameter :: THREE = 3.0_real_type, &
                                           TWO = 2.0_real_type

      ! Calc the term on the outside of the paranthesis
      outside_term_1 = sqrt(THREE) / ( 2.0 * J2**(1.5) )

      outside_term_2 = 1/sqrt( 1 - (THREE * sqrt(THREE)/TWO * J3/J2**1.5)**2 )

       ! Calc dJ2/dSigma
      dJ2_dSigma = calc_dJ2_to_dSigma(dev)

      ! Calc the term inside 
      inside = THREE / TWO * J3/J2 * dJ2_dSigma - dJ3_dSigma
      
      ! Calc the total term
      dtheta_dSigma = outside_term_1 * outside_term_2 * inside
   end function calc_dtheta_to_dSigma_2

end module mod_stress_invar_deriv

