module mod_test_stress_invariants_suite
   ! Local imports
   use mod_csm_kinds, only: wp
   use mod_strain_invariants, only : calc_eps_inv
   use mod_stress_invariants, only : calc_sig_inv, calc_lode_inv, calc_J2_inv, &
                                     calc_J3_inv, calc_p_inv
   use stdlib_linalg        , only : det
   use mod_voigt_utils  , only : calc_voigt_to_matrix, calc_dev_stress

   ! Testdrive imports
   use testdrive, only : new_unittest, unittest_type, error_type, check

   implicit none

   private
   public :: collect_stress_invariants_suite
   ! Note the convention used in incremental driver is
contains


   subroutine collect_stress_invariants_suite(testsuite)
      ! Collection of tests
      ! Inidividual tests are stored in unitest_type

      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      ! Make the test suite using the subroutines defined below
      testsuite = [ &
         new_unittest("mean_stress_check"     , test_mean_stress), &
         new_unittest("J2_check"              , test_J2_invariant) , &
         new_unittest("J3 check"              , test_J3_invariant) , &
         new_unittest("lode_angle_check"      , test_lode_angle )     , &
         new_unittest("stress_invariant_check", test_main_stress_invariant)  &

         ]

   end subroutine collect_stress_invariants_suite
   
   subroutine test_mean_stress(error)
      type(error_type), allocatable, intent(out) :: error

      ! Local variables
      real(kind=wp) :: exp_p, p
      real(kind=wp) :: stress(6)

      stress = [100.0_wp, 50.0_wp, 25.0_wp, 10.0_wp, 5.0_wp, 2.0_wp]

      p = calc_p_inv(stress)

      exp_p = sum(stress(1:3)) / 3.0_wp

      call check(error, p, exp_p, more = "Indiv. Mean Stress test")
      if(allocated(error)) return

   end subroutine test_mean_stress

   subroutine test_J2_invariant(error)
      type(error_type), allocatable, intent(out) :: error

      ! Local variables
      real(kind=wp) :: exp_J2, J2, stress(6), dev_stress(6)
      real(kind=wp) :: mean_stress

      stress = [100.0_wp, 50.0_wp, 25.0_wp, 10.0_wp, 5.0_wp, 2.0_wp]
      mean_stress = calc_p_inv(stress)

      dev_stress = calc_dev_stress(stress, mean_stress)
      
      J2 = calc_J2_inv(dev_stress)

      ! Square all the terms
      dev_stress = dev_stress**2

      ! Calc J2 (Need to double the shear terms)
      exp_J2 = 0.5_wp * ( sum(dev_stress(1:3)) + sum(2 * dev_stress(4:6)) )

      call check(error, J2, exp_J2, more = "Indiv. J2 test")
      if(allocated(error)) return

   end subroutine test_J2_invariant

   subroutine test_J3_invariant(error)
      type(error_type), allocatable, intent(out) :: error

      ! Local variables
      real(kind=wp) :: J3, exp_J3, stress(6), dev_stress(6), dev_matrix(3,3)

      stress = [100.0_wp, 50.0_wp, 25.0_wp, 10.0_wp, 5.0_wp, 2.0_wp]

      dev_stress = calc_dev_stress( stress, calc_p_inv(stress) )

      ! Form the full matrix
      dev_matrix = calc_voigt_to_matrix(dev_stress)

      ! Calc the determinant of the matrix
      exp_J3 = det(dev_matrix)

      J3 = calc_J3_inv(dev_stress)

      call check(error, J3, exp_J3, more = "Indiv. J3 test")
      if(allocated(error)) return

   end subroutine test_J3_invariant

   subroutine test_lode_angle(error)
      ! Testing the Lode angle calculation for different stress conditions
      type(error_type), allocatable, intent(out) :: error
   
      ! Define stress states and the value of PI
      real(kind=wp), parameter :: &
         trx_compression(6) = [77.0_wp, 16.0_wp, 16.0_wp, 0.0_wp, 0.0_wp, 0.0_wp], &
         trx_extension(6)   = [5.0_wp, 5.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp], &
         shear(6)           = [1.0_wp, 2.0_wp, 3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp], &
         PI = 4.0_wp * atan(1.0_wp)   ! Value of PI
   
      real(kind=wp) :: compress_lode, exten_lode, shear_lode
      real(kind=wp) :: dev(6), mean_stress
   
      ! Test for triaxial compression (expected Lode angle is -pi/6)
      mean_stress = calc_p_inv(trx_compression)
      dev = calc_dev_stress(trx_compression, mean_stress)
      compress_lode = calc_lode_inv(calc_J2_inv(dev), calc_J3_inv(dev))

      call check(error, compress_lode, -PI / 6.0_wp, more = "compression lode angle")
      if (allocated(error)) return
   
      ! Test for triaxial extension (expected Lode angle is pi/6)
      mean_stress = calc_p_inv(trx_extension)
      dev = calc_dev_stress(trx_extension, mean_stress)
      exten_lode = calc_lode_inv(calc_J2_inv(dev), calc_J3_inv(dev))

      call check(error, exten_lode, PI / 6.0_wp, more = "extension lode angle")
      if (allocated(error)) return
   
      ! Test for shear stress (expected Lode angle is 0)
      mean_stress = calc_p_inv(shear)
      dev = calc_dev_stress(shear, mean_stress)
      shear_lode = calc_lode_inv(calc_J2_inv(dev), calc_J3_inv(dev))

      call check(error, shear_lode, 0.0_wp, more = "shear lode angle")
      if (allocated(error)) return
   
   end subroutine test_lode_angle

   ! Test the calc_sig_inv functions. This is some what redundant due to the above tests but it makes it easier to test each of them individually
   subroutine test_main_stress_invariant(error)
      type(error_type), allocatable, intent(out) :: error

      ! Local varaibles
      real(kind=wp) :: stress(6), p, q, lode_angle
      real(kind=wp), parameter :: &
         exp_p     = 3.0_wp , &
         exp_q     = 32.07802986469088_wp

      real(kind=wp) :: exp_lode_angle, dev(6), J3, J2

      stress = [1.0_wp, 3.0_wp, 5.0_wp, 7.0_wp, 11.0_wp, 13.0_wp]

      ! Calc the stress invariants
      call calc_sig_inv(stress, p, q, lode_angle)

      ! Check the mean stress
      call check(error, p, exp_p, more = "Mean Stress Test^")
      if (allocated(error)) return

      ! Check the dev. stress, (q, von mises stress)
      call check(error, q, exp_q, more = "Dev. Stress Test^")
      if(allocated(error)) return

      ! Check the lode angle
      ! See function calc_lode_inv for more information about the lode angle used here

      dev = calc_dev_stress(stress, p)
      J3 = calc_J3_inv(dev)
      J2 = calc_J2_inv(dev)

      ! Formula from Potts and Zdravković
      exp_lode_angle = - asin(3.0_wp * sqrt(3.0_wp) / 2.0_wp * J3 / J2**1.5_wp) / 3.0_wp

      call check(error, lode_angle, exp_lode_angle)
      if(allocated(error)) return

   end subroutine test_main_stress_invariant

end module mod_test_stress_invariants_suite
