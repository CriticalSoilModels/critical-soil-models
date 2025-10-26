module mod_test_stress_invariants_suite
   ! Local imports
   use kind_precision_module, only : dp, i32
   use mod_strain_invariants, only : Get_strain_invariants
   use mod_stress_invariants, only : Get_invariants, calc_theta_s, calc_J2_invariant, &
                                     calc_inc_driver_J3_invariant, calc_mean_stress
   use stdlib_linalg        , only : det
   use mod_voigt_functions  , only : inc_driver_voigt_2_matrix, calc_dev_stess

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
      real(kind = dp) :: exp_p, p
      real(kind = dp) :: stress(6)

      call random_number(stress)
      
      p = calc_mean_stress(stress)

      exp_p = sum(stress(1:3)) / 3.0_dp

      call check(error, p, exp_p, more = "Indiv. Mean Stress test")
      if(allocated(error)) return

   end subroutine test_mean_stress

   subroutine test_J2_invariant(error)
      type(error_type), allocatable, intent(out) :: error

      ! Local variables
      real(kind = dp) :: exp_J2, J2, stress(6), dev_stress(6)
      real(kind = dp) :: mean_stress

      call random_number(stress)
      mean_stress = calc_mean_stress(stress)

      dev_stress = calc_dev_stess(stress, mean_stress)
      
      J2 = calc_J2_invariant(dev_stress)

      ! Square all the terms
      dev_stress = dev_stress**2

      ! Calc J2 (Need to double the shear terms)
      exp_J2 = 0.5_dp * ( sum(dev_stress(1:3)) + sum(2 * dev_stress(4:6)) )

      call check(error, J2, exp_J2, more = "Indiv. J2 test")
      if(allocated(error)) return

   end subroutine test_J2_invariant

   subroutine test_J3_invariant(error)
      type(error_type), allocatable, intent(out) :: error

      ! Local variables
      real(kind = dp) :: J3, exp_J3, stress(6), dev_stress(6), dev_matrix(3,3)

      call random_number(stress)

      dev_stress = calc_dev_stess( stress, calc_mean_stress(stress) )

      ! Form the full matrix
      dev_matrix = inc_driver_voigt_2_matrix(dev_stress)

      ! Calc the determinant of the matrix
      exp_J3 = det(dev_matrix)

      J3 = calc_inc_driver_J3_invariant(dev_stress)

      call check(error, J3, exp_J3, more = "Indiv. J3 test")
      if(allocated(error)) return

   end subroutine test_J3_invariant

   subroutine test_lode_angle(error)
      ! Testing the Lode angle calculation for different stress conditions
      type(error_type), allocatable, intent(out) :: error
   
      ! Define stress states and the value of PI
      real(kind = dp), parameter :: &
         trx_compression(6) = [77.0_dp, 16.0_dp, 16.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
         trx_extension(6)   = [5.0_dp, 5.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
         shear(6)           = [1.0_dp, 2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
         PI = 4.0_dp * atan(1.0_dp)   ! Value of PI
   
      real(kind = dp) :: compress_lode, exten_lode, shear_lode
      real(kind = dp) :: dev(6), mean_stress
   
      ! Test for triaxial compression (expected Lode angle is -pi/6)
      mean_stress = calc_mean_stress(trx_compression)
      dev = calc_dev_stess(trx_compression, mean_stress)
      compress_lode = calc_theta_s(calc_J2_invariant(dev), calc_inc_driver_J3_invariant(dev))

      call check(error, compress_lode, -PI / 6.0_dp, more = "compression lode angle")
      if (allocated(error)) return
   
      ! Test for triaxial extension (expected Lode angle is pi/6)
      mean_stress = calc_mean_stress(trx_extension)
      dev = calc_dev_stess(trx_extension, mean_stress)
      exten_lode = calc_theta_s(calc_J2_invariant(dev), calc_inc_driver_J3_invariant(dev))

      call check(error, exten_lode, PI / 6.0_dp, more = "extension lode angle")
      if (allocated(error)) return
   
      ! Test for shear stress (expected Lode angle is 0)
      mean_stress = calc_mean_stress(shear)
      dev = calc_dev_stess(shear, mean_stress)
      shear_lode = calc_theta_s(calc_J2_invariant(dev), calc_inc_driver_J3_invariant(dev))

      call check(error, shear_lode, 0.0_dp, more = "shear lode angle")
      if (allocated(error)) return
   
   end subroutine test_lode_angle

   ! Test the Get_invariants functions. This is some what redundant due to the above tests but it makes it easier to test each of them individually
   subroutine test_main_stress_invariant(error)
      type(error_type), allocatable, intent(out) :: error

      ! Local varaibles
      real(kind = dp) :: stress(6), p, q, theta
      real(kind = dp), parameter :: &
         exp_p     = 3.0_dp , &
         exp_q     = 32.07802986469088_dp

      real(kind = dp) :: exp_theta, dev(6), J3, J2

      stress = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, 11.0_dp, 13.0_dp]

      ! Calc the stress invariants
      call Get_invariants(stress, p, q, theta)

      ! Check the mean stress
      call check(error, p, exp_p, more = "Mean Stress Test^")
      if (allocated(error)) return

      ! Check the dev. stress, (q, von mises stress)
      call check(error, q, exp_q, more = "Dev. Stress Test^")
      if(allocated(error)) return

      ! Check the lode angle
      ! See function calc_theta_s for more information about the lode angle used here
   
      dev = calc_dev_stess(stress, p)
      J3 = calc_inc_driver_J3_invariant(dev)
      J2 = calc_J2_invariant(dev)

      ! Formula from Potts and Zdravković
      exp_theta = - asin(3.0_dp * sqrt(3.0_dp) / 2.0_dp * J3 / J2**1.5_dp) / 3.0_dp
      
      call check(error, theta, exp_theta)
      if(allocated(error)) return

   end subroutine test_main_stress_invariant

end module mod_test_stress_invariants_suite
