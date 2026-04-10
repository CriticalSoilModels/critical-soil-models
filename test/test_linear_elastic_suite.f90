module mod_test_linear_elastic_suite
   use mod_csm_kinds,             only: wp
   use mod_linear_elastic_model,  only: linear_elastic_model_t, le_from_props
   use mod_euler_substep,         only: euler_substep, integrator_params_t
   use mod_tensor_value_checker,  only: check_tensor_values
   use testdrive, only: new_unittest, unittest_type, error_type, check

   implicit none
   private
   public :: collect_linear_elastic_suite

   ! G = 10000, nu = 0.25 => lame_1 = 3G = 30000, lame_2 = G = 10000
   real(wp), parameter :: G_TEST   = 10000.0_wp
   real(wp), parameter :: NU_TEST  = 0.25_wp
   real(wp), parameter :: FTOL     = 1.0e-9_wp
   real(wp), parameter :: STOL     = 1.0e-6_wp
   real(wp), parameter :: TOL      = 1.0e-10_wp

contains

   subroutine collect_linear_elastic_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("volumetric compression", test_volumetric), &
         new_unittest("pure shear",             test_pure_shear),  &
         new_unittest("uniaxial strain",        test_uniaxial)     &
      ]
   end subroutine collect_linear_elastic_suite

   subroutine test_volumetric(error)
      !! deps = e*[1,1,1,0,0,0] => dsig = 3K*e * [1,1,1,0,0,0]
      !! With nu=0.25: K = 5G/3, so 3K = 5G = 50000
      type(error_type), allocatable, intent(out) :: error

      type(linear_elastic_model_t) :: model
      real(wp) :: sig(6), deps(6), expected(6)
      logical  :: passed

      model = le_from_props([G_TEST, NU_TEST])

      sig  = 0.0_wp
      deps = [1.0e-4_wp, 1.0e-4_wp, 1.0e-4_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      ! 3K = 5G = 50000; dsig_ii = 50000 * 1e-4 = 5.0
      expected = [5.0_wp, 5.0_wp, 5.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      call euler_substep(model, sig, deps, &
                         integrator_params_t(ftol=FTOL, stol=STOL, dt_min=1.0e-9_wp, max_iters=100))

      call check_tensor_values(sig, expected, TOL, passed)
      call check(error, passed, .true., more="volumetric: sig /= D_e * deps")

   end subroutine test_volumetric

   subroutine test_pure_shear(error)
      !! deps = [0,0,0,g,0,0] => dsig = [0,0,0,G*g,0,0]
      type(error_type), allocatable, intent(out) :: error

      type(linear_elastic_model_t) :: model
      real(wp) :: sig(6), deps(6), expected(6)
      logical  :: passed

      model = le_from_props([G_TEST, NU_TEST])

      sig  = 0.0_wp
      deps = [0.0_wp, 0.0_wp, 0.0_wp, 1.0e-4_wp, 0.0_wp, 0.0_wp]

      ! dsig_12 = G * g = 10000 * 1e-4 = 1.0
      expected = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp]

      call euler_substep(model, sig, deps, &
                         integrator_params_t(ftol=FTOL, stol=STOL, dt_min=1.0e-9_wp, max_iters=100))

      call check_tensor_values(sig, expected, TOL, passed)
      call check(error, passed, .true., more="pure shear: sig /= D_e * deps")

   end subroutine test_pure_shear

   subroutine test_uniaxial(error)
      !! deps = [e,0,0,0,0,0] => dsig = [lame_1*e, lame_2*e, lame_2*e, 0,0,0]
      !! lame_1 = 3G = 30000, lame_2 = G = 10000
      type(error_type), allocatable, intent(out) :: error

      type(linear_elastic_model_t) :: model
      real(wp) :: sig(6), deps(6), expected(6)
      logical  :: passed

      model = le_from_props([G_TEST, NU_TEST])

      sig  = 0.0_wp
      deps = [1.0e-4_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      ! dsig = [30000, 10000, 10000, 0, 0, 0] * 1e-4 = [3.0, 1.0, 1.0, 0, 0, 0]
      expected = [3.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      call euler_substep(model, sig, deps, &
                         integrator_params_t(ftol=FTOL, stol=STOL, dt_min=1.0e-9_wp, max_iters=100))

      call check_tensor_values(sig, expected, TOL, passed)
      call check(error, passed, .true., more="uniaxial: sig /= D_e * deps")

   end subroutine test_uniaxial

end module mod_test_linear_elastic_suite
