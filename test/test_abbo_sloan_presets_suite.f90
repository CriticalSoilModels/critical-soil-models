module mod_test_abbo_sloan_presets_suite
   use mod_csm_kinds,          only: wp
   use mod_mcss_types,         only: abbo_sloan_params_t, DEFAULT_SMOOTH_COEFF
   use mod_abbo_sloan_presets, only: abbo_sloan_preset
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none
   private
   public :: collect_abbo_sloan_presets_suite

contains

   subroutine collect_abbo_sloan_presets_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("as_preset: 25.0 deg constants",  test_preset_25),   &
         new_unittest("as_preset: 29.5 deg constants",  test_preset_295),  &
         new_unittest("as_preset: unknown -> 29.5 deg", test_preset_default) &
         ]
   end subroutine collect_abbo_sloan_presets_suite

   subroutine test_preset_25(error)
      type(error_type), allocatable, intent(out) :: error
      type(abbo_sloan_params_t) :: asp
      real(wp), parameter :: tol = 1.0e-12_wp

      asp = abbo_sloan_preset(25.0_wp)

      call check(error, abs(asp%lode_tr - 0.436332312998582_wp) < tol, .True., &
                 more = "25 deg: lode_tr mismatch")
      if (allocated(error)) return
      call check(error, abs(asp%A1 - 1.432052062044227_wp) < tol, .True., &
                 more = "25 deg: A1 mismatch")
      if (allocated(error)) return
      call check(error, abs(asp%A2 - 0.406941858374615_wp) < tol, .True., &
                 more = "25 deg: A2 mismatch")
      if (allocated(error)) return
      call check(error, abs(asp%B1 - 0.544290524902313_wp) < tol, .True., &
                 more = "25 deg: B1 mismatch")
      if (allocated(error)) return
      call check(error, abs(asp%B2 - 0.673903324498392_wp) < tol, .True., &
                 more = "25 deg: B2 mismatch")
      if (allocated(error)) return
      call check(error, abs(asp%smooth_coeff - DEFAULT_SMOOTH_COEFF) < tol, .True., &
                 more = "25 deg: smooth_coeff mismatch")
   end subroutine test_preset_25

   subroutine test_preset_295(error)
      type(error_type), allocatable, intent(out) :: error
      type(abbo_sloan_params_t) :: asp
      real(wp), parameter :: tol = 1.0e-12_wp

      asp = abbo_sloan_preset(29.5_wp)

      call check(error, abs(asp%lode_tr - 0.514872129338327_wp) < tol, .True., &
                 more = "29.5 deg: lode_tr mismatch")
      if (allocated(error)) return
      call check(error, abs(asp%A1 - 7.138654723242414_wp) < tol, .True., &
                 more = "29.5 deg: A1 mismatch")
      if (allocated(error)) return
      call check(error, abs(asp%A2 - 6.112267270920612_wp) < tol, .True., &
                 more = "29.5 deg: A2 mismatch")
      if (allocated(error)) return
      call check(error, abs(asp%B1 - 6.270447753139589_wp) < tol, .True., &
                 more = "29.5 deg: B1 mismatch")
      if (allocated(error)) return
      call check(error, abs(asp%B2 - 6.398760841429403_wp) < tol, .True., &
                 more = "29.5 deg: B2 mismatch")
   end subroutine test_preset_295

   subroutine test_preset_default(error)
      !! An unrecognised angle (e.g. 27.0°) should fall back to the 29.5° set.
      type(error_type), allocatable, intent(out) :: error
      type(abbo_sloan_params_t) :: asp
      real(wp), parameter :: tol = 1.0e-12_wp

      asp = abbo_sloan_preset(27.0_wp)

      call check(error, abs(asp%lode_tr - 0.514872129338327_wp) < tol, .True., &
                 more = "default: lode_tr should match 29.5 deg set")
      if (allocated(error)) return
      call check(error, abs(asp%A1 - 7.138654723242414_wp) < tol, .True., &
                 more = "default: A1 should match 29.5 deg set")
   end subroutine test_preset_default

end module mod_test_abbo_sloan_presets_suite
