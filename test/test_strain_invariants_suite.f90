module mod_test_strain_invariants_suite
    ! Local imports
    use mod_csm_kinds, only: wp
    use mod_strain_invariants, only: calc_eps_inv, calc_eps_q_inv, &
                                     calc_eps_vol_inv, calc_dev_strain

    ! Testdrive imports
    use testdrive, only : new_unittest, unittest_type, error_type, check

    implicit none

    private
    public :: collect_strain_invariants_suite

contains

    subroutine collect_strain_invariants_suite(testsuite)
       type(unittest_type), allocatable, intent(out) :: testsuite(:)

       testsuite = [ &
          new_unittest("strain_invariant_check",        test_strain_invariant),    &
          new_unittest("eps_q_matches_reference_impl",  test_eps_q_vs_reference)   &
          ]

    end subroutine collect_strain_invariants_suite

    subroutine test_strain_invariant(error)
       type(error_type), allocatable, intent(out) :: error

       real(kind=wp) :: strain(6), eps_v, eps_q
       real(kind=wp), parameter :: exp_eps_v = 9.0_wp,               &
          exp_eps_q = 10.878112581387148_wp,                          &
          tol       = epsilon(eps_v)

       strain = [1.0_wp, 3.0_wp, 5.0_wp, 7.0_wp, 11.0_wp, 13.0_wp]

       call calc_eps_inv(strain, eps_v, eps_q)

       call check(error, eps_v, exp_eps_v, more="Eps Volumetric test")
       if (allocated(error)) return

       call check(error, eps_q, exp_eps_q, more="Eps q test")

    end subroutine test_strain_invariant

    subroutine test_eps_q_vs_reference(error)
       !! Verify calc_eps_q_inv against an independent reference formula
       !! (Python implementation, ported to Fortran):
       !!
       !!   eps_q = 1/3 * sqrt( 2*((e22-e33)^2 + (e33-e11)^2 + (e11-e22)^2)
       !!                       + 3*(e12^2 + e13^2 + e23^2) )
       !!
       !! Algebraically equivalent to sqrt(2/3 * dev:dev) for engineering shear
       !! strains in Voigt notation. Voigt order: [11,22,33,12,13,23].
       type(error_type), allocatable, intent(out) :: error

       real(kind=wp) :: strain(6), eps_v, dev(6)
       real(kind=wp) :: eps_q_lib, eps_q_ref
       real(kind=wp), parameter :: tol = 100.0_wp * epsilon(eps_q_lib)

       strain = [1.0_wp, 3.0_wp, 5.0_wp, 7.0_wp, 11.0_wp, 13.0_wp]

       eps_v     = calc_eps_vol_inv(strain)
       dev       = calc_dev_strain(strain, eps_v)
       eps_q_lib = calc_eps_q_inv(dev)
       eps_q_ref = calc_eps_q_reference(strain)

       call check(error, abs(eps_q_lib - eps_q_ref) < tol, .true., &
          more="calc_eps_q_inv must match reference formula")

    end subroutine test_eps_q_vs_reference

    pure function calc_eps_q_reference(strain) result(eps_q)
       !! Reference implementation ported from Python.
       !! strain(4:6) are engineering shear strains (Voigt, order 12,13,23).
       real(kind=wp), intent(in) :: strain(6)
       real(kind=wp) :: eps_q

       real(kind=wp) :: normal_terms, shear_term

       normal_terms = (strain(2) - strain(3))**2 &   ! (e22 - e33)^2
                    + (strain(3) - strain(1))**2 &   ! (e33 - e11)^2
                    + (strain(1) - strain(2))**2      ! (e11 - e22)^2

       shear_term   = strain(4)**2 + strain(5)**2 + strain(6)**2

       eps_q = (1.0_wp / 3.0_wp) * sqrt(2.0_wp * normal_terms + 3.0_wp * shear_term)

    end function calc_eps_q_reference

end module mod_test_strain_invariants_suite
