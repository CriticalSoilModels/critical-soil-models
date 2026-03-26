module mod_test_stress_invar_deriv_suite
    ! Imports
    use mod_csm_kinds, only: wp
    use mod_stress_invariants, only: calc_p_inv, calc_q_inv, calc_J2_inv, &
                                     calc_J3_inv, calc_lode_inv
    use mod_stress_invar_deriv, only: calc_dp_by_dsig, calc_dq_by_dsig, calc_dJ2_by_dsig, calc_dJ3_by_dsig, &
                                      calc_dlode_angle_by_dsig
    use mod_stress_invar_refs, only: calc_dJ3_by_dsig_full
    use mod_voigt_utils   , only: calc_dev_stress
    use mod_tensor_value_checker, only: check_tensor_values

    ! Testdrive imports
    use testdrive, only : new_unittest, unittest_type, error_type, check

    implicit none

    private
    public :: collect_stress_invar_deriv_suite
contains

    subroutine collect_stress_invar_deriv_suite(testsuite)
        type(unittest_type), allocatable, intent(out) ::  testsuite(:)

        ! Make the test suite using the subroutine defined below
        testsuite = [ &
            new_unittest("dp/dSigma" , test_dMean_Stress_to_dSigma),  &
            new_unittest("dq/dSigma" , test_dq_to_dSigma),  &
            new_unittest("dJ2/dSigma", test_dJ2_to_dSigma), &
            new_unittest("dJ3/dSigma", test_dJ3_to_dSigma), &
            new_unittest("dlode_angle_dSigma", test_dlode_angle_to_dSigma) &
            ]

    end subroutine collect_stress_invar_deriv_suite

    subroutine test_dMean_Stress_to_dSigma(error)
        type(error_type), allocatable, intent(out) :: error

        ! Local variables
        real(kind=wp) :: exp_dmean_dSigma(6), dmean_dSigma(6)
        real(kind=wp), parameter :: tol = 1e-9
        logical :: passed = .False.

        ! Calc the value using the function
        dmean_dSigma = calc_dp_by_dsig()

        ! Calc the value using another method
        exp_dmean_dSigma(:) = 0.0
        exp_dmean_dSigma(1:3) = 1.0_wp / 3.0_wp

        ! Check that all the values are the same within a tolerance
        call check_tensor_values(dmean_dSigma, exp_dmean_dSigma, tol, passed)

        ! Check that the check passed
        call check(error, passed, .True., more = "Indiv. dP/dSigma test")
        if(allocated(error)) return

    end subroutine test_dMean_Stress_to_dSigma

    subroutine test_dq_to_dSigma(error)
        type(error_type), allocatable, intent(out) ::  error

        ! Local variables
        real(kind=wp) :: exp_dq_dSigma(6), dq_dSigma(6)
        real(kind=wp) :: stress(6), dev_stress(6)
        real(kind=wp) :: mean_stress, J2, q
        real(kind=wp), parameter :: tol = 1e-9
        logical :: passed = .False.

        stress = [100.0_wp, 50.0_wp, 25.0_wp, 10.0_wp, 5.0_wp, 2.0_wp]
        mean_stress = calc_p_inv(stress)

        dev_stress = calc_dev_stress(stress, mean_stress)

        J2 = calc_J2_inv(dev_stress)

        q = calc_q_inv(J2)

        ! Calc the value using the function
        dq_dSigma = calc_dq_by_dsig(dev_stress, q)

        ! Calc the value using another method
        exp_dq_dSigma = 3.0_wp/(2.0_wp * q) * dev_stress

        ! Double the shear terms
        exp_dq_dSigma(4:6) = 2.0_wp * exp_dq_dSigma(4:6)

        ! Check that all the values are the same within a tolerance
        call check_tensor_values(dq_dSigma, exp_dq_dSigma, tol, passed)

        ! Check that the check passed
        call check(error, passed, .True., more = "Indiv. dq/dSigma test")
        if(allocated(error)) return

    end subroutine test_dq_to_dSigma

    subroutine test_dJ2_to_dSigma(error)
        type(error_type), allocatable, intent(out) ::  error

        ! Local variables
        real(kind=wp) :: exp_dJ2_dSigma(6), dJ2_dSigma(6)
        real(kind=wp) :: stress(6), dev_stress(6), mean_stress
        real(kind=wp), parameter :: tol = 1e-9
        logical :: passed = .False.

        stress = [100.0_wp, 50.0_wp, 25.0_wp, 10.0_wp, 5.0_wp, 2.0_wp]
        mean_stress = calc_p_inv(stress)

        dev_stress = calc_dev_stress(stress, mean_stress)

        ! Calc the value using the function
        dJ2_dSigma = calc_dJ2_by_dsig(dev_stress)

        ! Calc the value using another method
        exp_dJ2_dSigma(1) = dev_stress(1)
        exp_dJ2_dSigma(2) = dev_stress(2)
        exp_dJ2_dSigma(3) = dev_stress(3)
        exp_dJ2_dSigma(4) = 2.0_wp * dev_stress(4)
        exp_dJ2_dSigma(5) = 2.0_wp * dev_stress(5)
        exp_dJ2_dSigma(6) = 2.0_wp * dev_stress(6)

        ! Check that all the values are the same within a tolerance
        call check_tensor_values(dJ2_dSigma, exp_dJ2_dSigma, tol, passed)

        ! Check that the check passed
        call check(error, passed, .True., more = "Indiv. dJ2/dSigma test")
        if(allocated(error)) return

    end subroutine test_dJ2_to_dSigma

    subroutine test_dJ3_to_dSigma(error)
        type(error_type), allocatable, intent(out) ::  error

        ! Local variables
        real(kind=wp) :: exp_dJ3_dSigma(6), dJ3_dSigma(6)
        real(kind=wp) :: stress(6), dev(6), mean_stress
        real(kind=wp), parameter :: tol = 1e-9
        logical :: passed = .False.

        stress = [100.0_wp, 50.0_wp, 25.0_wp, 10.0_wp, 5.0_wp, 2.0_wp]
        mean_stress = calc_p_inv(stress)

        dev = calc_dev_stress(stress, mean_stress)

        ! Calc the value using the function
        dJ3_dSigma = calc_dJ3_by_dsig(dev)

        ! Calc the value using another method
        exp_dJ3_dSigma = calc_dJ3_by_dsig_full(stress)

        ! Check that all the values are the same within a tolerance
        call check_tensor_values(dJ3_dSigma, exp_dJ3_dSigma , tol, passed)

        ! Check that the check passed
        call check(error, passed, .True., more = "Indiv. dJ3/dSigma test")
        if(allocated(error)) return


    end subroutine test_dJ3_to_dSigma

    subroutine test_dlode_angle_to_dSigma(error)
        type(error_type), allocatable, intent(out) ::  error

        ! Local variables
        real(kind=wp) :: exp_dlode_angle_by_dsig(6), dlode_angle_by_dsig(6)
        real(kind=wp) :: stress(6), dev(6), mean_stress
        real(kind=wp) :: J2, J3, lode_angle, dJ3_dSigma(6)
        real(kind=wp), parameter :: tol = 1e-9
        logical :: passed = .False.

        stress = [100.0_wp, 50.0_wp, 25.0_wp, 10.0_wp, 5.0_wp, 2.0_wp]
        mean_stress = calc_p_inv(stress)

        dev = calc_dev_stress(stress, mean_stress)

        J2 = calc_J2_inv(dev)
        J3 = calc_J3_inv(dev)
        lode_angle = calc_lode_inv(J2, J3)
        dJ3_dSigma = calc_dJ3_by_dsig(dev)

        ! Calc the value using the function
        dlode_angle_by_dsig = calc_dlode_angle_by_dsig(dJ3_dSigma, dev, J3, J2, lode_angle)

        ! Calc the value using the local reference function
        exp_dlode_angle_by_dsig = dlode_angle_reference(dJ3_dSigma, dev, J3, J2)

        ! Check that all the values are the same within a tolerance
        call check_tensor_values(dlode_angle_by_dsig, exp_dlode_angle_by_dsig, tol, passed)

        ! Check that the check passed
        call check(error, passed, .True., more = "Indiv. dlode_angle/dSigma test")
        if(allocated(error)) return

    end subroutine test_dlode_angle_to_dSigma

    pure function dlode_angle_reference(dJ3_dSigma, dev, J3, J2) result(dlode_angle_by_dsig)
       real(wp), intent(in) :: dJ3_dSigma(6), dev(6), J3, J2
       real(wp) :: dlode_angle_by_dsig(6)
       real(wp) :: outside_term_1, outside_term_2, inside(6)
       real(wp) :: dJ2_dSigma(6)
       real(wp), parameter :: THREE = 3.0_wp, TWO = 2.0_wp
       outside_term_1 = sqrt(THREE) / ( 2.0 * J2**(1.5) )
       outside_term_2 = 1/sqrt( 1 - (THREE * sqrt(THREE)/TWO * J3/J2**1.5)**2 )
       dJ2_dSigma = calc_dJ2_by_dsig(dev)
       inside = THREE / TWO * J3/J2 * dJ2_dSigma - dJ3_dSigma
       dlode_angle_by_dsig = outside_term_1 * outside_term_2 * inside
    end function dlode_angle_reference

end module mod_test_stress_invar_deriv_suite
