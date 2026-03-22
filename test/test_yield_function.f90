module mod_test_yield_function
    use stdlib_kinds, only: dp, i32 => int32
    use mod_yield_function, only : calc_dF_by_dsig
    use mod_tensor_value_checker, only: check_tensor_values
    use ieee_arithmetic, only: ieee_is_nan

    use testdrive, only : new_unittest, unittest_type, error_type, check

    implicit none

    private
    public :: collect_yield_function_suite

contains

    subroutine collect_yield_function_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [ &
            new_unittest("dF/dSigma no NaN", test_dF_by_dsig) &
        ]

    end subroutine collect_yield_function_suite

    subroutine test_dF_by_dsig(error)
        ! Smoke test: verify calc_dF_by_dsig produces no NaN values
        ! TODO: add hand-calculated expected values
        type(error_type), allocatable, intent(out) :: error

        real(dp) :: M_tc, eta_y, sig(6), n_vec(6)
        logical :: no_nans

        M_tc  = 1.0_dp
        eta_y = 1.5_dp
        sig   = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, 11.0_dp, 13.0_dp]

        call calc_dF_by_dsig(M_tc, eta_y, sig, n_vec)

        no_nans = all(ieee_is_nan(n_vec) .eqv. .false.)

        call check(error, no_nans, .true., more = "calc_dF_by_dsig produced NaN values")
        if (allocated(error)) return

    end subroutine test_dF_by_dsig

end module mod_test_yield_function
