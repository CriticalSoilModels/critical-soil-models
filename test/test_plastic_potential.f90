! Module for holding the platic potential unit tests

module mod_test_plastic_potential_suite
    use mod_csm_kinds, only: wp
    use ieee_arithmetic, only: ieee_is_nan
    use mod_stress_invariants, only: calc_q_inv, calc_p_inv, calc_J2_inv
    use mod_voigt_utils, only: calc_dev_stress
    use mod_stress_invar_deriv, only: calc_dq_by_dsig, calc_dp_by_dsig
    use mod_plastic_potential, only: calc_dg_plas_by_dsig
    use mod_tensor_value_checker, only: check_tensor_values

    use testdrive, only : new_unittest, unittest_type, error_type, check

    implicit none

    private
    public :: collect_plastic_potential_suite
contains
    
    subroutine collect_plastic_potential_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        ! Mkae the test suite using the subroutine defined below
        testsuite = [&
            new_unittest("dP/dSigma", test_dP_to_dSigma) &
        ]
    end subroutine collect_plastic_potential_suite

    subroutine test_dP_to_dSigma(error)
        type(error_type), allocatable, intent(out) :: error

        ! Local variables
        real(wp) :: m_vec(6), exp_m_vec(6)
        real(wp) :: dilatancy, stress(6), dev(6), q, mean_stress, J2
        logical :: passed
        real(wp), parameter :: tol = 1e-9_wp
        stress = [1.0_wp, 3.0_wp, 5.0_wp, 7.0_wp, 11.0_wp, 13.0_wp]

        dilatancy = -0.1
        
        ! Calc the value
        call calc_dg_plas_by_dsig(dilatancy, stress, m_vec)
        
        mean_stress = calc_p_inv(stress)
        dev = calc_dev_stress(stress, mean_stress)
        J2 = calc_J2_inv(dev)
        q = calc_q_inv(J2)

        ! Set the expected value from a hand calc
        exp_m_vec =  -dilatancy * calc_dp_by_dsig() + calc_dq_by_dsig(dev, q)

        if (all(ieee_is_nan(m_vec) .eqv. .False.)) then
            call check_tensor_values(m_vec, exp_m_vec, tol, passed)
        else
            ! Force the error
            passed = .False.
            print *, "There is a NaN in the result"
            print *, m_vec
        end if

        call check(error, passed, .True., more = "dPlasPoten/dSigma")

        if(allocated(error)) return

    end subroutine test_dP_to_dSigma
end module mod_test_plastic_potential_suite