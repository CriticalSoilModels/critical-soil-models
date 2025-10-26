! Module for holding the platic potential unit tests

module mod_test_plastic_potential_suite
    use stdlib_kinds, only: dp
    use ieee_arithmetic, only: ieee_is_nan
    use mod_stress_invariants, only: calc_q_invariant, calc_mean_stress, calc_J2_invariant
    use mod_voight_funcs, only: calc_dev_stress 
    use mod_stress_invar_deriv, only: calc_dq_to_dSigma, calc_mean_stress_to_dSigma
    use mod_plastic_potential, only: Get_dP_to_dSigma
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
        real(dp) :: m_vec(6), exp_m_vec(6)
        real(dp) :: dilatancy, stress(6), dev(6), q, mean_stress, J2
        logical :: passed
        real(dp), parameter :: tol = 1e-9_dp
        stress = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, 11.0_dp, 13.0_dp]

        dilatancy = -0.1
        
        ! Calc the value
        call Get_dP_to_dSigma(dilatancy, stress, m_vec)
        
        mean_stress = calc_mean_stress(stress)
        dev = calc_dev_stress(stress, mean_stress)
        J2 = calc_J2_invariant(dev)
        q = calc_q_invariant(J2)

        ! Set the expected value from a hand calc
        exp_m_vec =  -dilatancy * calc_mean_stress_to_dSigma() + calc_dq_to_dSigma(dev, q)

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