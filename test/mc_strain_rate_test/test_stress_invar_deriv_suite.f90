module mod_test_stress_invar_deriv_suite
    ! Imports
    use kind_precision_module, only: dp, i32
    use mod_stress_invariants, only: calc_mean_stress, calc_q_invariant, calc_J2_invariant, &
                                     calc_inc_driver_J3_invariant, calc_theta_s
    use mod_stress_invar_deriv, only: calc_mean_stress_to_dSigma, calc_dq_to_dSigma, calc_dJ2_to_dSigma, calc_dJ3_to_dSigma, &
                                      calc_inc_driver_dJ3_to_dSigma, calc_dtheta_to_dSigma, calc_dtheta_to_dSigma_2
    use mod_voigt_functions   , only: calc_dev_stess
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
            new_unittest("dtheta_dSigma", test_dtheta_to_dSigma) &
            ]
        
    end subroutine collect_stress_invar_deriv_suite

    subroutine test_dMean_Stress_to_dSigma(error)
        type(error_type), allocatable, intent(out) :: error

        ! Local variables
        real(kind = dp) :: exp_dmean_dSigma(6), dmean_dSigma(6)
        real(kind = dp), parameter :: tol = 1e-9
        logical :: passed = .False.
        
        ! Calc the value using the function
        dmean_dSigma = calc_mean_stress_to_dSigma()
        
        ! Calc the value using another method
        exp_dmean_dSigma(:) = 0.0
        exp_dmean_dSigma(1:3) = 1.0_dp / 3.0_dp
        
        ! Check that all the values are the same within a tolerance
        call check_tensor_values(dmean_dSigma, exp_dmean_dSigma, tol, passed)
        
        ! Check that the check passed
        call check(error, passed, .True., more = "Indiv. dP/dSigma test")
        if(allocated(error)) return
        
    end subroutine test_dMean_Stress_to_dSigma

    subroutine test_dq_to_dSigma(error)
        type(error_type), allocatable, intent(out) ::  error
    
        ! Local variables
        real(kind = dp) :: exp_dq_dSigma(6), dq_dSigma(6)
        real(kind = dp) :: stress(6), dev_stress(6)
        real(kind = dp) :: mean_stress, J2, q
        real(kind = dp), parameter :: tol = 1e-9
        logical :: passed = .False.
        
        call random_number(stress)
        
        stress = stress - [0.5, 0.0, 0.5, 0.0, 0.5, 0.0]
        mean_stress = calc_mean_stress(stress)

        dev_stress = calc_dev_stess(stress, mean_stress)
        
        J2 = calc_J2_invariant(dev_stress)

        q = calc_q_invariant(J2)

        ! Calc the value using the function
        dq_dSigma = calc_dq_to_dSigma(dev_stress, q)

        ! Calc the value using another method
        exp_dq_dSigma = 3.0_dp/(2.0_dp * q) * dev_stress
        
        ! Double the shear terms
        exp_dq_dSigma(4:6) = 2.0_dp * exp_dq_dSigma(4:6) 
        
        ! Check that all the values are the same within a tolerance
        call check_tensor_values(dq_dSigma, exp_dq_dSigma, tol, passed)
        
        ! Check that the check passed
        call check(error, passed, .True., more = "Indiv. dq/dSigma test")
        if(allocated(error)) return

    end subroutine test_dq_to_dSigma

    subroutine test_dJ2_to_dSigma(error)
        type(error_type), allocatable, intent(out) ::  error
    
        ! Local variables
        real(kind = dp) :: exp_dJ2_dSigma(6), dJ2_dSigma(6)
        real(kind = dp) :: stress(6), dev_stress(6), mean_stress
        real(kind = dp), parameter :: tol = 1e-9
        logical :: passed = .False.

        call random_number(stress)

        mean_stress = calc_mean_stress(stress)

        dev_stress = calc_dev_stess(stress, mean_stress)
        
        ! Calc the value using the function
        dJ2_dSigma = calc_dJ2_to_dSigma(dev_stress)

        ! Calc the value using another method
        exp_dJ2_dSigma(1) = dev_stress(1)
        exp_dJ2_dSigma(2) = dev_stress(2)
        exp_dJ2_dSigma(3) = dev_stress(3) 
        exp_dJ2_dSigma(4) = 2.0_dp * dev_stress(4)
        exp_dJ2_dSigma(5) = 2.0_dp * dev_stress(5)
        exp_dJ2_dSigma(6) = 2.0_dp * dev_stress(6)

        ! Check that all the values are the same within a tolerance
        call check_tensor_values(dJ2_dSigma, exp_dJ2_dSigma, tol, passed)
        
        ! Check that the check passed
        call check(error, passed, .True., more = "Indiv. dJ2/dSigma test")
        if(allocated(error)) return
        
    end subroutine test_dJ2_to_dSigma

    subroutine test_dJ3_to_dSigma(error)
        type(error_type), allocatable, intent(out) ::  error
    
        ! Local variables
        real(kind = dp) :: exp_dJ3_dSigma(6), dJ3_dSigma(6)
        real(kind = dp) :: stress(6), dev(6), mean_stress
        real(kind = dp), parameter :: tol = 1e-9
        logical :: passed = .False.
        
        call random_number(stress)

        mean_stress = calc_mean_stress(stress)

        dev = calc_dev_stess(stress, mean_stress)

        ! Calc the value using the function
        dJ3_dSigma = calc_dJ3_to_dSigma(dev)

        ! Calc the value using another method
        exp_dJ3_dSigma = calc_inc_driver_dJ3_to_dSigma(stress)

        ! Check that all the values are the same within a tolerance
        call check_tensor_values(dJ3_dSigma, exp_dJ3_dSigma , tol, passed)
        
        ! Check that the check passed
        call check(error, passed, .True., more = "Indiv. dJ3/dSigma test")
        if(allocated(error)) return
        
    
    end subroutine test_dJ3_to_dSigma

    subroutine test_dtheta_to_dSigma(error)
        type(error_type), allocatable, intent(out) ::  error
    
        ! Local variables
        real(kind = dp) :: exp_dtheta_dSigma(6), dtheta_dSigma(6)
        real(kind = dp) :: stress(6), dev(6), mean_stress
        real(kind = dp) :: J2, J3, theta, dJ3_dSigma(6)
        real(kind = dp), parameter :: tol = 1e-9
        logical :: passed = .False.
                
        call random_number(stress)

        mean_stress = calc_mean_stress(stress)

        dev = calc_dev_stess(stress, mean_stress)
        
        J2 = calc_J2_invariant(dev)
        J3 = calc_inc_driver_J3_invariant(dev)
        theta = calc_theta_s(J2, J3)
        dJ3_dSigma = calc_dJ3_to_dSigma(dev)

        ! Calc the value using the function
        dtheta_dSigma = calc_dtheta_to_dSigma(dJ3_dSigma, dev, J3, J2, theta)
        
        ! Calc the value using another method
        exp_dtheta_dSigma = calc_dtheta_to_dSigma_2(dJ3_dSigma, dev, J3, J2)

        ! Check that all the values are the same within a tolerance
        call check_tensor_values(dtheta_dSigma, exp_dtheta_dSigma, tol, passed)
        
        ! Check that the check passed
        call check(error, passed, .True., more = "Indiv. dtheta/dSigma test")
        if(allocated(error)) return
        
    end subroutine test_dtheta_to_dSigma

end module mod_test_stress_invar_deriv_suite