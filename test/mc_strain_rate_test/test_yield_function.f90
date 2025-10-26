module mod_test_yield_function
    ! Local imports
    use kind_precision_module, only: real_type => dp, i32
    use mod_yield_function, only : Get_dF_to_dSigma, Get_dF_to_dSigma_3
    use mod_shape_checker , only : check_matrix_shape
    use mod_tensor_value_checker, only: check_tensor_values
    use ieee_arithmetic, only: ieee_is_nan
    
    ! Testdrive imports
    use testdrive, only : new_unittest, unittest_type, error_type, check
 
    implicit none
 
    private
    public :: collect_yield_function_suite
    ! Note the convention used in incremental driver is
 contains
 
 
     subroutine collect_yield_function_suite(testsuite)
         ! Collection of tests
         ! Inidividual tests are stored in unitest_type
 
         type(unittest_type), allocatable, intent(out) :: testsuite(:)
 
         ! Make the test suite using the subroutines defined below
         testsuite = [ &
         new_unittest("dF_to_dSigma_check", test_dF_to_dSigma)  &
 
         ]
 
     end subroutine collect_yield_function_suite
    
    subroutine test_dF_to_dSigma(error)
        ! Check the derivatives of the yield function
        ! TODO: Need to do the hand calc for this function
        type(error_type), allocatable, intent(out) :: error
        
        ! Local variables
        real(kind = real_type) :: &
            M_tc  = 1.0  , &
            eta_y = 1.5  , &
            Sig(6)       , &
            n_vec(6), &
            tol = 1e-9
        real(kind = real_type) :: exp_n_vec(6)
        logical :: passed

        Sig = [1.0, 3.0, 5.0, 7.0, 11.0, 13.0]

        ! The expected n_vector
        
        ! Call the original df_dsigma
        call Get_dF_to_dSigma(M_tc, eta_y, Sig, exp_n_vec)
        
        ! Call the updated df_dsigma
        call Get_dF_to_dSigma_3(M_tc, eta_y, Sig, n_vec)
        
        ! Check if ther are any NaN values in the calculation 
        if (all(ieee_is_nan(n_vec) .eqv. .False.)) then
            ! There are no NaNs
            call check_tensor_values(n_vec, exp_n_vec, tol, passed)
        else
            ! Force an error
            passed = .False. ! Fail the test
            print *, 'Theres a NaN in the results'
            print *, n_vec
        end if

        ! Check for error (logical)
        call check(error, passed,.False., more = "Need to do the hand calc. Currently not valid comparison")
        if ( allocated(error) ) return 
        
    end subroutine test_dF_to_dSigma

 end module mod_test_yield_function
 