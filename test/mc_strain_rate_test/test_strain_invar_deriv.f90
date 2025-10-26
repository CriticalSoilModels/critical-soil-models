module mod_test_strain_invar_deriv_suite
    ! Local imports
    use kind_precision_module, only : real_type => dp, i32
    use mod_strain_invar_deriv, only : Get_dEpsq_to_dEps, calc_inc_driver_dEpsq_to_dEps
    use mod_strain_invariants, only: calc_eps_q_invariant, calc_dev_strain, Get_strain_invariants
    use mod_check_NaN_and_tensor_value, only : check_NaN_and_tensor_value

    ! Testdrive imports
    use testdrive, only : new_unittest, unittest_type, error_type, check
 
    implicit none
 
    private
    public :: collect_strain_invar_deriv_suite
    ! Note the convention used in incremental driver is
 contains
 
 
    subroutine collect_strain_invar_deriv_suite(testsuite)
       ! Collection of tests
       ! Inidividual tests are stored in unitest_type
 
       type(unittest_type), allocatable, intent(out) :: testsuite(:)
 
       ! Make the test suite using the subroutines defined below
       testsuite = [ &
          new_unittest("dEpsq/dEps", test_dEpsq_dEps)&

          ]
    end subroutine collect_strain_invar_deriv_suite
    
    subroutine test_dEpsq_dEps(error)
        type(error_type), allocatable, intent(out) :: error

        ! Local variables
        real(kind = real_type) :: exp_dEq_dEpsq(6), dEq_dEpsq(6)
        real(kind = real_type) :: Eps(6), Eps_v, Eps_q, dev_strain(6), test_arr(6)
        logical :: passed
        real(kind = real_type) :: a
        real(kind = real_type), parameter :: tol = 1e-9_real_type
        Eps = [1.0, 3.0, 5.0, 7.0, 11.0, 13.0]
        
        call Get_strain_invariants(Eps, Eps_v, Eps_q)
    
        ! Calc the value
        call Get_dEpsq_to_dEps(Eps_q, Eps, dEq_dEpsq)
        
        exp_dEq_dEpsq = calc_inc_driver_dEpsq_to_dEps(Eps)

        ! Check the values and make sure there isn't a NaN
        call check_NaN_and_tensor_value(dEq_dEpsq, exp_dEq_dEpsq, tol, passed)

        call check(error, passed, .True., more = "dEpsq_dEps")
    end subroutine test_dEpsq_dEps

 
 
 end module mod_test_strain_invar_deriv_suite
 