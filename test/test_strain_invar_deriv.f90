module mod_test_strain_invar_deriv_suite
    ! Local imports
    use mod_csm_kinds, only: wp
    use mod_strain_invar_deriv, only : calc_deps_q_by_deps
    use mod_strain_invar_refs, only: calc_deps_q_by_deps_full
    use mod_strain_invariants, only: calc_eps_q_inv, calc_dev_strain, calc_eps_inv
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
        real(wp) :: exp_dEq_dEpsq(6), dEq_dEpsq(6)
        real(wp) :: Eps(6), Eps_v, Eps_q, dev_strain(6), test_arr(6)
        logical :: passed
        real(wp) :: a
        real(wp), parameter :: tol = 1e-9_wp
        Eps = [1.0, 3.0, 5.0, 7.0, 11.0, 13.0]
        
        call calc_eps_inv(Eps, Eps_v, Eps_q)
    
        ! Calc the value
        dEq_dEpsq = calc_deps_q_by_deps(Eps_q, Eps)
        
        exp_dEq_dEpsq = calc_deps_q_by_deps_full(Eps)

        ! Check the values and make sure there isn't a NaN
        call check_NaN_and_tensor_value(dEq_dEpsq, exp_dEq_dEpsq, tol, passed)

        call check(error, passed, .True., more = "dEpsq_dEps")
    end subroutine test_dEpsq_dEps

 
 
 end module mod_test_strain_invar_deriv_suite
 