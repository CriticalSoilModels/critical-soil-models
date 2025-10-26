module mod_test_strain_invariants_suite
    ! Local imports
    use kind_precision_module, only : dp, i32
    use mod_strain_invariants, only : Get_strain_invariants
    
    ! Testdrive imports
    use testdrive, only : new_unittest, unittest_type, error_type, check
 
    implicit none
 
    private
    public :: collect_strain_invariants_suite
    ! Note the convention used in incremental driver is
 contains
 
 
    subroutine collect_strain_invariants_suite(testsuite)
       ! Collection of tests
       ! Inidividual tests are stored in unitest_type
 
       type(unittest_type), allocatable, intent(out) :: testsuite(:)
 
       ! Make the test suite using the subroutines defined below
       testsuite = [ &
          new_unittest("strain_invariant_check", test_strain_invariant)&

          ]
 
    end subroutine collect_strain_invariants_suite
    
    subroutine test_strain_invariant(error)
       type(error_type), allocatable, intent(out) :: error
 
       ! Local varaibles
       real(kind = dp) :: strain(6), eps_v, eps_q
       real(kind = dp), parameter :: exp_eps_v = 9.0_dp               , &
          exp_eps_q = 10.878112581387148_dp, &
          tol       = epsilon(eps_v)
 
       strain = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp, 11.0_dp, 13.0_dp]
 
       ! Calc the strain invariant for the predifed strain (in voigt notation)
       call Get_strain_invariants(strain, eps_v, eps_q)
 
       ! Check if the strain invariants meet the expected values
       call check(error, eps_v, exp_eps_v, more = "Eps Volumetric test")
       if(allocated(error)) return
 
       call check(error, eps_q, exp_eps_q, more = "Eps q test")
       if(allocated(error)) return
 
    end subroutine test_strain_invariant
 
 
 end module mod_test_strain_invariants_suite
 