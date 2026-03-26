module mod_test_strain_invariants_suite
    ! Local imports
    use mod_csm_kinds, only: wp
    use mod_strain_invariants, only : calc_eps_inv
    
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
       real(kind=wp) :: strain(6), eps_v, eps_q
       real(kind=wp), parameter :: exp_eps_v = 9.0_wp               , &
          exp_eps_q = 10.878112581387148_wp, &
          tol       = epsilon(eps_v)
 
       strain = [1.0_wp, 3.0_wp, 5.0_wp, 7.0_wp, 11.0_wp, 13.0_wp]
 
       ! Calc the strain invariant for the predifed strain (in voigt notation)
       call calc_eps_inv(strain, eps_v, eps_q)
 
       ! Check if the strain invariants meet the expected values
       call check(error, eps_v, exp_eps_v, more = "Eps Volumetric test")
       if(allocated(error)) return
 
       call check(error, eps_q, exp_eps_q, more = "Eps q test")
       if(allocated(error)) return
 
    end subroutine test_strain_invariant
 
 
 end module mod_test_strain_invariants_suite
 