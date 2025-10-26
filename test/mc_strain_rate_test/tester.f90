!> Driver for unit testing
program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type, &
      & select_suite, run_selected, get_argument, test_failed, init_color_output
  
    use mod_test_stress_invariants_suite, only : collect_stress_invariants_suite
    use mod_test_stress_invar_deriv_suite, only: collect_stress_invar_deriv_suite
    use mod_test_strain_invariants_suite, only: collect_strain_invariants_suite
    use mod_test_strain_invar_deriv_suite, only: collect_strain_invar_deriv_suite
    ! use mod_test_yield_function, only : collect_yield_function_suite
    use mod_test_plastic_potential_suite, only: collect_plastic_potential_suite
    implicit none
    integer :: stat, is
    
    character(len=:), allocatable :: suite_name, test_name
    ! Testsuites are stored in testsuite_type
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'
  
    ! init varialbe to track if any of the tests have failed
    stat = 0
  
    ! Make a list of the test suites from the differnt modules
    testsuites = [ &
      new_testsuite("test_stress_invar_suite", collect_stress_invariants_suite), &
      new_testsuite("test_stress_invar_deriv_suite", collect_stress_invar_deriv_suite), &
      new_testsuite("test_strain_invar_suite", collect_strain_invariants_suite),&
      new_testsuite("test_strain_invar_deriv_suite", collect_strain_invar_deriv_suite), &
      ! new_testsuite("test_yield_function", collect_yield_function_suite), &
      new_testsuite("test_plastic_potential", collect_plastic_potential_suite) &
      ]
    ! Make the output colorful
    call init_color_output(.True.)
  
    ! Allows input of a single suite and test name for running individual tests
    call get_argument(1, suite_name)
    call get_argument(2, test_name)
  
    ! If a suite name is entered only run that suite of tests
    if (allocated(suite_name)) then
  
      ! Get the index for the test suite from the testsuites list
      is = select_suite(testsuites, suite_name)
      
      ! If suite is inside of the testsuites...
      if (is > 0 .and. is <= size(testsuites)) then
        ! Check if a test name inside of the entered suite was entered
        if (allocated(test_name)) then
          ! Write the Suite and the test name
          write(error_unit, fmt) "Suite:", testsuites(is)%name
          ! Run the selected test
          call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
          
          ! If got an error stop
          if (stat < 0) then
            error stop 1
          end if
        
        else
          ! If a test name isn't passed run all the tests in the suite
          write(error_unit, fmt) "Testing:", testsuites(is)%name ! Writes the name of the test
          
          ! Runs all the tests in a suite
          call run_testsuite(testsuites(is)%collect, error_unit, stat)
        end if
      
      ! If invalid testsuite name passed...
      else
        ! Print the available testsuites
        write(error_unit, fmt) "Available testsuites"
        ! Loop over all testsuites stored in testsuites list
        do is = 1, size(testsuites)
          write(error_unit, fmt) "-", testsuites(is)%name
        end do
        ! Throw an error because an invalid testsuite was entered
        error stop 1
      end if
  
    ! If no testsuites (and test name) entered then run all of the tests
    else
      ! Otherwise run all of the tests
      do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
      end do
    end if
  
    ! Check if any tests failed...
    if (stat > 0) then
      ! If any failed print the number that failed
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      
      ! Throw an error stopping the program
      error stop 1
    end if
  
  
  end program tester