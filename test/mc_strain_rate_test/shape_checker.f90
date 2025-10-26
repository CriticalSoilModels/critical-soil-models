module mod_shape_checker
    use kind_precision_module, only: dp
    implicit none
    private
    public :: check_matrix_shape

    interface check_matrix_shape
        module procedure check_real_matrix_shape
        module procedure check_integer_matrix_shape
        module procedure check_complex_matrix_shape
    end interface check_matrix_shape

contains

    ! Subroutine for real matrices
    subroutine check_real_matrix_shape(actual, expected, passed)
        real(kind=dp), intent(in) :: actual(:,:), expected(:,:)
        logical, intent(out) :: passed

        passed = (size(actual, 1) == size(expected, 1)) .and. &
                 (size(actual, 2) == size(expected, 2))

        if (.not. passed) then
            print *, "Error: Real matrices have different shapes."
            print *, "Actual shape: ", size(actual, 1), "x", size(actual, 2)
            print *, "Expected shape: ", size(expected, 1), "x", size(expected, 2)
        end if
    end subroutine check_real_matrix_shape

    ! Subroutine for integer matrices
    subroutine check_integer_matrix_shape(actual, expected, passed)
        integer, intent(in) :: actual(:,:), expected(:,:)
        logical, intent(out) :: passed

        passed = (size(actual, 1) == size(expected, 1)) .and. &
                 (size(actual, 2) == size(expected, 2))

        if (.not. passed) then
            print *, "Error: Integer matrices have different shapes."
            print *, "Actual shape: ", size(actual, 1), "x", size(actual, 2)
            print *, "Expected shape: ", size(expected, 1), "x", size(expected, 2)
        end if
    end subroutine check_integer_matrix_shape

    ! Subroutine for complex matrices
    subroutine check_complex_matrix_shape(actual, expected, passed)
        complex(kind=dp), intent(in) :: actual(:,:), expected(:,:)
        logical, intent(out) :: passed

        passed = (size(actual, 1) == size(expected, 1)) .and. &
                 (size(actual, 2) == size(expected, 2))

        if (.not. passed) then
            print *, "Error: Complex matrices have different shapes."
            print *, "Actual shape: ", size(actual, 1), "x", size(actual, 2)
            print *, "Expected shape: ", size(expected, 1), "x", size(expected, 2)
        end if
    end subroutine check_complex_matrix_shape

end module mod_shape_checker
