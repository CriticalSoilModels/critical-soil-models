module mod_tensor_value_checker
   use kind_precision_module, only: dp, sp
   implicit none

   private

   public :: check_tensor_values

   interface check_tensor_values
      module procedure check_real_sp_scalar_value, check_real_dp_scalar_value
      module procedure check_real_sp_vector_values
      module procedure check_real_dp_vector_values
      module procedure check_real_sp_matrix_values
      module procedure check_real_dp_matrix_values
   end interface check_tensor_values

contains

   ! Subroutine for checking values of a real double-precision scalar
   subroutine check_real_sp_scalar_value(actual, expected, tol, passed)
      real(kind= sp), intent(in) :: actual, expected, tol
      logical, intent(out) :: passed
      real(kind = sp), allocatable :: diff

      ! Calculate the difference between the values
      diff = abs(actual - expected)

      ! Check that the difference between all the elements is less than the tolerance
      if (diff < tol) then
         passed = .True.
      else
         print *, "Error: Real(sp) scalar values do not match."
         print *, "Actual: ", actual
         print *, "Expected: ", expected
         print *, "Difference: ", diff
         passed = .False.
      end if
   end subroutine check_real_sp_scalar_value

   ! Subroutine for checking values of a real double-precision scalar
   subroutine check_real_dp_scalar_value(actual, expected, tol, passed)
      real(kind= dp), intent(in) :: actual, expected, tol
      logical, intent(out) :: passed
      real(kind = dp), allocatable :: diff

      ! Calculate the difference between the values
      diff = abs(actual - expected)

      ! Check that the difference between all the elements is less than the tolerance
      if (diff < tol) then
         passed = .True.
      else
         print *, "Error: Real(dp) scalar values do not match."
         print *, "Actual: ", actual
         print *, "Expected: ", expected
         print *, "Difference: ", diff
         passed = .False.
      end if
   end subroutine check_real_dp_scalar_value

   ! Subroutine for checking values of a real single-precision vector
   subroutine check_real_sp_vector_values(actual, expected, tol, passed)
      real(kind= sp), intent(in) :: actual(:), expected(:), tol
      logical, intent(out) :: passed
      real(kind = sp), allocatable :: diff(:)

      ! Calculate the difference between the values
      diff = abs(actual - expected)

      ! Check that the difference between all the elements is less than the tolerance
      if (all(diff < tol)) then
         passed = .True.
      else
         print *, "Error: Real(sp) vector values do not match."
         print *, "Actual: ", actual
         print *, "Expected: ", expected
         print *, "Difference: ", diff
         passed = .False.
      end if
   end subroutine check_real_sp_vector_values

   ! Subroutine for checking values of a real single-precision vector
   subroutine check_real_dp_vector_values(actual, expected, tol, passed)
    real(kind= dp), intent(in) :: actual(:), expected(:), tol
    logical, intent(out) :: passed
    real(kind = dp), allocatable :: diff(:)

    ! Calculate the difference between the values
    diff = abs(actual - expected)

    ! Check that the difference between all the elements is less than the tolerance
    if (all(diff < tol)) then
       passed = .True.
    else
       print *, "Error: Real(dp) vector values do not match."
       print *, "Actual: ", actual
       print *, "Expected: ", expected
       print *, "Difference: ", diff
       passed = .False.
    end if
 end subroutine check_real_dp_vector_values

   ! Subroutine for checking values of a real single-precision matrix
   subroutine check_real_sp_matrix_values(actual, expected, tol, passed)
      real(kind= sp), intent(in) :: actual(:,:), expected(:,:), tol
      logical, intent(out) :: passed
      real(kind = sp), allocatable :: diff(:,:)

      ! Calculate the difference between the values
      diff = abs(actual - expected)

      ! Check that the difference between all the elements is less than the tolerance
      if (all(diff < tol)) then
         passed = .True.
      else
         print *, "Error: Real(sp) matrix values do not match."
         print *, "Actual: ", actual
         print *, "Expected: ", expected
         print *, "Difference: ", diff
         passed = .False.
      end if
   end subroutine check_real_sp_matrix_values

   ! Subroutine for checking values of a real double-precision matrix
   subroutine check_real_dp_matrix_values(actual, expected, tol, passed)
      real(kind= dp), intent(in) :: actual(:,:), expected(:,:), tol
      logical, intent(out) :: passed
      real(kind = dp), allocatable :: diff(:,:)

      ! Calculate the difference between the values
      diff = abs(actual - expected)

      ! Check that the difference between all the elements is less than the tolerance
      if (all(diff < tol)) then
         passed = .True.
      else
         print *, "Error: Real(dp) matrix values do not match."
         print *, "Actual: ", actual
         print *, "Expected: ", expected
         print *, "Difference: ", diff
         passed = .False.
      end if
   end subroutine check_real_dp_matrix_values
end module mod_tensor_value_checker

