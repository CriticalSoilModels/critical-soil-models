module mod_check_NaN_and_tensor_value
   use kind_precision_module, only: dp, sp
   use ieee_arithmetic, only: ieee_is_nan
   use mod_tensor_value_checker, only: check_tensor_values
   implicit none

   private

   public :: check_NaN_and_tensor_value

   interface check_NaN_and_tensor_value
      module procedure check_NaN_and_tensor_value_sp_actual, check_NaN_and_tensor_value_dp_actual
      module procedure check_NaN_and_tensor_value_sp_1d, check_NaN_and_tensor_value_dp_1d
      module procedure check_NaN_and_tensor_value_sp_2d, check_NaN_and_tensor_value_dp_2d
   end interface check_NaN_and_tensor_value

contains

   ! For single precision (real(kind = sp)) 1D array
   subroutine check_NaN_and_tensor_value_sp_1d(actual, expected, tol, passed)
      real(kind = sp), intent(in) :: actual(:)
      real(kind = sp), intent(in) :: expected(:)
      real(kind = sp), intent(in) :: tol
      logical, intent(out) :: passed

      if (all(ieee_is_nan(actual) .eqv. .false.)) then
         call check_tensor_values(actual, expected, tol, passed)
      else
         passed = .false.
         print *, "There is a NaN in the result (double precision 1D array)"
         print *, actual
      end if
   end subroutine check_NaN_and_tensor_value_sp_1d

   ! For double precision (real(kind = dp)) 1D array
   subroutine check_NaN_and_tensor_value_dp_1d(actual, expected, tol, passed)
      real(kind = dp), intent(in) :: actual(:)
      real(kind = dp), intent(in) :: expected(:)
      real(kind = dp), intent(in) :: tol
      logical, intent(out) :: passed

      if (all(ieee_is_nan(actual) .eqv. .false.)) then
         call check_tensor_values(actual, expected, tol, passed)
      else
         passed = .false.
         print *, "There is a NaN in the result (double precision 1D array)"
         print *, actual
      end if
   end subroutine check_NaN_and_tensor_value_dp_1d

   ! For single precision (real(kind = sp)) 2D array
   subroutine check_NaN_and_tensor_value_sp_2d(actual, expected, tol, passed)
      real(kind = sp), intent(in) :: actual(:,:)
      real(kind = sp), intent(in) :: expected(:,:)
      real(kind = sp), intent(in) :: tol
      logical, intent(out) :: passed

      if (all(ieee_is_nan(actual) .eqv. .false.)) then
         call check_tensor_values(actual, expected, tol, passed)
      else
         passed = .false.
         print *, "There is a NaN in the result (single precision 2D array)"
         print *, actual
      end if
   end subroutine check_NaN_and_tensor_value_sp_2d

   ! For double precision (real(kind = dp)) 2D array
   subroutine check_NaN_and_tensor_value_dp_2d(actual, expected, tol, passed)
      real(kind = dp), intent(in) :: actual(:,:)
      real(kind = dp), intent(in) :: expected(:,:)
      real(kind = dp), intent(in) :: tol
      logical, intent(out) :: passed

      if (all(ieee_is_nan(actual) .eqv. .false.)) then
         call check_tensor_values(actual, expected, tol, passed)
      else
         passed = .false.
         print *, "There is a NaN in the result (double precision 2D array)"
         print *, actual
      end if
   end subroutine check_NaN_and_tensor_value_dp_2d

   ! For single precision (real(kind = sp)) actual
   subroutine check_NaN_and_tensor_value_sp_actual(actual, expected, tol, passed)
      real(kind = sp), intent(in) :: actual
      real(kind = sp), intent(in) :: expected
      real(kind = sp), intent(in) :: tol
      logical, intent(out) :: passed

      if (.not. ieee_is_nan(actual)) then
         call check_tensor_values(actual, expected, tol, passed)
      else
         passed = .false.
         print *, "There is a NaN in the result (single precision scalar)"
         print *, actual
      end if
   end subroutine check_NaN_and_tensor_value_sp_actual

   ! For double precision (real(kind = dp)) actual
   subroutine check_NaN_and_tensor_value_dp_actual(actual, expected, tol, passed)
      real(kind = dp), intent(in) :: actual
      real(kind = dp), intent(in) :: expected
      real(kind = dp), intent(in) :: tol
      logical, intent(out) :: passed

      if (.not. ieee_is_nan(actual)) then
         call check_tensor_values(actual, expected, tol, passed)
      else
         passed = .false.
         print *, "There is a NaN in the result (double precision scalar)"
         print *, actual
      end if
   end subroutine check_NaN_and_tensor_value_dp_actual

end module mod_check_NaN_and_tensor_value
