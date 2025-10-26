module mod_array_helper
   use stdlib_kinds, only: dp, i32 =>int32

   implicit none
   private
   public :: reorder_real_array, all_unique_integer

contains

   ! Check that all the elements in an array are unique
    function all_unique_integer(arr) result(is_unique)
        integer, intent(in) :: arr(:)
        logical :: is_unique
        integer :: i, j

        is_unique = .true.

        do i = 1, size(arr) - 1
            do j = i + 1, size(arr)
                if (arr(i) == arr(j)) then
                    is_unique = .false.
                    return
                end if
            end do
        end do
    end function all_unique_integer

   function reorder_real_array(input_arr, new_order) result(output_arr)
      real(kind = dp)    , intent(in)  :: input_arr(:)
      integer(kind = i32), intent(in)  :: new_order(:)
      real(kind = dp)    , allocatable :: output_arr(:)

      ! Local variables
      integer(kind=i32) :: i, n
      logical :: is_valid_order

      n = size(input_arr)

      ! Check that the length of the input_arr and new_order is the same
      if (size(input_arr) /= size(new_order)) then
         print *, "Error: The length of the input array and the new order isn't the same"
         print *, "The length of the input arr is: ", size(input_arr)
         print *, "The length of the new order arr is: ", size(new_order)
         return  ! output_arr remains unallocated
      end if

      ! Check if new_order contains all unique values between 1 and n
      is_valid_order = all(new_order >= 1 .and. new_order <= n) .and. all_unique_integer(new_order)
      if (.not. is_valid_order) then
         print *, "Error: Invalid reordering array"
         return  ! output_arr remains unallocated
      end if

      ! Allocate the new array
      allocate(output_arr(n))

      ! Reorder the array
      do i = 1, n         
         output_arr(i) = input_arr( new_order(i))
      end do

   end function reorder_real_array



end module mod_array_helper
