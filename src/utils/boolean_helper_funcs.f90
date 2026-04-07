!! Boolean / logical conversion utilities.
!!
!! Bridges legacy code that represents booleans as `real(wp)` (0.0 = false, 1.0 = true)
!! with standard Fortran `logical` types.

module mod_bool_helper
   use mod_csm_kinds, only: wp
   implicit none

contains

   subroutine dbltobool(A, B)
      !! Convert a double-precision flag (0.0 or 1.0) to a logical.
      real(wp), intent(in)  :: A  !! Input flag: 0.0 → false, anything ≥ 1.0 → true
      logical,  intent(out) :: B  !! Output logical
      if (A < 1.0) then
         B = .false.
      else
         B = .true.
      end if
   end subroutine dbltobool

   real(wp) function logic2dbl(a)
      !! Convert a logical to a double-precision flag (false → 0.0, true → 1.0).
      logical, intent(in) :: a  !! Input logical
      if (a) then
         logic2dbl = 1.d0
      else
         logic2dbl = 0.d0
      end if
   end function logic2dbl

end module mod_bool_helper
