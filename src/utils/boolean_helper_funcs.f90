! Functions that help working with booleans (True or False values)

module mod_bool_helper
   use kind_precision_module, only: real_type => dp
   implicit none

contains
   subroutine dbltobool(A,B)
      !******************************************************************
      ! Takes a double which values are either 1.0 or 0.0 and returns a *
      ! Boolean
      !******************************************************************
      implicit none
      double precision, intent(in):: A
      logical, intent(out):: B
      if (A<1.0) then
         B=.false.
      else
         B=.true.
      endif
   end subroutine dbltobool

   real(real_type) function logic2dbl(a)
      logical, intent(in) :: a

      if (a) then
         logic2dbl = 1.d0
      else
         logic2dbl = 0.d0
      end if
   end function logic2dbl

end module mod_bool_helper
