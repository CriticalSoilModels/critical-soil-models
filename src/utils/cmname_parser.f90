! Utilities for parsing the Abaqus CMNAME string.
!
! Convention:  MATERIALNAME_INTEG_INTEGRATORNAME
!
! Examples:
!   "LINEAR_ELASTIC"              → material="LINEAR_ELASTIC", integrator="euler" (default)
!   "MCSS_INTEG_EULER"            → material="MCSS",           integrator="euler"
!   "MCSS_INTEG_CPA"              → material="MCSS",           integrator="cpa"
!
! Abaqus always uppercases CMNAME before passing it.
! material_name returns the uppercase prefix as-is.
! integrator_name lowercases the suffix to match the method string constants.

module mod_cmname_parser
   implicit none
   private
   public :: material_name, integrator_name

   character(len=*), parameter :: INTEG_SEP     = "_INTEG_"
   character(len=*), parameter :: DEFAULT_INTEG = "euler"

contains

   pure function material_name(cmname) result(name)
      !! Returns the material model prefix of CMNAME (before _INTEG_, or whole string).
      character(len=80), intent(in) :: cmname
      character(len=80)             :: name
      integer :: pos

      pos = index(cmname, INTEG_SEP)
      if (pos > 0) then
         name = cmname(1:pos-1)
      else
         name = trim(cmname)
      end if
   end function material_name

   pure function integrator_name(cmname) result(name)
      !! Returns the integrator suffix of CMNAME, lowercased.
      !! Returns "euler" (default) if no _INTEG_ separator is present.
      character(len=80), intent(in) :: cmname
      character(len=80)             :: name
      integer :: pos

      pos = index(cmname, INTEG_SEP)
      if (pos > 0) then
         name = to_lower(trim(cmname(pos + len(INTEG_SEP):)))
      else
         name = DEFAULT_INTEG
      end if
   end function integrator_name

   ! ---------------------------------------------------------------------------
   ! Private helper
   ! ---------------------------------------------------------------------------

   pure function to_lower(str) result(lower)
      character(len=*), intent(in) :: str
      character(len=len(str))      :: lower
      integer :: i, c

      do i = 1, len(str)
         c = iachar(str(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) then
            lower(i:i) = achar(c + 32)
         else
            lower(i:i) = str(i:i)
         end if
      end do
   end function to_lower

end module mod_cmname_parser
