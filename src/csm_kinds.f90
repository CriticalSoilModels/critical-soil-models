! Working precision kinds for the critical-soil-models library.
!
! To change precision across the entire library, update wp here.
! All source files import wp from this module rather than using stdlib_kinds directly.

module mod_csm_kinds
   use stdlib_kinds, only: dp, sp
   implicit none
   private

   public :: sp                              !! Single precision (re-exported from stdlib)
   integer, parameter, public :: wp = dp    !! Working precision (currently double)

end module mod_csm_kinds
