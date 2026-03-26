! POD types for the linear elastic model.
! No logic here — just data layout.
module mod_le_types
   use mod_csm_kinds, only: wp
   implicit none
   private

   public :: le_params_t, le_state_t

   type :: le_params_t
      real(wp) :: G   !! Shear modulus [kPa]
      real(wp) :: nu  !! Poisson's ratio [-]
   end type le_params_t

   type :: le_state_t
      !! Linear elastic has no internal state variables.
      !! Empty type preserves the uniform (params, state, ...) call signature.
   end type le_state_t

end module mod_le_types
