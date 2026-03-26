! POD types for the Mohr-Coulomb Strain Softening model.
! No logic here — just data layout.
module mod_mcss_types
   use mod_csm_kinds, only: wp
   implicit none
   private

   public :: mcss_params_t, mcss_state_t

   type :: mcss_params_t
      real(wp) :: G           !! Shear modulus [kPa]
      real(wp) :: nu          !! Poisson's ratio [-]
      real(wp) :: c_peak      !! Peak cohesion [kPa]
      real(wp) :: c_res       !! Residual cohesion [kPa]
      real(wp) :: phi_peak    !! Peak friction angle [rad]
      real(wp) :: phi_res     !! Residual friction angle [rad]
      real(wp) :: psi_peak    !! Peak dilation angle [rad]
      real(wp) :: psi_res     !! Residual dilation angle [rad]
      real(wp) :: factor      !! Softening rate (shape parameter) [-]
   end type mcss_params_t

   type :: mcss_state_t
      real(wp) :: c           !! Current cohesion [kPa]
      real(wp) :: phi         !! Current friction angle [rad]
      real(wp) :: psi         !! Current dilation angle [rad]
      real(wp) :: eps_p(6)    !! Accumulated plastic strain (Voigt) [-]
   end type mcss_state_t

end module mod_mcss_types
