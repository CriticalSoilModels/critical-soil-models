! POD types for the NorSand constitutive model.
! No logic here — just data layout.
module mod_norsand_types
   use mod_csm_kinds, only: wp
   implicit none
   private

   public :: norsand_params_t, norsand_state_t

   type :: norsand_params_t
      real(wp) :: G_0        !! Reference shear modulus [kPa]
      real(wp) :: p_ref      !! Reference mean stress [kPa] (positive)
      real(wp) :: nG         !! Shear modulus exponent [-]
      real(wp) :: nu         !! Poisson's ratio [-]
      real(wp) :: e_o        !! Initial void ratio [-]
      real(wp) :: Gamma      !! CSL altitude at p=1 kPa [-]
      real(wp) :: lambda_c   !! CSL slope (ln scale) [-]
      real(wp) :: R          !! OCR (used to initialise p_i) [-]
      real(wp) :: M_tc       !! Critical friction ratio (triaxial compression) [-]
      real(wp) :: N          !! Nova volumetric coupling coefficient [-]
      !! TODO: verify sign convention and role of N in yield surface shape (N=0 recovers standard NorSand)
      real(wp) :: chi_tc     !! Dilatancy coefficient (triaxial compression) [-]
      real(wp) :: H_0        !! Hardening modulus intercept [-]
      real(wp) :: H_y        !! Hardening modulus slope [-]
   end type norsand_params_t

   type :: norsand_state_t
      real(wp) :: G              !! Current shear modulus [kPa] (from pre_step)
      real(wp) :: K              !! Current bulk modulus [kPa] (from pre_step)
      real(wp) :: p              !! Current mean effective stress [kPa] (from pre_step, negative compression)
      real(wp) :: e              !! Void ratio [-] (evolves via volumetric plastic strain)
      real(wp) :: psi            !! State parameter ψ = e − e_c [-]
      real(wp) :: chi_tce        !! Current dilatancy coefficient [-]
      real(wp) :: p_i            !! Image mean stress [kPa] (main hardening variable)
      real(wp) :: M_i            !! Image stress ratio [-]
      logical  :: switch_yield   !! Yielding flag
      real(wp) :: eps_p(6)       !! Accumulated plastic strain (Voigt) [-]
   end type norsand_state_t

end module mod_norsand_types
