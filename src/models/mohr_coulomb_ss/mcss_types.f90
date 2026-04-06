! POD types for the Mohr-Coulomb Strain Softening model.
! No logic here — just data layout.
module mod_mcss_types
   use mod_csm_kinds, only: wp
   implicit none
   private

   public :: mcss_params_t, mcss_state_t
   public :: abbo_sloan_params_t
   public :: DEFAULT_AS_PARAMS, DEFAULT_SMOOTH_COEFF

   ! ---------------------------------------------------------------------------
   ! Default Abbo & Sloan (1995) smoothing constants — LodeT = 29.5 degrees.
   ! Override by supplying PROPS(14:19); zero or missing entries use these values.
   ! ---------------------------------------------------------------------------
   real(wp), parameter :: DEFAULT_LODE_TR      = 0.514872129338327_wp  !! Transition Lode angle [rad]
   real(wp), parameter :: DEFAULT_A1           = 7.138654723242414_wp
   real(wp), parameter :: DEFAULT_A2           = 6.112267270920612_wp
   real(wp), parameter :: DEFAULT_B1           = 6.270447753139589_wp
   real(wp), parameter :: DEFAULT_B2           = 6.398760841429403_wp
   real(wp), parameter :: DEFAULT_SMOOTH_COEFF = 0.0005_wp             !! Tip smoothing: a = coeff*c*cot(phi)

   type :: abbo_sloan_params_t
      real(wp) :: lode_tr      !! Transition Lode angle [rad]
      real(wp) :: A1           !! Corner-rounding constant
      real(wp) :: A2           !! Corner-rounding constant
      real(wp) :: B1           !! Corner-rounding constant
      real(wp) :: B2           !! Corner-rounding constant
      real(wp) :: smooth_coeff !! Tip smoothing coefficient [-]
   end type abbo_sloan_params_t

   !> Pre-built default instance — use as a fallback or initialiser.
   type(abbo_sloan_params_t), parameter :: DEFAULT_AS_PARAMS = abbo_sloan_params_t( &
      lode_tr      = DEFAULT_LODE_TR,      &
      A1           = DEFAULT_A1,           &
      A2           = DEFAULT_A2,           &
      B1           = DEFAULT_B1,           &
      B2           = DEFAULT_B2,           &
      smooth_coeff = DEFAULT_SMOOTH_COEFF  &
   )

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
      type(abbo_sloan_params_t) :: as_params  !! Abbo & Sloan smoothing constants
   end type mcss_params_t

   type :: mcss_state_t
      real(wp) :: c           !! Current cohesion [kPa]
      real(wp) :: phi         !! Current friction angle [rad]
      real(wp) :: psi         !! Current dilation angle [rad]
      real(wp) :: eps_p(6)    !! Accumulated plastic strain (Voigt) [-]
   end type mcss_state_t

end module mod_mcss_types
