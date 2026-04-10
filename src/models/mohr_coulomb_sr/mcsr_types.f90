!! Types for the Mohr-Coulomb Strain Rate (MCSR) model.
!!
!! `mcsr_params_t` holds fixed parameters (unpacked from PROPS once per call).
!! `mcsr_state_t` holds the evolving internal state (round-tripped through STATEV).
!!
!! ### STATEV layout (14 entries)
!!
!! | Index | Symbol      | Description                          | Units |
!! |-------|-------------|--------------------------------------|-------|
!! | 1     | G           | Current shear modulus                | kPa   |
!! | 2     | K           | Current bulk modulus                 | kPa   |
!! | 3     | eta_y       | Friction ratio                       | —     |
!! | 4     | dilation    | Dilation                             | —     |
!! | 5     | I_coeff     | Inertial coefficient                 | —     |
!! | 6     | switch_yield| Yielding flag (0.0 / 1.0)            | —     |
!! | 7–12  | eps_p       | Accumulated plastic strain (Voigt)   | —     |
!! | 13    | N_i         | Strain rate smoothing counter        | —     |
!! | 14    | SUM_rate    | Strain rate sum for smoothing        | —     |
!!
!! `M_lode` is a transient field set by `mcsr_update_rate_state` at the start
!! of each timestep and is NOT persisted in STATEV.

module mod_mcsr_types
   use mod_csm_kinds, only: wp
   implicit none
   private

   public :: mcsr_params_t, mcsr_state_t

   type :: mcsr_params_t
      !! Fixed model parameters (from PROPS, does not evolve).
      real(wp) :: G_0          !! Initial shear modulus [kPa]
      real(wp) :: nu           !! Poisson's ratio [-]
      real(wp) :: M_tc         !! Critical stress ratio, triaxial compression [-]
      real(wp) :: N            !! Nova's volumetric coupling coefficient [-]
                               !! TODO: verify role — appears in eta_y = M - dil*(1-N)
      real(wp) :: D_min        !! Minimum dilation [-]
      real(wp) :: h            !! Dilation hardening parameter [-]
      real(wp) :: alpha_G      !! Shear modulus rate exponent [-]
      real(wp) :: alpha_K      !! Bulk modulus rate exponent [-]
      real(wp) :: alpha_D      !! Dilation rate exponent [-]
      real(wp) :: D_part       !! Particle diameter [mm]
      real(wp) :: G_s          !! Specific gravity [-]
      real(wp) :: ref_e_rate   !! Reference strain rate [1/s]
      logical  :: switch_smooth   !! Activate strain rate smoothing
      integer  :: N_S             !! Degree of smoothing [-]
      real(wp) :: yield_tol    !! Yield surface tolerance [-]
      integer  :: max_iters    !! Maximum substep iterations
   end type mcsr_params_t

   type :: mcsr_state_t
      !! Evolving internal state — updated each timestep.
      real(wp) :: G            !! Current shear modulus [kPa]
      real(wp) :: K            !! Current bulk modulus [kPa]
      real(wp) :: eta_y        !! Current friction ratio [-]
      real(wp) :: dilation     !! Current dilation [-]
      real(wp) :: I_coeff      !! Inertial coefficient (from previous step) [-]
      real(wp) :: M_lode       !! Lode-corrected M — computed pre-step, NOT in STATEV [-]
      logical  :: switch_yield !! Yielding flag
      real(wp) :: eps_p(6)     !! Accumulated plastic strain (Voigt) [-]
      integer  :: N_i          !! Strain rate smoothing counter
      real(wp) :: SUM_rate     !! Accumulated strain rate sum for smoothing
   end type mcsr_state_t

end module mod_mcsr_types
