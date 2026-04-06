module mod_abbo_sloan_presets
   !! Factory function for Abbo & Sloan (1995) corner-rounding constants.
   !!
   !! The Abbo-Sloan yield surface uses five constants (lode_tr, A1, A2, B1, B2)
   !! that depend on the chosen transition Lode angle. This module provides
   !! pre-calibrated sets for the two published angles (25° and 29.5°).
   !! The 30° case is excluded — the formula is singular there (cos(90°) = 0).
   use mod_csm_kinds,  only: wp
   use mod_mcss_types, only: abbo_sloan_params_t, DEFAULT_SMOOTH_COEFF
   implicit none
   private
   public :: abbo_sloan_preset

contains

   pure function abbo_sloan_preset(lode_tr_deg) result(asp)
      !! Return the Abbo & Sloan (1995) smoothing constants for a given
      !! transition Lode angle in degrees.
      !!
      !! Supported values: 25.0, 29.5.
      !! Any unrecognised angle falls back to the 29.5° set.
      !! `smooth_coeff` (tip smoothing) is always DEFAULT_SMOOTH_COEFF —
      !! it is independent of the corner-rounding preset.
      real(wp), intent(in) :: lode_tr_deg
      type(abbo_sloan_params_t) :: asp

      ! Match on tenths-of-degrees (integer) to avoid floating-point equality.
      select case (nint(lode_tr_deg * 10.0_wp))

      case (250)   ! 25.0°
         asp = abbo_sloan_params_t( &
            lode_tr      = 0.436332312998582_wp,  &
            A1           = 1.432052062044227_wp,   &
            A2           = 0.406941858374615_wp,   &
            B1           = 0.544290524902313_wp,   &
            B2           = 0.673903324498392_wp,   &
            smooth_coeff = DEFAULT_SMOOTH_COEFF    )

      case default ! 29.5° — and fallback for any unrecognised angle
         asp = abbo_sloan_params_t( &
            lode_tr      = 0.514872129338327_wp,  &
            A1           = 7.138654723242414_wp,   &
            A2           = 6.112267270920612_wp,   &
            B1           = 6.270447753139589_wp,   &
            B2           = 6.398760841429403_wp,   &
            smooth_coeff = DEFAULT_SMOOTH_COEFF    )

      end select
   end function abbo_sloan_preset

end module mod_abbo_sloan_presets
