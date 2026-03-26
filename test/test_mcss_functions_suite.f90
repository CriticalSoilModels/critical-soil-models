module mod_test_mcss_functions_suite
   use mod_csm_kinds,      only: wp
   use mod_mcss_types,     only: mcss_params_t, mcss_state_t, DEFAULT_AS_PARAMS
   use mod_mcss_functions, only: mcss_yield_fn, mcss_flow_rule, &
                                  mcss_plastic_potential, mcss_update_hardening
   use ieee_arithmetic,    only: ieee_is_nan
   use testdrive,          only: new_unittest, unittest_type, error_type, check
   implicit none
   private
   public :: collect_mcss_functions_suite

   real(wp), parameter :: PI         = 4.0_wp * atan(1.0_wp)
   real(wp), parameter :: DEG_TO_RAD = PI / 180.0_wp

   ! Shared material parameters
   real(wp), parameter :: G_TEST   = 10000.0_wp               !! [kPa]
   real(wp), parameter :: NU_TEST  = 0.3_wp
   real(wp), parameter :: C_PEAK   = 10.0_wp                  !! [kPa]
   real(wp), parameter :: C_RES    = 2.0_wp                   !! [kPa]
   real(wp), parameter :: PHI_PEAK = 30.0_wp * DEG_TO_RAD     !! [rad]
   real(wp), parameter :: PHI_RES  = 20.0_wp * DEG_TO_RAD     !! [rad]
   real(wp), parameter :: PSI_PEAK = 10.0_wp * DEG_TO_RAD     !! [rad]
   real(wp), parameter :: PSI_RES  = 0.0_wp                   !! [rad]
   real(wp), parameter :: FACTOR   = 100.0_wp

contains

   subroutine collect_mcss_functions_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("yield_fn: hydrostatic (elastic)", test_yield_fn_hydrostatic),  &
         new_unittest("yield_fn: on surface at lode=0",  test_yield_fn_on_surface),   &
         new_unittest("yield_fn: plastic state",          test_yield_fn_plastic),      &
         new_unittest("hardening: zero plastic strain",   test_hardening_peak),        &
         new_unittest("hardening: exponential softening", test_hardening_softening),   &
         new_unittest("flow_rule: no NaN",                test_flow_rule_no_nan),      &
         new_unittest("plastic_potential: no NaN",        test_plastic_potential_no_nan) &
      ]
   end subroutine collect_mcss_functions_suite

   ! ---------------------------------------------------------------------------
   ! Helpers
   ! ---------------------------------------------------------------------------

   pure function make_params() result(p)
      type(mcss_params_t) :: p
      p%G         = G_TEST
      p%nu        = NU_TEST
      p%c_peak    = C_PEAK
      p%c_res     = C_RES
      p%phi_peak  = PHI_PEAK
      p%phi_res   = PHI_RES
      p%psi_peak  = PSI_PEAK
      p%psi_res   = PSI_RES
      p%factor    = FACTOR
      p%as_params = DEFAULT_AS_PARAMS
   end function make_params

   pure function make_peak_state() result(s)
      !! State at peak strength with no accumulated plastic strain.
      type(mcss_state_t) :: s
      s%c     = C_PEAK
      s%phi   = PHI_PEAK
      s%psi   = PSI_PEAK
      s%eps_p = 0.0_wp
   end function make_peak_state

   ! ---------------------------------------------------------------------------
   ! Yield function tests
   ! ---------------------------------------------------------------------------

   subroutine test_yield_fn_hydrostatic(error)
      !! Hydrostatic compression sits inside the yield surface: F < 0.
      !! p = -100 kPa, J = 0  =>  F = p*sin(phi) + a*sin(phi) - c*cos(phi) ≈ -58.7 kPa
      type(error_type), allocatable, intent(out) :: error

      type(mcss_params_t) :: params
      type(mcss_state_t)  :: state
      real(wp) :: F

      params = make_params()
      state  = make_peak_state()

      F = mcss_yield_fn(params, state, [-100.0_wp, -100.0_wp, -100.0_wp, 0.0_wp, 0.0_wp, 0.0_wp])

      call check(error, F < 0.0_wp, .true., more="hydrostatic: expected F < 0")
   end subroutine test_yield_fn_hydrostatic

   subroutine test_yield_fn_on_surface(error)
      !! Construction: at lode = 0, sig = [p+J, p, p-J, 0,0,0] with
      !! J = c*cos(phi) - p*sin(phi) places the stress exactly on the yield surface.
      !! The sig2 = p component ensures Sy = 0 => J3 = 0 => lode = 0.
      !! Tip smoothing introduces O(a^2 / J) ≈ 1e-7 kPa error; tolerance is 1e-4 kPa.
      type(error_type), allocatable, intent(out) :: error

      type(mcss_params_t) :: params
      type(mcss_state_t)  :: state
      real(wp) :: sig(6), F
      real(wp) :: p_val, J_yield

      params  = make_params()
      state   = make_peak_state()
      p_val   = -200.0_wp
      J_yield = C_PEAK * cos(PHI_PEAK) - p_val * sin(PHI_PEAK)

      sig = [p_val + J_yield, p_val, p_val - J_yield, 0.0_wp, 0.0_wp, 0.0_wp]
      F   = mcss_yield_fn(params, state, sig)

      call check(error, abs(F) < 1.0e-4_wp, .true., more="on-surface: expected |F| < 1e-4 kPa")
   end subroutine test_yield_fn_on_surface

   subroutine test_yield_fn_plastic(error)
      !! Scaling the on-surface deviatoric stress by 1.5x pushes outside: F > 0.
      type(error_type), allocatable, intent(out) :: error

      type(mcss_params_t) :: params
      type(mcss_state_t)  :: state
      real(wp) :: sig(6), F
      real(wp) :: p_val, J_yield

      params  = make_params()
      state   = make_peak_state()
      p_val   = -200.0_wp
      J_yield = C_PEAK * cos(PHI_PEAK) - p_val * sin(PHI_PEAK)

      sig = [p_val + 1.5_wp*J_yield, p_val, p_val - 1.5_wp*J_yield, 0.0_wp, 0.0_wp, 0.0_wp]
      F   = mcss_yield_fn(params, state, sig)

      call check(error, F > 0.0_wp, .true., more="plastic: expected F > 0")
   end subroutine test_yield_fn_plastic

   ! ---------------------------------------------------------------------------
   ! Hardening / softening tests
   ! ---------------------------------------------------------------------------

   subroutine test_hardening_peak(error)
      !! Zero plastic strain increment from zero accumulated strain => state at peak.
      type(error_type), allocatable, intent(out) :: error

      type(mcss_params_t) :: params
      type(mcss_state_t)  :: state
      real(wp), parameter :: TOL = 1.0e-12_wp

      params = make_params()
      state  = make_peak_state()

      call mcss_update_hardening(params, state, [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp])

      call check(error, abs(state%c   - C_PEAK)   < TOL, .true., more="peak: c /= c_peak")
      if (allocated(error)) return
      call check(error, abs(state%phi - PHI_PEAK) < TOL, .true., more="peak: phi /= phi_peak")
      if (allocated(error)) return
      call check(error, abs(state%psi - PSI_PEAK) < TOL, .true., more="peak: psi /= psi_peak")
   end subroutine test_hardening_peak

   subroutine test_hardening_softening(error)
      !! Pure deviatoric increment deps_p = [0.5, 0.5, -1, 0,0,0] gives eps_p_eq = 1.0.
      !! With factor = 1, expected: val = res + (peak - res)*exp(-1).
      !! Verification: eps_pm = 0, eps_p_eq = sqrt(2/3*(0.25+0.25+1.0)) = 1.0
      type(error_type), allocatable, intent(out) :: error

      type(mcss_params_t) :: params
      type(mcss_state_t)  :: state
      real(wp) :: c_exp, phi_exp, psi_exp
      real(wp), parameter :: FACTOR_UNIT = 1.0_wp
      real(wp), parameter :: EPS_P_EQ    = 1.0_wp
      real(wp), parameter :: TOL         = 1.0e-12_wp

      params        = make_params()
      params%factor = FACTOR_UNIT
      state         = make_peak_state()

      call mcss_update_hardening(params, state, &
         [0.5_wp, 0.5_wp, -1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp])

      c_exp   = C_RES   + (C_PEAK   - C_RES)   * exp(-EPS_P_EQ)
      phi_exp = PHI_RES + (PHI_PEAK - PHI_RES) * exp(-EPS_P_EQ)
      psi_exp = PSI_RES + (PSI_PEAK - PSI_RES) * exp(-EPS_P_EQ)

      call check(error, abs(state%c   - c_exp)   < TOL, .true., more="softening: c mismatch")
      if (allocated(error)) return
      call check(error, abs(state%phi - phi_exp) < TOL, .true., more="softening: phi mismatch")
      if (allocated(error)) return
      call check(error, abs(state%psi - psi_exp) < TOL, .true., more="softening: psi mismatch")
   end subroutine test_hardening_softening

   ! ---------------------------------------------------------------------------
   ! Gradient smoke tests
   ! ---------------------------------------------------------------------------

   subroutine test_flow_rule_no_nan(error)
      !! dF/dsig must not produce NaN at a general stress state or at J=0 (hydrostatic).
      type(error_type), allocatable, intent(out) :: error

      type(mcss_params_t) :: params
      type(mcss_state_t)  :: state
      real(wp) :: n(6)

      params = make_params()
      state  = make_peak_state()

      n = mcss_flow_rule(params, state, &
         [-100.0_wp, -200.0_wp, -300.0_wp, 50.0_wp, 30.0_wp, 20.0_wp])
      call check(error, all(.not. ieee_is_nan(n)), .true., &
         more="flow_rule: NaN at general stress state")
      if (allocated(error)) return

      n = mcss_flow_rule(params, state, &
         [-200.0_wp, -200.0_wp, -200.0_wp, 0.0_wp, 0.0_wp, 0.0_wp])
      call check(error, all(.not. ieee_is_nan(n)), .true., &
         more="flow_rule: NaN at hydrostatic (J=0) state")
   end subroutine test_flow_rule_no_nan

   subroutine test_plastic_potential_no_nan(error)
      !! dG/dsig must not produce NaN at a general stress state or at J=0 (hydrostatic).
      !! The J=0 case also exercises the psi-to-phi interpolation branch.
      type(error_type), allocatable, intent(out) :: error

      type(mcss_params_t) :: params
      type(mcss_state_t)  :: state
      real(wp) :: m(6)

      params = make_params()
      state  = make_peak_state()

      m = mcss_plastic_potential(params, state, &
         [-100.0_wp, -200.0_wp, -300.0_wp, 50.0_wp, 30.0_wp, 20.0_wp])
      call check(error, all(.not. ieee_is_nan(m)), .true., &
         more="plastic_potential: NaN at general stress state")
      if (allocated(error)) return

      m = mcss_plastic_potential(params, state, &
         [-200.0_wp, -200.0_wp, -200.0_wp, 0.0_wp, 0.0_wp, 0.0_wp])
      call check(error, all(.not. ieee_is_nan(m)), .true., &
         more="plastic_potential: NaN at hydrostatic (J=0) state")
   end subroutine test_plastic_potential_no_nan

end module mod_test_mcss_functions_suite
