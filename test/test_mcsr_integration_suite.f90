!! Integration tests for the MCSR model comparing euler_substep and cpa_step.
!!
!! Both integrators operate through the same csm_model_t interface; these tests
!! verify that:
!!   1. On a purely elastic step, both give identical stress (they must — no
!!      plastic flow occurs so the result is the elastic predictor).
!!   2. On a plastic step, both return stress on/inside the yield surface.
!!   3. For a small plastic overshoot both methods agree within 1%.
!!
!! ### Initial state
!!
!! Isotropic stress at p = -100 kPa (sig = [-100,-100,-100,0,0,0]).
!! At q = 0, calc_lode_inv returns θ = π/6, giving:
!!   M_lode = M_tc * (1 + 0.25 * cos(π/2)^1.2) = M_tc
!! With zero dilation and N = 0: eta_y = M_tc.
!! Yield radius: q_yield = M_tc * |p| = 1.2 * 100 = 120 kPa.
!! Elastic yield strain ≈ 120 / (sqrt(3) * 2*G) ≈ 3.46e-3.

module mod_test_mcsr_integration_suite
   use mod_csm_kinds,   only: wp
   use mod_mcsr_model,  only: mcsr_model_t, mcsr_from_props, mcsr_load_state
   use mod_euler_substep, only: euler_substep, integrator_params_t
   use mod_cpa,          only: cpa_step
   use testdrive,         only: new_unittest, unittest_type, error_type, check
   implicit none
   private
   public :: collect_mcsr_integration_suite

   real(wp), parameter :: PI = 4.0_wp * atan(1.0_wp)

   ! Model parameters
   real(wp), parameter :: G_0_TEST  = 10000.0_wp  !! [kPa]
   real(wp), parameter :: NU_TEST   = 0.3_wp
   real(wp), parameter :: M_TC_TEST = 1.2_wp
   real(wp), parameter :: YIELD_TOL = 1.0e-6_wp
   real(wp), parameter :: K_0_TEST  = &
      2.0_wp * G_0_TEST * (1.0_wp + NU_TEST) / (3.0_wp * (1.0_wp - 2.0_wp * NU_TEST))

   ! Initial stress: isotropic at p = -100 kPa
   real(wp), parameter :: SIG_INIT(6) = &
      [-100.0_wp, -100.0_wp, -100.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

   ! Note: integrator_params_t is NOT declared as a module-level parameter.
   ! gfortran miscompiles derived-type parameter initialisations at module scope
   ! when the type mixes real and integer fields, causing a heap corruption at
   ! program startup. Declare iparams as a local variable in each test instead.

contains

   subroutine collect_mcsr_integration_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("elastic step: euler and cpa identical",     test_elastic_agreement),    &
         new_unittest("plastic step: euler on yield surface",        test_euler_on_surface),     &
         new_unittest("plastic step: cpa on yield surface",         test_cpa_on_surface),      &
         new_unittest("small plastic step: euler and cpa agree",    test_stress_comparison)     &
      ]
   end subroutine collect_mcsr_integration_suite

   ! ---------------------------------------------------------------------------
   ! Helper: build MCSR model at isotropic initial state.
   !
   ! At q=0, lode angle = pi/6, so M_lode = M_tc (see module doc).
   ! With zero dilation and N=0: eta_y = M_tc.
   ! M_lode is set explicitly — it is transient (not in STATEV) and would
   ! otherwise be zero after mcsr_load_state, causing update_hardening to
   ! incorrectly zero eta_y on the first plastic increment.
   ! ---------------------------------------------------------------------------
   subroutine make_model(model)
      type(mcsr_model_t), intent(out) :: model
      real(wp) :: props(18), statev(14)

      props = [G_0_TEST, NU_TEST, M_TC_TEST, &  ! G_0, nu, M_tc
               0.0_wp,                        &  ! N (Nova coupling)
               0.3_wp, 2.0_wp,                &  ! D_min, h
               0.2_wp, 0.5_wp, 0.2_wp,        &  ! alpha_G, alpha_K, alpha_D
               1.0e-3_wp, 2.65_wp, 1.0e-3_wp, &  ! D_part, G_s, ref_e_rate
               0.0_wp, 5.0_wp,                &  ! switch_smooth, N_S
               1.0_wp,                        &  ! switch_original (Wang)
               YIELD_TOL, 100.0_wp, 0.0_wp]      ! yield_tol, max_iters, switch_integ

      statev       = 0.0_wp
      statev(1)    = G_0_TEST   ! G  [kPa]
      statev(2)    = K_0_TEST   ! K  [kPa]
      statev(3)    = M_TC_TEST  ! eta_y = M_tc at zero dilation, N=0, lode=pi/6
      statev(4)    = 0.0_wp     ! dilation
      statev(5)    = 1.0e-3_wp  ! I_coeff (at reference rate)

      model = mcsr_from_props(props)
      call mcsr_load_state(model, statev)

      ! Initialise transient field — not in STATEV but needed by update_hardening
      model%state%M_lode = M_TC_TEST

   end subroutine make_model

   ! ---------------------------------------------------------------------------
   ! Test 1: elastic step — both integrators must give identical stress.
   ! For a purely elastic increment neither integrator does any plastic work;
   ! both return the elastic predictor, so results are bit-for-bit identical.
   ! ---------------------------------------------------------------------------
   subroutine test_elastic_agreement(error)
      type(error_type), allocatable, intent(out) :: error

      type(mcsr_model_t) :: model_e, model_c
      real(wp) :: sig_e(6), sig_c(6), deps(6)
      type(integrator_params_t) :: iparams

      iparams = integrator_params_t(ftol=YIELD_TOL, stol=1.0e-4_wp, dt_min=1.0e-9_wp, max_iters=100)
      call make_model(model_e)
      call make_model(model_c)

      sig_e = SIG_INIT
      sig_c = SIG_INIT
      deps  = 1.0e-5_wp * [1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      call euler_substep(model_e, sig_e, deps, iparams)
      call cpa_step(model_c, sig_c, deps, iparams)

      call check(error, norm2(sig_e - sig_c) < 1.0e-10_wp, .true., &
                 more="elastic step: euler and cpa must give identical stress")
   end subroutine test_elastic_agreement

   ! ---------------------------------------------------------------------------
   ! Test 2: large plastic step — euler returns stress on yield surface.
   ! ---------------------------------------------------------------------------
   subroutine test_euler_on_surface(error)
      type(error_type), allocatable, intent(out) :: error

      type(mcsr_model_t) :: model
      real(wp) :: sig(6), deps(6)
      type(integrator_params_t) :: iparams

      iparams = integrator_params_t(ftol=YIELD_TOL, stol=1.0e-4_wp, dt_min=1.0e-9_wp, max_iters=100)
      call make_model(model)

      sig  = SIG_INIT
      deps = [6.0e-3_wp, 0.0_wp, -6.0e-3_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      call euler_substep(model, sig, deps, iparams)

      call check(error, model%yield_fn(sig) <= 10.0_wp * YIELD_TOL, .true., &
                 more="euler: stress should be on/inside yield surface after plastic step")
      if (allocated(error)) return

      call check(error, any(model%state%eps_p /= 0.0_wp), .true., &
                 more="euler: plastic strain should be non-zero")
   end subroutine test_euler_on_surface

   ! ---------------------------------------------------------------------------
   ! Test 3: large plastic step — cpa returns stress on yield surface.
   ! ---------------------------------------------------------------------------
   subroutine test_cpa_on_surface(error)
      type(error_type), allocatable, intent(out) :: error

      type(mcsr_model_t) :: model
      real(wp) :: sig(6), deps(6)
      type(integrator_params_t) :: iparams

      iparams = integrator_params_t(ftol=YIELD_TOL, stol=1.0e-4_wp, dt_min=1.0e-9_wp, max_iters=100)
      call make_model(model)

      sig  = SIG_INIT
      deps = [6.0e-3_wp, 0.0_wp, -6.0e-3_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      call cpa_step(model, sig, deps, iparams)

      call check(error, model%yield_fn(sig) <= 10.0_wp * YIELD_TOL, .true., &
                 more="cpa: stress should be on/inside yield surface after plastic step")
      if (allocated(error)) return

      call check(error, any(model%state%eps_p /= 0.0_wp), .true., &
                 more="cpa: plastic strain should be non-zero")
   end subroutine test_cpa_on_surface

   ! ---------------------------------------------------------------------------
   ! Test 4: small plastic overshoot — euler and cpa agree within 1%.
   ! Both methods are consistent for small increments near the yield surface.
   ! Increment is ~15% above the elastic limit so both methods are clearly
   ! plastic but the nonlinearity is small enough that they agree closely.
   ! ---------------------------------------------------------------------------
   subroutine test_stress_comparison(error)
      type(error_type), allocatable, intent(out) :: error

      type(mcsr_model_t) :: model_e, model_c
      real(wp) :: sig_e(6), sig_c(6), deps(6)
      real(wp) :: rel_diff
      type(integrator_params_t) :: iparams

      iparams = integrator_params_t(ftol=YIELD_TOL, stol=1.0e-4_wp, dt_min=1.0e-9_wp, max_iters=100)
      call make_model(model_e)
      call make_model(model_c)

      sig_e = SIG_INIT
      sig_c = SIG_INIT
      ! ~15% above elastic yield limit (q_yield = M_tc*100 = 120 kPa)
      deps = [4.0e-3_wp, 0.0_wp, -4.0e-3_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      call euler_substep(model_e, sig_e, deps, iparams)
      call cpa_step(model_c, sig_c, deps, iparams)

      rel_diff = norm2(sig_e - sig_c) / max(norm2(sig_e), 1.0e-12_wp)

      call check(error, rel_diff < 0.01_wp, .true., &
                 more="small plastic step: euler and cpa stress should agree within 1%")
   end subroutine test_stress_comparison

end module mod_test_mcsr_integration_suite
