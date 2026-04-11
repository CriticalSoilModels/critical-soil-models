!! Diagnostic sweep: integrate_stress (Euler substep vs CPA) agreement for NorSand.
!!
!! Sweeps deviatoric strain increments from just above the elastic yield limit
!! to a large overshoot. For each increment prints:
!!   - increment size (a)
!!   - estimated % overshoot of the elastic limit
!!   - relative stress difference between Euler and CPA
!!   - yield function value after each integrator
!!
!! Output is written to norsand_comparison.txt in the working directory.
!!
!! Usage:
!!   fpm run
!!
!! The elastic yield limit is:
!!   a_elastic = q_y / (2*G*sqrt(3))  where q_y = M_i * |p|
!!   At p=-100, G=10000, M_i = M_tc*3/(3+M_tc) ≈ 0.857:
!!   a_elastic ≈ 85.7 / (2*10000*1.732) ≈ 2.474e-3

program norsand_comparison
   use mod_csm_kinds,          only: wp
   use mod_norsand_model,      only: norsand_model_t, norsand_from_props, norsand_load_state
   use mod_integrate_stress,   only: integrate_stress, integrator_params_t, &
                                     INTEGRATION_EULER, INTEGRATION_CPA
   use mod_stress_invariants,  only: calc_sig_inv
   use stdlib_math,            only: linspace
   implicit none

   ! Model parameters
   real(wp), parameter :: G_0       = 10000.0_wp
   real(wp), parameter :: P_REF     = 100.0_wp
   real(wp), parameter :: NG        = 0.5_wp
   real(wp), parameter :: NU        = 0.3_wp
   real(wp), parameter :: M_TC      = 1.2_wp
   real(wp), parameter :: GAMMA     = 0.9_wp
   real(wp), parameter :: LAMBDA_C  = 0.05_wp
   real(wp), parameter :: YIELD_TOL = 1.0e-8_wp

   integer,  parameter :: N_STEPS    = 20
   real(wp), parameter :: A_MIN_FRAC = 1.01_wp   ! 1% above elastic limit
   real(wp), parameter :: A_MAX_FRAC = 2.50_wp   ! 150% above elastic limit

   integer,  parameter :: OUT_UNIT = 10
   character(len=*), parameter :: OUT_FILE = 'example/outputs/norsand_comparison.txt'

   real(wp), parameter :: SIG_INIT(6) = [-100.0_wp, -100.0_wp, -100.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

   type(norsand_model_t)     :: model_e, model_c
   type(integrator_params_t) :: iparams
   real(wp) :: props(13), statev(15)
   real(wp) :: sig_e(6), sig_c(6), deps(6)
   real(wp) :: G_init, K_init, p_i_init, e_init, M_i_init
   real(wp) :: p_init, q_init, lode_init
   real(wp) :: q_yield, a_elastic, a, overshoot_pct, rel_diff
   real(wp) :: F_e, F_c
   real(wp) :: a_vals(N_STEPS)
   integer  :: i

   ! ---------------------------------------------------------------------------
   ! Open output file
   ! ---------------------------------------------------------------------------
   open(unit=OUT_UNIT, file=OUT_FILE, status='replace', action='write')

   ! ---------------------------------------------------------------------------
   ! Compute invariants of initial stress state
   ! ---------------------------------------------------------------------------
   call calc_sig_inv(SIG_INIT, p_init, q_init, lode_init)

   ! ---------------------------------------------------------------------------
   ! Build initial model state
   ! ---------------------------------------------------------------------------
   props = [G_0, P_REF, NG, NU, 0.7_wp, GAMMA, LAMBDA_C, &
            1.0_wp, M_TC, 0.0_wp, 0.0_wp, 100.0_wp, 0.0_wp]

   G_init   = G_0 * (abs(p_init) / P_REF)**NG
   K_init   = G_init * 2.0_wp * (1.0_wp + NU) / (3.0_wp * (1.0_wp - 2.0_wp * NU))
   p_i_init = p_init / exp(1.0_wp)
   e_init   = GAMMA - LAMBDA_C * log(-p_i_init)
   M_i_init = M_TC * 3.0_wp / (3.0_wp + M_TC)

   statev    = 0.0_wp
   statev(1) = G_init
   statev(2) = K_init
   statev(3) = p_init
   statev(4) = e_init
   statev(5) = 0.0_wp
   statev(6) = 0.0_wp
   statev(7) = p_i_init
   statev(8) = M_i_init

   iparams = integrator_params_t(ftol=YIELD_TOL, stol=1.0e-4_wp, dt_min=1.0e-9_wp, max_iters=500)

   q_yield   = M_i_init * abs(p_init)
   a_elastic = q_yield / (2.0_wp * G_init * sqrt(3.0_wp))

   ! ---------------------------------------------------------------------------
   ! Pre-compute sweep values via stdlib linspace
   ! ---------------------------------------------------------------------------
   a_vals = linspace(a_elastic * A_MIN_FRAC, a_elastic * A_MAX_FRAC, N_STEPS)

   ! ---------------------------------------------------------------------------
   ! Header
   ! ---------------------------------------------------------------------------
   write(OUT_UNIT, '(a)') "NorSand: euler_substep vs cpa_step agreement sweep"
   write(OUT_UNIT, '(a,es10.3,a)') "Elastic limit a_elastic = ", a_elastic, " [-]"
   write(OUT_UNIT, '(a,f8.4,a)')   "Yield radius  q_yield   = ", q_yield, " kPa"
   write(OUT_UNIT, '(a,f8.2,a)')   "Initial p               = ", p_init,  " kPa"
   write(OUT_UNIT, '(a)')
   write(OUT_UNIT, '(a6, a12, a12, a14, a12, a12)') &
      "step", "a (x1e-3)", "overshoot%", "rel_diff", "F_euler", "F_cpa"
   write(OUT_UNIT, '(a)') repeat("-", 68)

   ! ---------------------------------------------------------------------------
   ! Sweep
   ! ---------------------------------------------------------------------------
   do i = 1, N_STEPS
      a = a_vals(i)

      model_e = norsand_from_props(props)
      call norsand_load_state(model_e, statev)
      call model_e%pre_step(SIG_INIT)

      model_c = norsand_from_props(props)
      call norsand_load_state(model_c, statev)
      call model_c%pre_step(SIG_INIT)

      sig_e = SIG_INIT
      sig_c = SIG_INIT
      deps  = [a, 0.0_wp, -a, 0.0_wp, 0.0_wp, 0.0_wp]

      call integrate_stress(model_e, sig_e, deps, INTEGRATION_EULER, iparams)
      call integrate_stress(model_c, sig_c, deps, INTEGRATION_CPA,   iparams)

      overshoot_pct = (a / a_elastic - 1.0_wp) * 100.0_wp
      rel_diff      = norm2(sig_e - sig_c) / max(norm2(sig_e), 1.0e-12_wp)
      F_e           = model_e%yield_fn(sig_e)
      F_c           = model_c%yield_fn(sig_c)

      write(OUT_UNIT, '(i6, f12.4, f12.2, es14.4, es12.3, es12.3)') &
         i, a * 1000.0_wp, overshoot_pct, rel_diff, F_e, F_c
   end do

   ! ---------------------------------------------------------------------------
   ! Close output file and notify user
   ! ---------------------------------------------------------------------------
   close(OUT_UNIT)
   write(*, '(a)') "Output written to " // OUT_FILE

end program norsand_comparison
