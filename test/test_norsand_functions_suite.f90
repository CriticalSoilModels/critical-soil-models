!! Unit tests for `mod_norsand_functions`.
!!
!! ### Test model setup
!!
!! Parameters use N=0 (no Nova coupling) and chi_tc=0 throughout, which gives
!! M_i = M_theta independent of state. This decouples the Lode-angle interpolation
!! from state evolution and makes expected values analytically tractable.
!!
!! H_y=0 removes psi-dependence from the hardening modulus, further simplifying
!! expected values for the FD gradient check.
!!
!! ### Yield surface geometry at isotropic stress (q = 0)
!!
!! At isotropic stress SIG_ISO = [-100,-100,-100,0,0,0]:
!!   p = -100, q = 0, lode = pi/6 (library convention for q=0)
!!   M_theta = M_tc * 3/(3+M_tc)
!!
!!   F = 0 when 1 + ln(p_i/p) = 0 => p_i = p/e = -100/e
!!   F < 0 when p_i/p > 1/e  (e.g. p_i = -50 kPa)
!!   F > 0 when p_i/p < 1/e  (e.g. p_i = -20 kPa)
!!
!! ### Deviatoric test state for hardening-modulus FD check
!!
!! sig = [-100, -200, -300, 0, 0, 0]: p = -200, q = 100*sqrt(3), lode = 0 (pure shear).
!!   M_theta(lode=0) = M_tc - M_tc^2/(3+M_tc) * cos(pi/4)
!!   p_i chosen so that F = 0: p_i = p * exp(q / (|p| * M_theta) - 1)
!!   psi = 0 throughout => e = Gamma - lambda_c * ln(-p_i)

module mod_test_norsand_functions_suite
   use mod_csm_kinds,         only: wp
   use mod_norsand_types,     only: norsand_params_t, norsand_state_t
   use mod_norsand_functions, only: norsand_yield_fn, norsand_flow_rule, &
                                     norsand_plastic_potential, norsand_update_hardening, &
                                     norsand_elastic_stiffness, norsand_hardening_modulus, &
                                     norsand_pre_step
   use ieee_arithmetic,       only: ieee_is_nan
   use testdrive,             only: new_unittest, unittest_type, error_type, check
   implicit none
   private
   public :: collect_norsand_functions_suite

   real(wp), parameter :: PI = 4.0_wp * atan(1.0_wp)

   ! Model parameters — N=0 and chi_tc=0 decouple M_i from state evolution
   real(wp), parameter :: G_0_TEST      = 10000.0_wp    !! Reference shear modulus [kPa]
   real(wp), parameter :: P_REF_TEST    = 100.0_wp      !! Reference mean stress [kPa]
   real(wp), parameter :: NG_TEST       = 0.5_wp        !! Shear modulus exponent [-]
   real(wp), parameter :: NU_TEST       = 0.3_wp        !! Poisson's ratio [-]
   real(wp), parameter :: E_O_TEST      = 0.7_wp        !! Initial void ratio [-]
   real(wp), parameter :: GAMMA_TEST    = 0.9_wp        !! CSL altitude [-]
   real(wp), parameter :: LAMBDA_C_TEST = 0.05_wp       !! CSL slope [-]
   real(wp), parameter :: R_TEST        = 1.0_wp        !! OCR [-]
   real(wp), parameter :: M_TC_TEST     = 1.2_wp        !! Critical friction ratio [-]
   real(wp), parameter :: N_TEST        = 0.0_wp        !! Nova coupling (zero = standard NorSand)
   real(wp), parameter :: CHI_TC_TEST   = 0.0_wp        !! Dilatancy coefficient (zero simplifies chi_i=0)
   real(wp), parameter :: H_0_TEST      = 100.0_wp      !! Hardening modulus intercept [-]
   real(wp), parameter :: H_Y_TEST      = 0.0_wp        !! Hardening slope (zero removes psi-dependence)

   ! Isotropic stress state: p = -100 kPa, q = 0
   real(wp), parameter :: P_INIT        = -100.0_wp
   real(wp), parameter :: SIG_ISO(6)    = [P_INIT, P_INIT, P_INIT, 0.0_wp, 0.0_wp, 0.0_wp]

contains

   subroutine collect_norsand_functions_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("yield_fn: stress inside surface",     test_yield_fn_elastic),      &
         new_unittest("yield_fn: stress on surface (q=0)",   test_yield_fn_on_surface),   &
         new_unittest("yield_fn: stress outside surface",    test_yield_fn_plastic),      &
         new_unittest("yield_fn: guard at p=0",              test_yield_fn_p_zero),       &
         new_unittest("pre_step: G and K updated",           test_pre_step_moduli),       &
         new_unittest("pre_step: M_i updated",               test_pre_step_M_i),          &
         new_unittest("elastic_stiffness: non-zero diagonal",test_elastic_stiffness_diag),&
         new_unittest("flow_rule == plastic_potential",      test_associated_flow),       &
         new_unittest("flow_rule: no NaN",                   test_flow_rule_no_nan),      &
         new_unittest("hardening: p_i evolves",              test_hardening_pi_evolves),  &
         new_unittest("hardening: e evolves with vol strain",test_hardening_void_ratio),  &
         new_unittest("hardening_modulus: no NaN at q=0",   test_hmod_no_nan_hydrostatic),&
         new_unittest("hardening_modulus: FD gradient check",test_hardening_modulus_fd)   &
      ]
   end subroutine collect_norsand_functions_suite

   ! ---------------------------------------------------------------------------
   ! Helpers
   ! ---------------------------------------------------------------------------

   pure function make_params() result(p)
      !! Construct standard test parameters.
      type(norsand_params_t) :: p
      p%G_0      = G_0_TEST
      p%p_ref    = P_REF_TEST
      p%nG       = NG_TEST
      p%nu       = NU_TEST
      p%e_o      = E_O_TEST
      p%Gamma    = GAMMA_TEST
      p%lambda_c = LAMBDA_C_TEST
      p%R        = R_TEST
      p%M_tc     = M_TC_TEST
      p%N        = N_TEST
      p%chi_tc   = CHI_TC_TEST
      p%H_0      = H_0_TEST
      p%H_y      = H_Y_TEST
   end function make_params

   pure function make_state_iso(params, p_i_val) result(s)
      !! State for isotropic stress SIG_ISO, psi=0, N=0, chi_tc=0.
      !! p_i_val sets the image pressure — use P_INIT/e for on-surface.
      type(norsand_params_t), intent(in) :: params
      real(wp),               intent(in) :: p_i_val
      type(norsand_state_t) :: s

      real(wp) :: G_val, e_c

      G_val = params%G_0 * (abs(P_INIT) / params%p_ref)**params%nG
      e_c   = params%Gamma - params%lambda_c * log(-p_i_val)

      s%G           = G_val
      s%K           = G_val * 2.0_wp * (1.0_wp + params%nu) / (3.0_wp * (1.0_wp - 2.0_wp * params%nu))
      s%p           = P_INIT
      s%e           = e_c       ! psi = 0
      s%psi         = 0.0_wp
      s%chi_tce     = params%chi_tc
      s%p_i         = p_i_val
      ! N=0 => M_i = M_theta; at isotropic lode = pi/6:
      ! M_theta = M_tc - M_tc^2/(3+M_tc)*cos(0) = M_tc*3/(3+M_tc)
      s%M_i         = params%M_tc * 3.0_wp / (3.0_wp + params%M_tc)
      s%switch_yield = .false.
      s%eps_p        = 0.0_wp
   end function make_state_iso

   subroutine make_state_deviatoric(params, state, sig)
      !! Build a state on the yield surface at the pure-shear stress point
      !! sig = [-100,-200,-300,0,0,0] (p=-200, q=100*sqrt(3), lode=0).
      !! N=0, psi=0 throughout.
      type(norsand_params_t), intent(in)  :: params
      type(norsand_state_t),  intent(out) :: state
      real(wp),               intent(out) :: sig(6)

      real(wp) :: p_t, q_t, M_theta_t, p_i_t, e_c_t, G_t

      sig = [-100.0_wp, -200.0_wp, -300.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      p_t = -200.0_wp
      q_t = 100.0_wp * sqrt(3.0_wp)

      ! At lode = 0: M_theta = M_tc - M_tc^2/(3+M_tc)*cos(pi/4)
      M_theta_t = params%M_tc - params%M_tc**2 / (3.0_wp + params%M_tc) * cos(PI / 4.0_wp)

      ! F=0 => p_i = p * exp(q/(|p|*M_theta) - 1)
      p_i_t = p_t * exp(q_t / (abs(p_t) * M_theta_t) - 1.0_wp)

      ! psi=0 => e = e_c
      e_c_t = params%Gamma - params%lambda_c * log(-p_i_t)
      G_t   = params%G_0 * (abs(p_t) / params%p_ref)**params%nG

      state%p           = p_t
      state%G           = G_t
      state%K           = G_t * 2.0_wp * (1.0_wp + params%nu) / (3.0_wp * (1.0_wp - 2.0_wp * params%nu))
      state%e           = e_c_t
      state%psi         = 0.0_wp
      state%chi_tce     = params%chi_tc
      state%p_i         = p_i_t
      state%M_i         = M_theta_t    ! N=0
      state%switch_yield = .false.
      state%eps_p        = 0.0_wp
   end subroutine make_state_deviatoric

   ! ---------------------------------------------------------------------------
   ! Yield function tests
   ! ---------------------------------------------------------------------------

   subroutine test_yield_fn_elastic(error)
      !! p_i = -50 kPa > P_INIT/e at isotropic state: p_i/p = 0.5 > 1/e, F < 0.
      type(error_type), allocatable, intent(out) :: error

      type(norsand_params_t) :: params
      type(norsand_state_t)  :: state
      real(wp) :: F

      params = make_params()
      state  = make_state_iso(params, -50.0_wp)
      F      = norsand_yield_fn(params, state, SIG_ISO)

      call check(error, F < 0.0_wp, .true., more="elastic: F should be < 0 (inside yield surface)")
   end subroutine test_yield_fn_elastic

   subroutine test_yield_fn_on_surface(error)
      !! At isotropic state with p_i = P_INIT/e, 1+ln(p_i/p)=0 => F=0 exactly.
      type(error_type), allocatable, intent(out) :: error

      type(norsand_params_t) :: params
      type(norsand_state_t)  :: state
      real(wp) :: F, p_i_on

      params   = make_params()
      p_i_on   = P_INIT / exp(1.0_wp)      ! p_i/p = 1/e => 1 + ln(p_i/p) = 0
      state    = make_state_iso(params, p_i_on)
      F        = norsand_yield_fn(params, state, SIG_ISO)

      call check(error, abs(F) < 1.0e-10_wp, .true., &
                 more="on-surface: |F| should be < 1e-10 kPa at constructed on-surface state")
   end subroutine test_yield_fn_on_surface

   subroutine test_yield_fn_plastic(error)
      !! p_i = -20 kPa < P_INIT/e at isotropic state: p_i/p = 0.2 < 1/e, F > 0.
      type(error_type), allocatable, intent(out) :: error

      type(norsand_params_t) :: params
      type(norsand_state_t)  :: state
      real(wp) :: F

      params = make_params()
      state  = make_state_iso(params, -20.0_wp)
      F      = norsand_yield_fn(params, state, SIG_ISO)

      call check(error, F > 0.0_wp, .true., more="plastic: F should be > 0 (outside yield surface)")
   end subroutine test_yield_fn_plastic

   subroutine test_yield_fn_p_zero(error)
      !! Near-zero p: guard should return F = 0 without NaN or division by zero.
      type(error_type), allocatable, intent(out) :: error

      type(norsand_params_t) :: params
      type(norsand_state_t)  :: state
      real(wp) :: sig_near_zero(6), F

      params = make_params()
      state  = make_state_iso(params, P_INIT / exp(1.0_wp))

      sig_near_zero = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      F = norsand_yield_fn(params, state, sig_near_zero)

      call check(error, .not. ieee_is_nan(F), .true., more="p=0 guard: F must not be NaN")
      if (allocated(error)) return
      call check(error, F == 0.0_wp, .true., more="p=0 guard: F should be 0 when |p| < tolerance")
   end subroutine test_yield_fn_p_zero

   ! ---------------------------------------------------------------------------
   ! pre_step tests
   ! ---------------------------------------------------------------------------

   subroutine test_pre_step_moduli(error)
      !! pre_step must update G = G_0*(|p|/p_ref)^nG and K from G and nu.
      !! At p=-100 and p_ref=100: G_new = G_0 * 1^0.5 = G_0.
      !! At p=-200: G_new = G_0 * (200/100)^0.5 = G_0 * sqrt(2).
      type(error_type), allocatable, intent(out) :: error

      type(norsand_params_t) :: params
      type(norsand_state_t)  :: state
      real(wp) :: sig_p200(6), G_expected, K_expected
      real(wp), parameter :: TOL = 1.0e-10_wp

      params     = make_params()
      state      = make_state_iso(params, P_INIT / exp(1.0_wp))
      sig_p200   = [-200.0_wp, -200.0_wp, -200.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      call norsand_pre_step(params, state, sig_p200)

      G_expected = G_0_TEST * sqrt(2.0_wp)         ! (200/100)^0.5
      K_expected = G_expected * 2.0_wp * (1.0_wp + NU_TEST) / (3.0_wp * (1.0_wp - 2.0_wp * NU_TEST))

      call check(error, abs(state%G - G_expected) < TOL * G_expected, .true., &
                 more="pre_step: G not updated to G_0*(|p|/p_ref)^nG")
      if (allocated(error)) return
      call check(error, abs(state%K - K_expected) < TOL * K_expected, .true., &
                 more="pre_step: K not updated consistently with G")
      if (allocated(error)) return
      call check(error, abs(state%p - (-200.0_wp)) < TOL, .true., &
                 more="pre_step: state%p not updated to current mean stress")
   end subroutine test_pre_step_moduli

   subroutine test_pre_step_M_i(error)
      !! At lode = pi/6 (isotropic, q=0) with N=0:
      !! M_theta = M_tc * 3/(3+M_tc)
      !! After pre_step, state%M_i must equal M_theta.
      type(error_type), allocatable, intent(out) :: error

      type(norsand_params_t) :: params
      type(norsand_state_t)  :: state
      real(wp) :: M_theta_expected
      real(wp), parameter :: TOL = 1.0e-10_wp

      params = make_params()
      state  = make_state_iso(params, P_INIT / exp(1.0_wp))

      ! Zero out M_i to confirm pre_step sets it correctly
      state%M_i = 0.0_wp
      call norsand_pre_step(params, state, SIG_ISO)

      M_theta_expected = M_TC_TEST * 3.0_wp / (3.0_wp + M_TC_TEST)

      call check(error, abs(state%M_i - M_theta_expected) < TOL, .true., &
                 more="pre_step: M_i not updated to M_theta for N=0 at isotropic state")
   end subroutine test_pre_step_M_i

   ! ---------------------------------------------------------------------------
   ! Elastic stiffness test
   ! ---------------------------------------------------------------------------

   subroutine test_elastic_stiffness_diag(error)
      !! The elastic stiffness matrix must have positive diagonal entries.
      type(error_type), allocatable, intent(out) :: error

      type(norsand_params_t) :: params
      type(norsand_state_t)  :: state
      real(wp) :: D(6,6)
      integer  :: i

      params = make_params()
      state  = make_state_iso(params, P_INIT / exp(1.0_wp))
      D      = norsand_elastic_stiffness(params, state)

      do i = 1, 6
         call check(error, D(i,i) > 0.0_wp, .true., &
                    more="elastic_stiffness: non-positive diagonal entry")
         if (allocated(error)) return
      end do
   end subroutine test_elastic_stiffness_diag

   ! ---------------------------------------------------------------------------
   ! Flow rule tests (associated plasticity)
   ! ---------------------------------------------------------------------------

   subroutine test_associated_flow(error)
      !! For NorSand, flow rule == plastic potential (associated plasticity).
      !! Verify they return identical vectors at a non-trivial stress state.
      type(error_type), allocatable, intent(out) :: error

      type(norsand_params_t) :: params
      type(norsand_state_t)  :: state
      real(wp) :: sig(6), dF(6), dG(6)
      real(wp), parameter :: TOL = 0.0_wp   ! must be bit-for-bit identical

      params = make_params()
      call make_state_deviatoric(params, state, sig)

      dF = norsand_flow_rule(params, state, sig)
      dG = norsand_plastic_potential(params, state, sig)

      call check(error, all(dF == dG), .true., &
                 more="associated flow: flow_rule and plastic_potential must be identical")
   end subroutine test_associated_flow

   subroutine test_flow_rule_no_nan(error)
      !! dF/dsig must not produce NaN at a general stress state or at q=0.
      type(error_type), allocatable, intent(out) :: error

      type(norsand_params_t) :: params
      type(norsand_state_t)  :: state
      real(wp) :: n(6), sig_dev(6)

      params  = make_params()
      state   = make_state_iso(params, -50.0_wp)
      sig_dev = [-100.0_wp, -200.0_wp, -300.0_wp, 50.0_wp, 30.0_wp, 20.0_wp]

      ! General deviatoric stress
      n = norsand_flow_rule(params, state, sig_dev)
      call check(error, all(.not. ieee_is_nan(n)), .true., &
                 more="flow_rule: NaN at general stress state")
      if (allocated(error)) return

      ! Isotropic (q = 0)
      n = norsand_flow_rule(params, state, SIG_ISO)
      call check(error, all(.not. ieee_is_nan(n)), .true., &
                 more="flow_rule: NaN at isotropic (q=0) state")
   end subroutine test_flow_rule_no_nan

   ! ---------------------------------------------------------------------------
   ! Hardening tests
   ! ---------------------------------------------------------------------------

   subroutine test_hardening_pi_evolves(error)
      !! After a deviatoric plastic increment, p_i must change from its initial value.
      !! Uses the deviatoric test state which is on the yield surface.
      type(error_type), allocatable, intent(out) :: error

      type(norsand_params_t) :: params
      type(norsand_state_t)  :: state
      real(wp) :: sig(6), dg(6), p_i_before
      real(wp), parameter :: DEPS_SCALE = 1.0e-5_wp

      params = make_params()
      call make_state_deviatoric(params, state, sig)
      p_i_before = state%p_i

      dg = norsand_plastic_potential(params, state, sig)
      dg = DEPS_SCALE * dg
      call norsand_update_hardening(params, state, dg)

      call check(error, state%p_i /= p_i_before, .true., &
                 more="hardening: p_i should evolve after a plastic strain increment")
   end subroutine test_hardening_pi_evolves

   subroutine test_hardening_void_ratio(error)
      !! A purely volumetric plastic increment must update the void ratio.
      !! deps_p = eps_vol * [1, 1, 1, 0, 0, 0] / 3 with eps_vol < 0 (compaction).
      !! Expected: e_new = e_old + eps_vol * (1 + e_old).
      type(error_type), allocatable, intent(out) :: error

      type(norsand_params_t) :: params
      type(norsand_state_t)  :: state
      real(wp) :: e_before, e_expected, deps_vol(6)
      real(wp), parameter :: EPS_VOL    = -1.0e-4_wp   ! compaction
      real(wp), parameter :: TOL        = 1.0e-12_wp

      params   = make_params()
      state    = make_state_iso(params, P_INIT / exp(1.0_wp))
      e_before = state%e
      deps_vol = [EPS_VOL/3.0_wp, EPS_VOL/3.0_wp, EPS_VOL/3.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      call norsand_update_hardening(params, state, deps_vol)

      e_expected = e_before + EPS_VOL * (1.0_wp + e_before)

      call check(error, abs(state%e - e_expected) < TOL, .true., &
                 more="hardening: void ratio e should update with volumetric plastic strain")
   end subroutine test_hardening_void_ratio

   ! ---------------------------------------------------------------------------
   ! Hardening modulus tests
   ! ---------------------------------------------------------------------------

   subroutine test_hmod_no_nan_hydrostatic(error)
      !! At q=0, dq/dsig = 0 so the flow direction is purely hydrostatic.
      !! The deviatoric part of dg is zero, giving deps_eq_by_dlambda = 0 => H = 0.
      !! No NaN should arise from the p ≈ 0 or p_i ≈ 0 guards.
      type(error_type), allocatable, intent(out) :: error

      type(norsand_params_t) :: params
      type(norsand_state_t)  :: state
      real(wp) :: dg(6), H

      params = make_params()
      state  = make_state_iso(params, P_INIT / exp(1.0_wp))
      dg     = norsand_plastic_potential(params, state, SIG_ISO)
      H      = norsand_hardening_modulus(params, state, SIG_ISO, dg)

      call check(error, .not. ieee_is_nan(H), .true., &
                 more="hardening_modulus: NaN at isotropic (q=0) state")
   end subroutine test_hmod_no_nan_hydrostatic

   subroutine test_hardening_modulus_fd(error)
      !! Finite-difference gradient check for norsand_hardening_modulus.
      !!
      !! Uses the deviatoric test state on the yield surface (lode=0, psi=0, N=0).
      !! Verifies that the analytical H matches a forward-difference approximation:
      !!
      !!   H_fd = (F(sig, state_perturbed) - F(sig, state_0)) / h
      !!
      !! where state_perturbed is obtained by calling norsand_update_hardening with
      !! h * dg_by_dsig (one unit of plastic multiplier scaled by h). Since F_0 ≈ 0
      !! (on surface), H_fd ≈ F_plus / h. The FD step h = 1e-7 gives relative
      !! truncation error well within the 1e-4 tolerance for smooth functions.
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: FD_STEP    = 1.0e-7_wp
      real(wp), parameter :: FD_TOL     = 1.0e-4_wp
      real(wp), parameter :: FD_REF_MIN = 1.0_wp

      type(norsand_params_t) :: params
      type(norsand_state_t)  :: state_0, state_plus
      real(wp) :: sig(6), dg(6)
      real(wp) :: H_analytical, H_fd, F_0, F_plus, rel_err

      params = make_params()
      call make_state_deviatoric(params, state_0, sig)

      ! Plastic flow direction at test state
      dg = norsand_plastic_potential(params, state_0, sig)

      ! Analytical hardening modulus
      H_analytical = norsand_hardening_modulus(params, state_0, sig, dg)

      ! Yield function on surface (should be ≈ 0 by construction)
      F_0 = norsand_yield_fn(params, state_0, sig)

      ! Forward-difference: apply FD_STEP * dg as plastic strain increment
      state_plus = state_0
      dg = FD_STEP * dg
      call norsand_update_hardening(params, state_plus, dg)
      F_plus = norsand_yield_fn(params, state_plus, sig)

      ! H = (F_plus - F_0) / FD_STEP  (codebase sign: H > 0 means softening)
      H_fd = (F_plus - F_0) / FD_STEP

      rel_err = abs(H_analytical - H_fd) / max(abs(H_analytical), FD_REF_MIN)

      call check(error, rel_err < FD_TOL, .true., &
                 more="hardening_modulus FD: analytical and finite-difference H disagree")
   end subroutine test_hardening_modulus_fd

end module mod_test_norsand_functions_suite
