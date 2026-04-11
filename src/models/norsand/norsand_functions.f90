! Pure constitutive functions for the NorSand model (Jefferies 1993).
!
! Yield surface (outer cap, compression negative):
!   F = q + p * M_i * (1 + ln(p_i / p))
!
! Lode-angle interpolation (Jefferies & Shuttle 2011):
!   M_theta = M_tc - M_tc^2/(3+M_tc) * cos(-3*theta/2 + pi/4)
!
! All procedures are pure (except norsand_pre_step) and operate on explicit
! (params, state, ...) arguments — no class(), no dynamic dispatch.
module mod_norsand_functions
   use mod_csm_kinds,         only: wp
   use mod_norsand_types,     only: norsand_params_t, norsand_state_t
   use mod_elastic_utils,     only: calc_stiffness_GK, calc_K_from_G_nu
   use mod_stress_invariants, only: calc_sig_inv, calc_J2_inv, calc_J3_inv
   use mod_stress_invar_deriv, only: calc_dp_by_dsig, calc_dJ_by_dsig, calc_dJ3_by_dsig, &
                                     calc_dlode_angle_by_dsig, calc_dJ2_by_dsig
   use mod_voigt_utils,       only: calc_dev_stress
   use mod_strain_invariants, only: calc_eps_vol_inv, calc_dev_strain, calc_eps_q_inv
   implicit none
   private

   public :: norsand_params_t, norsand_state_t
   public :: norsand_yield_fn
   public :: norsand_flow_rule
   public :: norsand_plastic_potential
   public :: norsand_elastic_stiffness
   public :: norsand_update_hardening
   public :: norsand_hardening_modulus
   public :: norsand_pre_step

   real(wp), parameter :: P_TINY      = 1.0e-12_wp  !! Guard for near-zero mean stress
   real(wp), parameter :: M_ITC_TINY  = 1.0e-14_wp  !! Guard for near-zero M_itc
   real(wp), parameter :: J_TINY      = 1.0e-4_wp   !! Guard for near-zero J in gradient
   real(wp), parameter :: PI          = 4.0_wp * atan(1.0_wp)
   real(wp), parameter :: LODE_MIN    = -PI / 6.0_wp
   real(wp), parameter :: LODE_MAX    =  PI / 6.0_wp

contains

   ! ---------------------------------------------------------------------------
   ! Yield function: F = q + p * M_i * (1 + ln(p_i / p))
   ! ---------------------------------------------------------------------------
   pure function norsand_yield_fn(params, state, sig) result(F)
      !! NorSand yield function (outer cap, Jefferies 1993).
      !! Returns 0 if |p| or |p_i| are below guard tolerance.
      type(norsand_params_t), intent(in) :: params
      type(norsand_state_t),  intent(in) :: state
      real(wp),               intent(in) :: sig(6)
      real(wp) :: F

      real(wp) :: p, q, lode, M_theta, chi_i, M_itc, M_i

      call calc_sig_inv(sig, p, q, lode)

      if (abs(p) < P_TINY .or. abs(state%p_i) < P_TINY) then
         F = 0.0_wp
         return
      end if

      M_theta = calc_M_theta(params%M_tc, lode)
      chi_i   = calc_chi_i(state%chi_tce, params%lambda_c, params%M_tc)
      M_itc   = calc_M_itc(params%M_tc, params%N, chi_i, state%psi)
      M_i     = M_theta * (1.0_wp - params%N * chi_i * abs(state%psi) / params%M_tc)

      F = q + p * M_i * (1.0_wp + log(state%p_i / p))

   end function norsand_yield_fn

   ! ---------------------------------------------------------------------------
   ! Flow rule: dF/dsig (associated plasticity)
   ! ---------------------------------------------------------------------------
   pure function norsand_flow_rule(params, state, sig) result(dF_by_dsig)
      !! NorSand flow rule — associated, delegates to shared gradient helper.
      type(norsand_params_t), intent(in) :: params
      type(norsand_state_t),  intent(in) :: state
      real(wp),               intent(in) :: sig(6)
      real(wp) :: dF_by_dsig(6)

      dF_by_dsig = norsand_dF_by_dsig(params, state, sig)

   end function norsand_flow_rule

   ! ---------------------------------------------------------------------------
   ! Plastic potential: identical to flow rule (associated)
   ! ---------------------------------------------------------------------------
   pure function norsand_plastic_potential(params, state, sig) result(dG_by_dsig)
      !! NorSand plastic potential gradient — associated, same as flow rule.
      type(norsand_params_t), intent(in) :: params
      type(norsand_state_t),  intent(in) :: state
      real(wp),               intent(in) :: sig(6)
      real(wp) :: dG_by_dsig(6)

      dG_by_dsig = norsand_dF_by_dsig(params, state, sig)

   end function norsand_plastic_potential

   ! ---------------------------------------------------------------------------
   ! Elastic stiffness: pressure-dependent, uses state%G and state%K
   ! ---------------------------------------------------------------------------
   pure function norsand_elastic_stiffness(params, state) result(stiff_e)
      !! Isotropic elastic stiffness from current G and K stored in state.
      type(norsand_params_t), intent(in) :: params
      type(norsand_state_t),  intent(in) :: state
      real(wp) :: stiff_e(6,6)

      stiff_e = calc_stiffness_GK(state%G, state%K)

   end function norsand_elastic_stiffness

   ! ---------------------------------------------------------------------------
   ! Hardening subroutine: update state given plastic strain increment deps_p(6)
   ! ---------------------------------------------------------------------------
   pure subroutine norsand_update_hardening(params, state, deps_p)
      !! Update NorSand hardening variables in-place.
      !!
      !! Updates: p_i, eps_p, e, psi, M_i.
      !! Uses state%p (set by pre_step) for hardening rate computation.
      type(norsand_params_t), intent(in)    :: params
      type(norsand_state_t),  intent(inout) :: state
      real(wp),               intent(in)    :: deps_p(6)

      real(wp) :: deps_p_eq, eps_vol_p
      real(wp) :: chi_i_old, M_itc_old, M_theta
      real(wp) :: H_eff, ratio, dp_i_deps
      real(wp) :: e_c, chi_i_new, M_itc_new
      real(wp) :: dev_deps_p(6), deps_vol_p

      ! Step 1: Equivalent deviatoric plastic strain increment
      deps_vol_p = calc_eps_vol_inv(deps_p)
      dev_deps_p    = calc_dev_strain(deps_p, deps_vol_p)
      deps_p_eq     = calc_eps_q_inv(dev_deps_p)

      ! Step 2: Recover M_theta from stored M_i using pre-hardening quantities
      chi_i_old = calc_chi_i(state%chi_tce, params%lambda_c, params%M_tc)
      M_itc_old = calc_M_itc(params%M_tc, params%N, chi_i_old, state%psi)

      if (M_itc_old > M_ITC_TINY) then
         M_theta = state%M_i / (1.0_wp - params%N * chi_i_old * abs(state%psi) / params%M_tc)
      else
         M_theta = state%M_i
      end if

      ! Step 3: Effective hardening modulus
      H_eff = params%H_0 - params%H_y * state%psi

      ! Step 4: Hardening rate dp_i/deps_eq
      if (abs(state%p) < P_TINY .or. abs(state%p_i) < P_TINY .or. M_itc_old < M_ITC_TINY) then
         dp_i_deps = 0.0_wp
      else
         ratio     = state%p / state%p_i
         dp_i_deps = H_eff * state%p_i * state%M_i * ratio**2 &
                   * (exp(-chi_i_old * state%psi / M_itc_old) - state%p_i / state%p) &
                   / M_itc_old
      end if

      ! Step 5: Update state variables
      state%p_i  = state%p_i + dp_i_deps * deps_p_eq
      state%eps_p = state%eps_p + deps_p

      ! Void ratio evolves with volumetric plastic strain
      eps_vol_p    = deps_p(1) + deps_p(2) + deps_p(3)
      state%e      = state%e + eps_vol_p * (1.0_wp + state%e)

      ! Guard: keep p_i negative (compression negative convention)
      if (-state%p_i < P_TINY) state%p_i = -P_TINY

      ! Update critical void ratio and state parameter
      e_c       = params%Gamma - params%lambda_c * log(-state%p_i)
      state%psi = state%e - e_c

      ! Update chi_i and M_i with new state
      chi_i_new   = calc_chi_i(state%chi_tce, params%lambda_c, params%M_tc)
      M_itc_new   = calc_M_itc(params%M_tc, params%N, chi_i_new, state%psi)
      state%M_i   = M_theta * (1.0_wp - params%N * chi_i_new * abs(state%psi) / params%M_tc)

   end subroutine norsand_update_hardening

   ! ---------------------------------------------------------------------------
   ! Hardening modulus: H = (∂F/∂p_i) * (dp_i/deps_eq) * (deps_eq/dλ)
   ! ---------------------------------------------------------------------------
   pure function norsand_hardening_modulus(params, state, sig, dg_by_dsig) result(H)
      !! Hardening modulus for the consistency condition.
      !! H > 0 corresponds to softening (shrinking yield surface) in the
      !! codebase sign convention.
      type(norsand_params_t), intent(in) :: params
      type(norsand_state_t),  intent(in) :: state
      real(wp),               intent(in) :: sig(6), dg_by_dsig(6)
      real(wp) :: H

      real(wp) :: p, q, lode
      real(wp) :: chi_i, M_itc, ratio
      real(wp) :: dF_dp_i, dp_i_deps, deps_eq_by_dlambda
      real(wp) :: H_eff, dev_dg(6), eps_vol_dg

      call calc_sig_inv(sig, p, q, lode)

      if (abs(p) < P_TINY .or. abs(state%p_i) < P_TINY) then
         H = 0.0_wp
         return
      end if

      chi_i = calc_chi_i(state%chi_tce, params%lambda_c, params%M_tc)
      M_itc = calc_M_itc(params%M_tc, params%N, chi_i, state%psi)

      if (M_itc < M_ITC_TINY) then
         H = 0.0_wp
         return
      end if

      ! ∂F/∂p_i = M_i * (p / p_i)
      dF_dp_i = state%M_i * (p / state%p_i)

      ! dp_i/deps_eq — same formula as in update_hardening
      H_eff = params%H_0 - params%H_y * state%psi
      ratio = state%p / state%p_i
      dp_i_deps = H_eff * state%p_i * state%M_i * ratio**2 &
                * (exp(-chi_i * state%psi / M_itc) - state%p_i / state%p) &
                / M_itc

      ! deps_eq/dλ = calc_eps_q_inv(dev(dg_by_dsig))
      eps_vol_dg         = calc_eps_vol_inv(dg_by_dsig)
      dev_dg             = calc_dev_strain(dg_by_dsig, eps_vol_dg)
      deps_eq_by_dlambda = calc_eps_q_inv(dev_dg)

      H = dF_dp_i * dp_i_deps * deps_eq_by_dlambda

   end function norsand_hardening_modulus

   ! ---------------------------------------------------------------------------
   ! Pre-step: update pressure-dependent G, K, p, and M_i from current stress
   ! NOT pure — modifies state
   ! ---------------------------------------------------------------------------
   subroutine norsand_pre_step(params, state, sig)
      !! Update stress-dependent state variables before each stress integration step.
      !! Called via csm_model_t%pre_step interface.
      type(norsand_params_t), intent(in)    :: params
      type(norsand_state_t),  intent(inout) :: state
      real(wp),               intent(in)    :: sig(6)

      real(wp) :: p, q, lode, G_new, M_theta, chi_i, M_itc

      call calc_sig_inv(sig, p, q, lode)
      state%p = p

      ! Pressure-dependent shear modulus: G = G_0 * (|p| / p_ref)^nG
      if (abs(p) < P_TINY) then
         G_new = params%G_0
      else
         G_new = params%G_0 * (abs(p) / params%p_ref)**params%nG
      end if
      state%G = G_new
      state%K = calc_K_from_G_nu(G_new, params%nu)

      ! Update M_i with current Lode angle
      M_theta   = calc_M_theta(params%M_tc, lode)
      chi_i     = calc_chi_i(state%chi_tce, params%lambda_c, params%M_tc)
      M_itc     = calc_M_itc(params%M_tc, params%N, chi_i, state%psi)
      state%M_i = M_theta * (1.0_wp - params%N * chi_i * abs(state%psi) / params%M_tc)

   end subroutine norsand_pre_step

   ! ===========================================================================
   ! Private helpers
   ! ===========================================================================

   pure function calc_M_theta(M_tc, lode) result(M_theta)
      !! Lode-angle interpolation of critical friction ratio.
      !! M_theta = M_tc - M_tc^2/(3+M_tc) * cos(-3*lode/2 + pi/4)
      !! Lode angle is clamped to [-pi/6, pi/6].
      real(wp), intent(in) :: M_tc   !! Critical friction ratio (triaxial compression) [-]
      real(wp), intent(in) :: lode   !! Lode angle [rad]
      real(wp) :: M_theta

      real(wp) :: lode_clamped

      lode_clamped = max(LODE_MIN, min(LODE_MAX, lode))
      M_theta = M_tc - M_tc**2 / (3.0_wp + M_tc) * cos(-1.5_wp * lode_clamped + PI / 4.0_wp)

   end function calc_M_theta

   pure function calc_dM_theta_by_dlode(M_tc, lode) result(dM_theta_by_dlode)
      !! Derivative of M_theta with respect to Lode angle.
      !! dM_theta/dlode = - M_tc^2/(3+M_tc) * sin(-3*lode/2 + pi/4) * (-3/2)
      !!                = 3/2 * M_tc^2/(3+M_tc) * sin(-3*lode/2 + pi/4)
      real(wp), intent(in) :: M_tc  !! Critical friction ratio [-]
      real(wp), intent(in) :: lode  !! Lode angle [rad]
      real(wp) :: dM_theta_by_dlode

      real(wp) :: lode_clamped

      lode_clamped = max(LODE_MIN, min(LODE_MAX, lode))
      dM_theta_by_dlode = 1.5_wp * M_tc**2 / (3.0_wp + M_tc) &
                        * sin(-1.5_wp * lode_clamped + PI / 4.0_wp)

   end function calc_dM_theta_by_dlode

   pure function calc_chi_i(chi_tce, lambda_c, M_tc) result(chi_i)
      !! CHI_i = chi_tce / (1 - lambda_c * chi_tce / M_tc)
      real(wp), intent(in) :: chi_tce   !! Current dilatancy coefficient [-]
      real(wp), intent(in) :: lambda_c  !! CSL slope [-]
      real(wp), intent(in) :: M_tc      !! Critical friction ratio [-]
      real(wp) :: chi_i

      real(wp) :: denom

      denom = 1.0_wp - lambda_c * chi_tce / M_tc
      if (abs(denom) < 1.0e-14_wp) denom = sign(1.0e-14_wp, denom)
      chi_i = chi_tce / denom

   end function calc_chi_i

   pure function calc_M_itc(M_tc, N, chi_i, psi) result(M_itc)
      !! M_itc = M_tc * (1 - N * chi_i * |psi| / M_tc)
      real(wp), intent(in) :: M_tc   !! Critical friction ratio [-]
      real(wp), intent(in) :: N      !! Nova volumetric coupling coefficient [-]
      real(wp), intent(in) :: chi_i  !! Dilatancy coefficient CHI_i [-]
      real(wp), intent(in) :: psi    !! State parameter [-]
      real(wp) :: M_itc

      M_itc = M_tc * (1.0_wp - N * chi_i * abs(psi) / M_tc)

   end function calc_M_itc

   pure function norsand_dF_by_dsig(params, state, sig) result(dF_by_dsig)
      !! Gradient ∂F/∂σ for F = q + p * M_i * (1 + ln(p_i/p)).
      !!
      !! ∂F/∂σ = (∂F/∂p) * ∂p/∂σ + (∂F/∂q) * ∂q/∂σ + (∂F/∂θ) * ∂θ/∂σ
      !!
      !! where q = sqrt(3*J2) and J = q/sqrt(3) = sqrt(J2).
      !!
      !! ∂F/∂p  = M_i * ln(p_i/p)
      !! ∂F/∂q  = 1
      !! ∂F/∂θ  = p * (1 + ln(p_i/p)) * (∂M_i/∂M_theta) * dM_theta/dθ
      !!
      !! The Lode derivative ∂θ/∂σ is computed via calc_dlode_angle_by_dsig.
      !! The deviatoric stress gradient ∂q/∂σ = (3/(2q)) * ∂J2/∂σ = (3/(2q)) * dJ2_dsig.
      type(norsand_params_t), intent(in) :: params
      type(norsand_state_t),  intent(in) :: state
      real(wp),               intent(in) :: sig(6)
      real(wp) :: dF_by_dsig(6)

      real(wp) :: p, q, lode
      real(wp) :: dev(6), J2, J3
      real(wp) :: M_theta, chi_i, M_itc, M_i
      real(wp) :: dF_dp, dF_dq, dF_dtheta
      real(wp) :: dM_i_by_dM_theta, dM_theta_dlode
      real(wp) :: dp_dsig(6), dq_dsig(6), dlode_dsig(6)
      real(wp) :: dJ3_dsig(6), dJ2_dsig(6)

      call calc_sig_inv(sig, p, q, lode)

      if (abs(p) < P_TINY .or. abs(state%p_i) < P_TINY) then
         dF_by_dsig = 0.0_wp
         return
      end if

      dev = calc_dev_stress(sig, p)

      ! J2 = (1/2) tr(s^2), J3 = det(s)
      J2 = calc_J2_inv(dev)

      ! Lode-dependent quantities
      M_theta = calc_M_theta(params%M_tc, lode)
      chi_i   = calc_chi_i(state%chi_tce, params%lambda_c, params%M_tc)
      M_itc   = calc_M_itc(params%M_tc, params%N, chi_i, state%psi)
      M_i     = M_theta * (1.0_wp - params%N * chi_i * abs(state%psi) / params%M_tc)

      ! Partial derivatives of F w.r.t. invariants
      dF_dp = M_i * log(state%p_i / p)
      dF_dq = 1.0_wp

      ! ∂M_i/∂M_theta = (1 - N * chi_i * |psi| / M_tc)
      dM_i_by_dM_theta  = 1.0_wp - params%N * chi_i * abs(state%psi) / params%M_tc
      dM_theta_dlode     = calc_dM_theta_by_dlode(params%M_tc, lode)

      ! ∂F/∂θ = ∂F/∂M_i * ∂M_i/∂M_theta * dM_theta/dθ
      !       = p*(1+ln(p_i/p)) * dM_i_by_dM_theta * dM_theta_dlode
      dF_dtheta = p * (1.0_wp + log(state%p_i / p)) * dM_i_by_dM_theta * dM_theta_dlode

      ! Stress-invariant gradients
      dp_dsig = calc_dp_by_dsig()

      ! ∂q/∂σ = 3/(2q) * ∂J2/∂σ   (guard q ≈ 0)
      dJ2_dsig = calc_dJ2_by_dsig(dev)
      if (q > J_TINY) then
         dq_dsig = (3.0_wp / (2.0_wp * q)) * dJ2_dsig
      else
         dq_dsig = 0.0_wp
      end if

      ! ∂θ/∂σ via existing helper (requires J3 explicitly)
      J3 = calc_J3_inv(dev)
      dJ3_dsig  = calc_dJ3_by_dsig(dev)
      dlode_dsig = calc_dlode_angle_by_dsig(dJ3_dsig, dev, J3, max(J2, P_TINY), lode)

      dF_by_dsig = dF_dp * dp_dsig + dF_dq * dq_dsig + dF_dtheta * dlode_dsig

   end function norsand_dF_by_dsig


end module mod_norsand_functions
