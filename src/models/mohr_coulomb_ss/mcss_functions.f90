! Pure constitutive functions for the Mohr-Coulomb Strain Softening model.
!
! Yield surface and plastic potential follow the smooth hyperbolic approximation
! of Abbo & Sloan (1995). All procedures are pure and operate on explicit
! (params, state, ...) arguments — no class(), no dynamic dispatch.
!
! Smoothing constants live in mcss_params_t%as_params (abbo_sloan_params_t).
! Default values (LodeT = 29.5 degrees) are in mod_mcss_types as DEFAULT_AS_PARAMS.
module mod_mcss_functions
   use mod_csm_kinds,          only: wp
   use mod_mcss_types,         only: mcss_params_t, mcss_state_t, abbo_sloan_params_t
   use mod_elastic_utils,      only: calc_stiffness_GK, calc_K_from_G_nu
   use mod_stress_invariants,  only: calc_sig_inv, calc_J_inv
   use mod_voigt_utils,        only: calc_dev_stress
   use mod_strain_invariants,  only: calc_eps_vol_inv, calc_dev_strain, calc_eps_q_inv
   use mod_stress_invar_deriv, only: calc_dp_by_dsig, calc_dJ3_by_dsig
   implicit none
   private

   public :: mcss_params_t, mcss_state_t
   public :: mcss_yield_fn
   public :: mcss_flow_rule
   public :: mcss_plastic_potential
   public :: mcss_update_hardening
   public :: mcss_elastic_stiffness
   public :: mcss_hardening_modulus

   real(wp), parameter :: INV_SQRT3   = 0.577350269189626_wp   !! 1/sqrt(3)
   real(wp), parameter :: SQRT3_OVER2 = 0.866025403784439_wp   !! sqrt(3)/2
   real(wp), parameter :: ONE_THIRD   = 1.0_wp / 3.0_wp
   real(wp), parameter :: J_TINY      = 1.0e-13_wp             !! Guard against J=0 in derivatives
   real(wp), parameter :: J0          = 0.001_wp               !! Psi interpolation range near J=0

contains

   ! ---------------------------------------------------------------------------
   ! Yield function: F = p*sin(phi) + sqrt(J^2*K^2 + a^2*sin^2(phi)) - c*cos(phi)
   ! ---------------------------------------------------------------------------
   pure function mcss_yield_fn(params, state, sig) result(F)
      type(mcss_params_t), intent(in) :: params
      type(mcss_state_t),  intent(in) :: state
      real(wp),            intent(in) :: sig(6)
      real(wp) :: F

      real(wp) :: p, q, lode, J, s3ta, K, a, a_sphi_sq

      call calc_sig_inv(sig, p, q, lode)
      J    = calc_J_inv(q)
      s3ta = sin(3.0_wp * lode)

      K        = calc_K(params%as_params, lode, s3ta, sin(state%phi))
      a        = calc_a_smooth(params%as_params, state%c, state%phi)
      a_sphi_sq = a * a * sin(state%phi) * sin(state%phi)

      F = p * sin(state%phi) + sqrt(J*J * K*K + a_sphi_sq) - state%c * cos(state%phi)
   end function mcss_yield_fn

   ! ---------------------------------------------------------------------------
   ! Flow rule: dF/dsig (normal to yield surface)
   ! ---------------------------------------------------------------------------
   pure function mcss_flow_rule(params, state, sig) result(dF_by_dsig)
      type(mcss_params_t), intent(in) :: params
      type(mcss_state_t),  intent(in) :: state
      real(wp),            intent(in) :: sig(6)
      real(wp) :: dF_by_dsig(6)

      real(wp) :: p, q, lode, J, s3ta, dev(6)

      call calc_sig_inv(sig, p, q, lode)
      J    = calc_J_inv(q)
      s3ta = sin(3.0_wp * lode)
      ! dp/dsig == dp/ds: mean stress is independent of deviatoric stress
      dev  = calc_dev_stress(sig, p)

      dF_by_dsig = calc_dF_dsig_abbo(params%as_params, J, dev, lode, s3ta, state%c, state%phi)
   end function mcss_flow_rule

   ! ---------------------------------------------------------------------------
   ! Plastic potential: dG/dsig (same form as yield, phi replaced by psi)
   ! Near J=0, psi is linearly interpolated toward phi to avoid singularity.
   ! ---------------------------------------------------------------------------
   pure function mcss_plastic_potential(params, state, sig) result(dG_by_dsig)
      type(mcss_params_t), intent(in) :: params
      type(mcss_state_t),  intent(in) :: state
      real(wp),            intent(in) :: sig(6)
      real(wp) :: dG_by_dsig(6)

      real(wp) :: p, q, lode, J, s3ta, dev(6), psi_eff

      call calc_sig_inv(sig, p, q, lode)
      J    = calc_J_inv(q)
      s3ta = sin(3.0_wp * lode)
      ! dp/dsig == dp/ds: mean stress is independent of deviatoric stress
      dev  = calc_dev_stress(sig, p)

      if (J < J0) then
         psi_eff = state%phi - J * (state%phi - state%psi) / J0
      else
         psi_eff = state%psi
      end if

      dG_by_dsig = calc_dF_dsig_abbo(params%as_params, J, dev, lode, s3ta, state%c, psi_eff)
   end function mcss_plastic_potential

   ! ---------------------------------------------------------------------------
   ! Hardening: accumulate plastic strain, then apply exponential softening
   ! ---------------------------------------------------------------------------
   pure subroutine mcss_update_hardening(params, state, deps_p)
      type(mcss_params_t), intent(in)    :: params
      type(mcss_state_t),  intent(inout) :: state
      real(wp),            intent(in)    :: deps_p(6)

      real(wp) :: eps_p_eq

      state%eps_p = state%eps_p + deps_p
      eps_p_eq    = calc_eps_q_inv(calc_dev_strain(state%eps_p, calc_eps_vol_inv(state%eps_p)))

      state%c   = params%c_res   + (params%c_peak   - params%c_res)   * exp(-params%factor * eps_p_eq)
      state%phi = params%phi_res + (params%phi_peak  - params%phi_res) * exp(-params%factor * eps_p_eq)
      state%psi = params%psi_res + (params%psi_peak  - params%psi_res) * exp(-params%factor * eps_p_eq)
   end subroutine mcss_update_hardening

   ! ---------------------------------------------------------------------------
   ! Hardening modulus: H = ∂F/∂κ · dκ/dλ
   !
   ! For the consistency condition dlambda = (n·De·dε) / (n·De·m - H).
   ! H > 0 for softening (yield surface shrinks with plastic flow).
   !
   ! Derivation:
   !   H = (∂F/∂c · dc/deps_eq + ∂F/∂φ · dφ/deps_eq) · deps_eq/dλ
   !
   !   dc/deps_eq  = -factor·(c - c_res)           [exponential softening rate]
   !   dφ/deps_eq  = -factor·(φ - φ_res)
   !   deps_eq/dλ  = calc_eps_q_inv(dev(dg/dσ))    [equivalent strain per unit λ]
   !
   !   ∂F/∂c = -cos(φ) + a²·sin²(φ) / (c·h_kernel)
   !   ∂F/∂φ = p·cos(φ) + c·sin(φ)
   !         + (J²·K·∂K/∂φ + a·sin²(φ)·∂a/∂φ + a²·sin(φ)·cos(φ)) / h_kernel
   ! ---------------------------------------------------------------------------
   pure function mcss_hardening_modulus(params, state, sig, dg_by_dsig) result(H)
      type(mcss_params_t), intent(in) :: params
      type(mcss_state_t),  intent(in) :: state
      real(wp),            intent(in) :: sig(6), dg_by_dsig(6)
      real(wp) :: H

      real(wp) :: p, q, lode, J, s3ta
      real(wp) :: dF_by_dc, dF_by_dphi
      real(wp) :: dc_by_deps_eq, dphi_by_deps_eq
      real(wp) :: deps_eq_by_dlambda, dev_dg(6)

      ! TODO: verify this hardening modulus implementation against a reference
      !       solution before using in production runs.

      call calc_sig_inv(sig, p, q, lode)
      J    = calc_J_inv(q)
      s3ta = sin(3.0_wp * lode)

      ! ∂F/∂κ — how yield surface moves with state params (c, φ)
      call calc_dF_by_dstate(params%as_params, state, p, J, lode, s3ta, dF_by_dc, dF_by_dphi)

      ! dκ/deps_eq — exponential softening rates
      call calc_dstate_by_deps_eq(params, state, dc_by_deps_eq, dphi_by_deps_eq)

      ! deps_eq/dλ — equivalent plastic strain per unit plastic multiplier
      dev_dg             = calc_dev_strain(dg_by_dsig, calc_eps_vol_inv(dg_by_dsig))
      deps_eq_by_dlambda = calc_eps_q_inv(dev_dg)

      H = (dF_by_dc * dc_by_deps_eq + dF_by_dphi * dphi_by_deps_eq) * deps_eq_by_dlambda

   end function mcss_hardening_modulus

   ! ---------------------------------------------------------------------------
   ! Elastic stiffness
   ! ---------------------------------------------------------------------------
   pure function mcss_elastic_stiffness(params, state) result(stiff_e)
      type(mcss_params_t), intent(in) :: params
      type(mcss_state_t),  intent(in) :: state
      real(wp) :: stiff_e(6,6)
      stiff_e = calc_stiffness_GK(params%G, calc_K_from_G_nu(params%G, params%nu))
   end function mcss_elastic_stiffness

   ! ===========================================================================
   ! Private helpers (Abbo & Sloan model-specific)
   ! ===========================================================================

   pure function calc_K(asp, lode, s3ta, sin_angle) result(K)
      !! Abbo-Sloan smoothing function K(θ, φ).
      !! Inner region: exact Mohr-Coulomb expression.
      !! Outer region: quadratic fit that avoids singularity at corners.
      type(abbo_sloan_params_t), intent(in) :: asp
      real(wp), intent(in) :: lode, s3ta, sin_angle
      real(wp) :: K, sgn, a_coeff, b_coeff

      if (abs(lode) < asp%lode_tr) then
         K = cos(lode) - INV_SQRT3 * sin_angle * sin(lode)
      else
         sgn     = sign(1.0_wp, lode)
         a_coeff = asp%A1 + asp%A2 * sgn * sin_angle
         b_coeff = asp%B1 * sgn + asp%B2 * sin_angle
         K       = a_coeff - b_coeff * s3ta
      end if
   end function calc_K

   pure function calc_dK_dlode(asp, lode, s3ta, sin_angle) result(dK_dlode)
      !! Derivative of Abbo-Sloan K with respect to lode angle.
      type(abbo_sloan_params_t), intent(in) :: asp
      real(wp), intent(in) :: lode, s3ta, sin_angle
      real(wp) :: dK_dlode, cos_lode, cos_3lode, sgn, b_coeff

      if (abs(lode) < asp%lode_tr) then
         cos_lode = cos(lode)
         dK_dlode = -sin(lode) - INV_SQRT3 * sin_angle * cos_lode
      else
         sgn      = sign(1.0_wp, lode)
         cos_lode = cos(lode)
         cos_3lode = cos_lode * (4.0_wp*cos_lode*cos_lode - 3.0_wp)
         b_coeff  = asp%B1 * sgn + asp%B2 * sin_angle
         dK_dlode = -3.0_wp * b_coeff * cos_3lode
      end if
   end function calc_dK_dlode

   pure function calc_a_smooth(asp, c, phi) result(a)
      !! Abbo-Sloan tip-smoothing parameter: a = smooth_coeff * c * cot(phi).
      type(abbo_sloan_params_t), intent(in) :: asp
      real(wp), intent(in) :: c, phi
      real(wp) :: a
      if (abs(phi) < 1.0e-14_wp) then
         a = asp%smooth_coeff * c
      else
         a = asp%smooth_coeff * c * cos(phi) / sin(phi)
      end if
   end function calc_a_smooth

   pure function calc_dF_dsig_abbo(asp, J, dev, lode, s3ta, c, angle) result(dF_by_dsig)
      !! Gradient of Abbo-Sloan yield/potential surface w.r.t. stress.
      !! dF/dsig = dF/dp * dp/dsig + dF/dJ * dJ/dsig + dF/dJ3 * dJ3/dsig
      !!
      !! Note: dp/dsig == dp/ds (mean stress is independent of deviatoric stress).
      type(abbo_sloan_params_t), intent(in) :: asp
      real(wp),                  intent(in) :: J, dev(6), lode, s3ta, c, angle
      real(wp) :: dF_by_dsig(6)

      real(wp) :: sin_angle, K, dK_dlode_val, a, a_sphi_sq, jk_over_H
      real(wp) :: cos_3lode, tan_3lode, j2, inv_2J
      real(wp) :: dF_dp, dF_dJ, dF_dJ3
      real(wp) :: dp_by_dsig(6), dJ_by_dsig(6), dJ3_by_dsig(6)

      sin_angle    = sin(angle)
      K            = calc_K(asp, lode, s3ta, sin_angle)
      dK_dlode_val = calc_dK_dlode(asp, lode, s3ta, sin_angle)
      a            = calc_a_smooth(asp, c, angle)
      a_sphi_sq    = a * a * sin_angle * sin_angle

      j2 = max(J * J, J_TINY)

      cos_3lode = cos(3.0_wp * lode)
      if (abs(cos_3lode) < 1.0e-14_wp) cos_3lode = sign(1.0e-14_wp, cos_3lode)
      tan_3lode = s3ta / cos_3lode

      ! dp/dsig == dp/ds: mean stress is independent of deviatoric stress
      dp_by_dsig = calc_dp_by_dsig()

      if (J > 1.0e-4_wp) then
         inv_2J = 0.5_wp / J
      else
         inv_2J = 0.0_wp
      end if
      dJ_by_dsig = [inv_2J*dev(1),         inv_2J*dev(2),         inv_2J*dev(3), &
                    2.0_wp*inv_2J*dev(4),   2.0_wp*inv_2J*dev(5),  2.0_wp*inv_2J*dev(6)]

      dJ3_by_dsig = calc_dJ3_by_dsig(dev)

      jk_over_H = J * K / sqrt(j2*K*K + a_sphi_sq)

      dF_dp  = sin_angle
      dF_dJ  = jk_over_H*K - tan_3lode * jk_over_H * dK_dlode_val
      dF_dJ3 = -SQRT3_OVER2 * dK_dlode_val * jk_over_H / (j2 * cos_3lode)

      dF_by_dsig = dF_dp*dp_by_dsig + dF_dJ*dJ_by_dsig + dF_dJ3*dJ3_by_dsig
   end function calc_dF_dsig_abbo

   ! ---------------------------------------------------------------------------
   ! ∂F/∂c and ∂F/∂φ: yield function derivatives w.r.t. state parameters.
   !
   ! ψ does not appear in F (only in G), so no ∂F/∂ψ term.
   ! If h_kernel < 1e-14 (J≈0 and a≈0), both derivatives are set to zero.
   ! ---------------------------------------------------------------------------
   pure subroutine calc_dF_by_dstate(asp, state, p, J, lode, s3ta, dF_by_dc, dF_by_dphi)
      type(abbo_sloan_params_t), intent(in)  :: asp
      type(mcss_state_t),        intent(in)  :: state
      real(wp),                  intent(in)  :: p, J, lode, s3ta
      real(wp),                  intent(out) :: dF_by_dc, dF_by_dphi

      real(wp) :: K, a, h_kernel, dK_by_dphi, da_by_dphi, sin_phi, cos_phi

      sin_phi  = sin(state%phi)
      cos_phi  = cos(state%phi)
      K        = calc_K(asp, lode, s3ta, sin_phi)
      a        = calc_a_smooth(asp, state%c, state%phi)
      h_kernel = sqrt(J*J * K*K + a*a * sin_phi*sin_phi)

      if (h_kernel < 1.0e-14_wp) then
         dF_by_dc   = 0.0_wp
         dF_by_dphi = 0.0_wp
         return
      end if

      ! ∂F/∂c = -cos(φ) + a²·sin²(φ) / (c·h_kernel)
      dF_by_dc = -cos_phi + a*a * sin_phi*sin_phi / (state%c * h_kernel)

      ! ∂K/∂φ through sin(φ) — inner vs outer Abbo-Sloan region
      if (abs(lode) < asp%lode_tr) then
         dK_by_dphi = -INV_SQRT3 * cos_phi * sin(lode)
      else
         dK_by_dphi = (asp%A2 * sign(1.0_wp, lode) - asp%B2 * s3ta) * cos_phi
      end if

      ! ∂a/∂φ: a = smooth_coeff·c·cos(φ)/sin(φ)  →  da/dφ = -smooth_coeff·c/sin²(φ)
      if (abs(state%phi) < 1.0e-14_wp) then
         da_by_dphi = 0.0_wp
      else
         da_by_dphi = -asp%smooth_coeff * state%c / (sin_phi*sin_phi)
      end if

      ! ∂F/∂φ = p·cos(φ) + c·sin(φ)
      !       + (J²·K·∂K/∂φ + a·sin²(φ)·∂a/∂φ + a²·sin(φ)·cos(φ)) / h_kernel
      dF_by_dphi = p * cos_phi + state%c * sin_phi                &
                 + (J*J * K * dK_by_dphi                          &
                    + a * sin_phi*sin_phi * da_by_dphi            &
                    + a*a * sin_phi * cos_phi) / h_kernel

   end subroutine calc_dF_by_dstate

   ! ---------------------------------------------------------------------------
   ! dc/deps_eq and dφ/deps_eq: softening rates from the exponential law.
   !   dκ/deps_eq = -factor·(κ - κ_res)
   ! ---------------------------------------------------------------------------
   pure subroutine calc_dstate_by_deps_eq(params, state, dc_by_deps_eq, dphi_by_deps_eq)
      type(mcss_params_t), intent(in)  :: params
      type(mcss_state_t),  intent(in)  :: state
      real(wp),            intent(out) :: dc_by_deps_eq, dphi_by_deps_eq

      dc_by_deps_eq   = -params%factor * (state%c   - params%c_res)
      dphi_by_deps_eq = -params%factor * (state%phi - params%phi_res)

   end subroutine calc_dstate_by_deps_eq

end module mod_mcss_functions
