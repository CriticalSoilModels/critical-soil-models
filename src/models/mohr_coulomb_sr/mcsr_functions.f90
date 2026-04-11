!! Pure constitutive functions for the MCSR (Mohr-Coulomb Strain Rate) model.
!!
!! All mathematics lives here — no class(), no dynamic dispatch, no I/O.
!! The functions operate on explicit (params, state, ...) arguments and are
!! suitable for GPU kernels via OpenACC/CUDA Fortran.
!!
!! ### Rate-state architecture
!!
!! Rate-dependent parameters (G, K, eta_y, dilation, M_lode) are pre-computed
!! **once per timestep** by `mcsr_update_rate_state` before the integrator runs.
!! Inside `euler_substep` those values are frozen; only dilation evolves with
!! accumulated plastic strain (via `mcsr_update_hardening`).
!!
!! ### Yield function
!!
!! F = q + eta_y * p,  where eta_y = M_lode - dilation*(1 - N)
!! and M_lode = M_tc * (1 + 0.25 * cos(1.5*theta + 0.25*pi)^1.2)
!!
!! Sign convention: compression is negative (UMAT standard).

module mod_mcsr_functions
   use mod_csm_kinds,          only: wp
   use mod_mcsr_types,         only: mcsr_params_t, mcsr_state_t
   use mod_elastic_utils,      only: calc_stiffness_GK
   use mod_stress_invariants,  only: calc_sig_inv, calc_J2_inv, calc_J3_inv
   use mod_stress_invar_deriv, only: calc_dp_by_dsig, calc_dq_by_dsig, &
                                     calc_dJ3_by_dsig, calc_dlode_angle_by_dsig
   use mod_voigt_utils,        only: calc_dev_stress
   use mod_strain_invariants,  only: calc_eps_vol_inv, calc_dev_strain, calc_eps_q_inv
   use mod_state_params,       only: Get_I_coeff, check4crossing, Update_GK, &
                                     get_dilation, Get_M
   use mod_state_params_deriv, only: calc_ddil_by_deps_p
   implicit none
   private

   public :: mcsr_params_t, mcsr_state_t
   public :: mcsr_yield_fn
   public :: mcsr_flow_rule
   public :: mcsr_plastic_potential
   public :: mcsr_update_hardening
   public :: mcsr_elastic_stiffness
   public :: mcsr_hardening_modulus
   public :: mcsr_update_rate_state

   real(wp), parameter :: PI = 4.0_wp * atan(1.0_wp)

contains

   ! ---------------------------------------------------------------------------
   ! Yield function: F = q + eta_y * p
   ! ---------------------------------------------------------------------------
   pure function mcsr_yield_fn(params, state, sig) result(F)
      type(mcsr_params_t), intent(in) :: params
      type(mcsr_state_t),  intent(in) :: state
      real(wp),            intent(in) :: sig(6)
      real(wp) :: F

      real(wp) :: p, q, lode

      call calc_sig_inv(sig, p, q, lode)
      F = q + state%eta_y * p
   end function mcsr_yield_fn

   ! ---------------------------------------------------------------------------
   ! Flow rule: dF/dsig = eta_y * dp/dsig + dq/dsig + dF/dtheta * dtheta/dsig
   ! ---------------------------------------------------------------------------
   pure function mcsr_flow_rule(params, state, sig) result(dF_by_dsig)
      type(mcsr_params_t), intent(in) :: params
      type(mcsr_state_t),  intent(in) :: state
      real(wp),            intent(in) :: sig(6)
      real(wp) :: dF_by_dsig(6)

      real(wp) :: p, q, lode, J2, J3
      real(wp) :: dev(6)
      real(wp) :: dp_dsig(6), dq_dsig(6), dJ3_dsig(6), dlode_dsig(6)
      real(wp) :: dF_by_dtheta

      call calc_sig_inv(sig, p, q, lode)
      dev      = calc_dev_stress(sig, p)
      J2       = calc_J2_inv(dev)
      J3       = calc_J3_inv(dev)

      dp_dsig   = calc_dp_by_dsig()
      dq_dsig   = calc_dq_by_dsig(dev, q)
      dJ3_dsig  = calc_dJ3_by_dsig(dev)
      dlode_dsig = calc_dlode_angle_by_dsig(dJ3_dsig, dev, J3, J2, lode)

      ! dF/dtheta from Lode-dependent M: M = M_tc*(1 + 0.25*cos(1.5*theta + 0.25*pi)^1.2)
      dF_by_dtheta = 0.45_wp * p * params%M_tc &
                     * cos(1.5_wp * lode + 0.25_wp * PI)**0.2_wp &
                     * sin(1.5_wp * lode + 0.25_wp * PI)

      dF_by_dsig = state%eta_y * dp_dsig + dq_dsig + dF_by_dtheta * dlode_dsig
   end function mcsr_flow_rule

   ! ---------------------------------------------------------------------------
   ! Plastic potential: dG/dsig = -dilation * dp/dsig + dq/dsig
   ! ---------------------------------------------------------------------------
   pure function mcsr_plastic_potential(params, state, sig) result(dG_by_dsig)
      type(mcsr_params_t), intent(in) :: params
      type(mcsr_state_t),  intent(in) :: state
      real(wp),            intent(in) :: sig(6)
      real(wp) :: dG_by_dsig(6)

      real(wp) :: p, q, lode, dev(6)

      call calc_sig_inv(sig, p, q, lode)
      dev        = calc_dev_stress(sig, p)
      dG_by_dsig = -state%dilation * calc_dp_by_dsig() + calc_dq_by_dsig(dev, q)
   end function mcsr_plastic_potential

   ! ---------------------------------------------------------------------------
   ! Hardening update: accumulate eps_p, update dilation, update eta_y.
   !
   ! Only dilation evolves inside the integrator. eta_y depends on dilation and
   ! M_lode, and M_lode is frozen for the timestep (computed pre-step).
   ! ---------------------------------------------------------------------------
   subroutine mcsr_update_hardening(params, state, deps_p)
      type(mcsr_params_t), intent(in)    :: params
      type(mcsr_state_t),  intent(inout) :: state
      real(wp),            intent(in)    :: deps_p(6)

      real(wp) :: eps_vol, dev(6), eps_q

      state%eps_p = state%eps_p + deps_p

      ! Accumulated deviatoric plastic strain
      eps_vol = calc_eps_vol_inv(state%eps_p)
      dev     = calc_dev_strain(state%eps_p, eps_vol)
      eps_q   = calc_eps_q_inv(dev)

      ! Update dilation with plastic strain hardening/softening
      call get_dilation(params%h, params%D_min, state%I_coeff, params%ref_e_rate, &
                        eps_q, params%alpha_D, .true., state%dilation)

      ! Update friction ratio from current dilation and M_lode
      state%eta_y = state%M_lode - state%dilation * (1.0_wp - params%N)
   end subroutine mcsr_update_hardening

   ! ---------------------------------------------------------------------------
   ! Elastic stiffness: isotropic from current rate-updated G and K
   ! ---------------------------------------------------------------------------
   pure function mcsr_elastic_stiffness(params, state) result(stiff_e)
      type(mcsr_params_t), intent(in) :: params
      type(mcsr_state_t),  intent(in) :: state
      real(wp) :: stiff_e(6,6)

      stiff_e = calc_stiffness_GK(state%G, state%K)
   end function mcsr_elastic_stiffness

   ! ---------------------------------------------------------------------------
   ! Hardening modulus: H = -p * (1-N) * dot(a, dg/dsig)
   !
   ! Derived from consistency condition via chain rule:
   !   dF/deta_y * deta_y/ddil * ddil/deps_p dot dg/dsig
   ! where dF/deta_y = p,  deta_y/ddil = -(1-N),  ddil/deps_p = a (6-vector).
   ! H > 0 for softening (dilation decreasing beyond peak).
   ! ---------------------------------------------------------------------------
   subroutine mcsr_hardening_modulus(params, state, sig, dg_by_dsig, H)
      type(mcsr_params_t), intent(in)  :: params
      type(mcsr_state_t),  intent(in)  :: state
      real(wp),            intent(in)  :: sig(6), dg_by_dsig(6)
      real(wp),            intent(out) :: H

      real(wp) :: p, q, lode
      real(wp) :: a(6)
      real(wp) :: eps_vol, eps_q

      call calc_sig_inv(sig, p, q, lode)

      eps_vol = calc_eps_vol_inv(state%eps_p)
      eps_q   = calc_eps_q_inv(calc_dev_strain(state%eps_p, eps_vol))

      call calc_ddil_by_deps_p(params%D_min, params%h, params%ref_e_rate, params%alpha_D, &
                                eps_q, eps_vol, state%eps_p, state%I_coeff, .true., a)

      ! H = dF/deta_y * deta_y/ddil * ddil/deps_p . dg/dsig
      !   = p * (-(1-N)) * a . dg/dsig
      H = -p * (1.0_wp - params%N) * dot_product(a, dg_by_dsig)
   end subroutine mcsr_hardening_modulus

   ! ---------------------------------------------------------------------------
   ! Pre-step rate state update (called once per timestep in UMAT wrapper,
   ! before calling integrate_stress).
   !
   ! Steps:
   !   1. Smooth strain rate (if switch_smooth)
   !   2. check4crossing: determine if rate effects apply
   !   3. Compute inertial coefficient I
   !   4. Update G, K from rate (Update_GK)
   !   5. Compute M_lode from M_tc and Lode angle of initial stress
   !   6. Compute dilation at current rate and plastic strain
   !   7. Update eta_y = M_lode - dilation*(1-N)
   ! ---------------------------------------------------------------------------
   subroutine mcsr_update_rate_state(params, state, e_rate_scalar, sig6)
      type(mcsr_params_t), intent(in)    :: params
      type(mcsr_state_t),  intent(inout) :: state
      real(wp),            intent(in)    :: e_rate_scalar  !! Scalar strain rate magnitude [1/s]
      real(wp),            intent(in)    :: sig6(6)        !! Current stress (6-vector, UMAT convention)

      real(wp) :: p, q, lode
      real(wp) :: e_rate_smoothed, e_rate_prev
      real(wp) :: d_e_rate, I_new, M_lode
      real(wp) :: eps_vol, eps_q
      logical  :: apply_rate

      ! --- 1. Strain rate smoothing ---
      if (params%switch_smooth) then
         state%N_i     = state%N_i + 1
         state%SUM_rate = state%SUM_rate + e_rate_scalar
         if (state%N_i >= params%N_S) then
            e_rate_smoothed = state%SUM_rate / real(state%N_i, wp)
            state%N_i       = 0
            state%SUM_rate  = 0.0_wp
         else
            e_rate_smoothed = e_rate_scalar
         end if
      else
         e_rate_smoothed = e_rate_scalar
      end if

      ! --- 2. check4crossing: detect rate-reference crossing ---
      e_rate_prev = state%I_coeff  ! temporarily reuse I as prev rate placeholder
      ! Note: check4crossing modifies its first two arguments in-place
      call check4crossing(e_rate_prev, e_rate_smoothed, d_e_rate, params%ref_e_rate, apply_rate)

      ! --- 3. Compute inertial coefficient I ---
      call calc_sig_inv(sig6, p, q, lode)
      call Get_I_coeff(params%D_part, params%G_s, p, e_rate_smoothed, I_new)
      state%I_coeff = I_new

      ! --- 4. Update elastic moduli ---
      call Update_GK(params%G_0, params%nu, I_new, params%ref_e_rate, &
                     params%alpha_G, params%alpha_K, state%G, state%K)

      ! --- 5. Compute Lode-corrected M ---
      call Get_M(params%M_tc, lode, M_lode)
      state%M_lode = M_lode

      ! --- 6. Update dilation at current rate and plastic strain ---
      eps_vol = calc_eps_vol_inv(state%eps_p)
      eps_q   = calc_eps_q_inv(calc_dev_strain(state%eps_p, eps_vol))

      call get_dilation(params%h, params%D_min, I_new, params%ref_e_rate, &
                        eps_q, params%alpha_D, apply_rate, state%dilation)

      ! --- 7. Update friction ratio ---
      state%eta_y = M_lode - state%dilation * (1.0_wp - params%N)

   end subroutine mcsr_update_rate_state

end module mod_mcsr_functions
