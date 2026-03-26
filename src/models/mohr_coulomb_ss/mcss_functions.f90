! Pure constitutive functions for the Mohr-Coulomb Strain Softening model.
!
! Yield surface and plastic potential follow the smooth hyperbolic approximation
! of Abbo & Sloan (1995). All procedures are pure and operate on explicit
! (params, state, ...) arguments — no class(), no dynamic dispatch.
!
! Smoothing constants (LodeT = 29.5 degrees):
!   A1, A2, B1, B2, LODE_T  — polynomial rounding at yield surface corners
!   SMOOTH_COEFF             — tip smoothing: a = SMOOTH_COEFF * c * cot(phi)
module mod_mcss_functions
   use mod_csm_kinds,  only: wp
   use mod_mcss_types, only: mcss_params_t, mcss_state_t
   use mod_elastic_utils, only: calc_stiffness_GK, calc_K_from_G_nu
   implicit none
   private

   public :: mcss_params_t, mcss_state_t
   public :: mcss_yield_fn
   public :: mcss_flow_rule
   public :: mcss_plastic_potential
   public :: mcss_update_hardening
   public :: mcss_elastic_stiffness

   ! ---------------------------------------------------------------------------
   ! Abbo & Sloan (1995) smoothing constants — LodeT = 29.5 degrees
   ! ---------------------------------------------------------------------------
   real(wp), parameter :: LODE_T      = 0.514872129338327_wp  !! Transition Lode angle [rad]
   real(wp), parameter :: A1          = 7.138654723242414_wp
   real(wp), parameter :: A2          = 6.112267270920612_wp
   real(wp), parameter :: B1          = 6.270447753139589_wp
   real(wp), parameter :: B2          = 6.398760841429403_wp
   real(wp), parameter :: SMOOTH_COEFF = 0.0005_wp             !! Tip smoothing: a = SMOOTH_COEFF*c*cot(phi)
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

      real(wp) :: p, J, lode, s3ta, K, a, asphi2

      call calc_invariants(sig, p, J, lode, s3ta)
      K      = calc_K(lode, s3ta, sin(state%phi))
      a      = calc_a_smooth(state%c, state%phi)
      asphi2 = a * a * sin(state%phi)**2

      F = p * sin(state%phi) + sqrt(J*J * K*K + asphi2) - state%c * cos(state%phi)
   end function mcss_yield_fn

   ! ---------------------------------------------------------------------------
   ! Flow rule: dF/dsig (normal to yield surface)
   ! ---------------------------------------------------------------------------
   pure function mcss_flow_rule(params, state, sig) result(dF_by_dsig)
      type(mcss_params_t), intent(in) :: params
      type(mcss_state_t),  intent(in) :: state
      real(wp),            intent(in) :: sig(6)
      real(wp) :: dF_by_dsig(6)

      real(wp) :: p, J, lode, s3ta
      call calc_invariants(sig, p, J, lode, s3ta)
      dF_by_dsig = calc_dF_dsig_abbo(sig, p, J, lode, s3ta, state%c, state%phi)
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

      real(wp) :: p, J, lode, s3ta, psi_eff

      call calc_invariants(sig, p, J, lode, s3ta)

      ! Interpolate psi -> phi near J=0 to remove singularity in derivatives
      if (J < J0) then
         psi_eff = state%phi - J * (state%phi - state%psi) / J0
      else
         psi_eff = state%psi
      end if

      if (psi_eff == state%phi) then
         dG_by_dsig = calc_dF_dsig_abbo(sig, p, J, lode, s3ta, state%c, state%phi)
      else
         dG_by_dsig = calc_dF_dsig_abbo(sig, p, J, lode, s3ta, state%c, psi_eff)
      end if
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
      eps_p_eq    = calc_eps_p_eq(state%eps_p)

      state%c   = params%c_res   + (params%c_peak   - params%c_res)   * exp(-params%factor * eps_p_eq)
      state%phi = params%phi_res + (params%phi_peak  - params%phi_res) * exp(-params%factor * eps_p_eq)
      state%psi = params%psi_res + (params%psi_peak  - params%psi_res) * exp(-params%factor * eps_p_eq)
   end subroutine mcss_update_hardening

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
   ! Private helpers
   ! ===========================================================================

   ! Stress invariants (p, J=sqrt(J2), Lode angle, sin(3*theta))
   pure subroutine calc_invariants(sig, p, J, lode, s3ta)
      real(wp), intent(in)  :: sig(6)
      real(wp), intent(out) :: p, J, lode, s3ta

      real(wp) :: Sx, Sy, Sz, J2, J3, h1, h2

      p  = ONE_THIRD * (sig(1) + sig(2) + sig(3))
      Sx = sig(1) - p
      Sy = sig(2) - p
      Sz = sig(3) - p

      J2 = (1.0_wp/6.0_wp) * ((sig(1)-sig(2))**2 + (sig(1)-sig(3))**2 + (sig(2)-sig(3))**2) &
           + sig(4)**2 + sig(5)**2 + sig(6)**2
      J  = sqrt(J2)

      J3 = Sx*Sy*Sz + 2.0_wp*sig(4)*sig(5)*sig(6) &
           - Sx*sig(5)**2 - Sy*sig(6)**2 - Sz*sig(4)**2

      if (J2 > 0.0_wp) then
         h1   = -3.0_wp / (2.0_wp * INV_SQRT3)
         h2   = J3 / J**3
         s3ta = h1 * h2
         s3ta = max(-1.0_wp, min(1.0_wp, s3ta))   ! clamp for asin
         lode = ONE_THIRD * asin(s3ta)
      else
         lode = 0.0_wp
         s3ta = 0.0_wp
      end if
   end subroutine calc_invariants

   ! Lode-angle shape function K(theta) — Abbo & Sloan (1995)
   pure function calc_K(lode, s3ta, sphi) result(K)
      real(wp), intent(in) :: lode, s3ta, sphi
      real(wp) :: K, A, B, sgn

      if (abs(lode) < LODE_T) then
         K = cos(lode) - INV_SQRT3 * sphi * sin(lode)
      else
         sgn = sign(1.0_wp, lode)
         A   = A1 + A2 * sgn * sphi
         B   = B1 * sgn + B2 * sphi
         K   = A - B * s3ta
      end if
   end function calc_K

   ! dK/dLode
   pure function calc_dK_dlode(lode, s3ta, sphi, cphi) result(dKdL)
      real(wp), intent(in) :: lode, s3ta, sphi, cphi
      real(wp) :: dKdL, cta, sta, c3ta, sgn, B

      if (abs(lode) < LODE_T) then
         cta  = cos(lode)
         sta  = s3ta / (4.0_wp*cta*cta - 1.0_wp)
         dKdL = -sta - INV_SQRT3 * sphi * cta
      else
         sgn  = sign(1.0_wp, lode)
         cta  = cos(lode)
         c3ta = cta * (4.0_wp*cta*cta - 3.0_wp)
         B    = B1 * sgn + B2 * sphi
         dKdL = -3.0_wp * B * c3ta
      end if
   end function calc_dK_dlode

   ! Tip smoothing parameter a = SMOOTH_COEFF * c * cot(phi)
   pure function calc_a_smooth(c, phi) result(a)
      real(wp), intent(in) :: c, phi
      real(wp) :: a
      if (abs(phi) < 1.0e-14_wp) then
         a = SMOOTH_COEFF * c   ! cot(phi) -> large; use c-only fallback
      else
         a = SMOOTH_COEFF * c * cos(phi) / sin(phi)
      end if
   end function calc_a_smooth

   ! Equivalent plastic strain: eps_p_eq = sqrt((2/3) * e_p:e_p)
   ! where e_p is the deviatoric plastic strain.
   ! Shear components carry factor 4/3 = 2*(2/3) due to engineering convention.
   pure function calc_eps_p_eq(eps_p) result(eps_p_eq)
      real(wp), intent(in) :: eps_p(6)
      real(wp) :: eps_p_eq

      real(wp) :: eps_pm, e1, e2, e3
      real(wp), parameter :: TWO_THIRDS  = 2.0_wp / 3.0_wp
      real(wp), parameter :: FOUR_THIRDS = 4.0_wp / 3.0_wp

      eps_pm = ONE_THIRD * (eps_p(1) + eps_p(2) + eps_p(3))
      e1 = eps_p(1) - eps_pm
      e2 = eps_p(2) - eps_pm
      e3 = eps_p(3) - eps_pm

      eps_p_eq = sqrt(TWO_THIRDS  * (e1*e1 + e2*e2 + e3*e3) + &
                      FOUR_THIRDS * (eps_p(4)**2 + eps_p(5)**2 + eps_p(6)**2))
   end function calc_eps_p_eq

   ! dF/dsig for the Abbo & Sloan yield/potential surface with a given angle `angle`
   ! (pass phi for yield function, psi_eff for plastic potential)
   pure function calc_dF_dsig_abbo(sig, p, J, lode, s3ta, c, angle) result(dFds)
      real(wp), intent(in) :: sig(6), p, J, lode, s3ta, c, angle
      real(wp) :: dFds(6)

      real(wp) :: sphi, cphi, K, dKdL, a, asphi2, D
      real(wp) :: C1, C2, C3, T3TA, C3TA, J2
      real(wp) :: Sx, Sy, Sz
      real(wp) :: dp_ds(6), dJ_ds(6), dJ3_ds(6)
      real(wp) :: i1, i2

      sphi  = sin(angle)
      cphi  = cos(angle)
      K     = calc_K(lode, s3ta, sphi)
      dKdL  = calc_dK_dlode(lode, s3ta, sphi, cphi)
      a     = calc_a_smooth(c, angle)
      asphi2 = a * a * sphi * sphi

      J2 = max(J * J, J_TINY)

      ! cos(3*theta) and tan(3*theta)
      C3TA = cos(lode) * (4.0_wp*cos(lode)**2 - 3.0_wp)
      if (abs(C3TA) < 1.0e-14_wp) C3TA = sign(1.0e-14_wp, C3TA)
      T3TA = s3ta / C3TA

      ! dp/dsig
      dp_ds = [ONE_THIRD, ONE_THIRD, ONE_THIRD, 0.0_wp, 0.0_wp, 0.0_wp]

      ! dJ/dsig
      Sx = sig(1) - p
      Sy = sig(2) - p
      Sz = sig(3) - p
      if (J > 1.0e-4_wp) then
         i1 = 0.5_wp / J
      else
         i1 = 0.0_wp
      end if
      dJ_ds = [i1*Sx, i1*Sy, i1*Sz, i1*2.0_wp*sig(4), i1*2.0_wp*sig(5), i1*2.0_wp*sig(6)]

      ! dJ3/dsig
      i2 = ONE_THIRD * J2
      dJ3_ds(1) = Sy*Sz - sig(5)**2 + i2
      dJ3_ds(2) = Sx*Sz - sig(6)**2 + i2
      dJ3_ds(3) = Sx*Sy - sig(4)**2 + i2
      dJ3_ds(4) = 2.0_wp*(sig(5)*sig(6) - Sz*sig(4))
      dJ3_ds(5) = 2.0_wp*(sig(6)*sig(4) - Sx*sig(5))
      dJ3_ds(6) = 2.0_wp*(sig(4)*sig(5) - Sy*sig(6))

      ! Scalar coefficients
      D  = J * K / sqrt(J2*K*K + asphi2)
      C1 = sphi
      C2 = D*K - T3TA * D * dKdL
      C3 = -SQRT3_OVER2 * dKdL * D / (J2 * C3TA)

      dFds = C1*dp_ds + C2*dJ_ds + C3*dJ3_ds
   end function calc_dF_dsig_abbo

end module mod_mcss_functions
