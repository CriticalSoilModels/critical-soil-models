! Concrete implementation of csm_model_t for Mohr-Coulomb Strain Softening.
! Pattern every model follows:
!   1. Extend csm_model_t with named param + state fields
!   2. Implement the seven deferred procedures
!   3. Provide from_props / load_state / save_state at the UMAT boundary

module mod_mcss_model
   use mod_csm_kinds,         only: wp
   use mod_csm_model,         only: csm_model_t
   use mod_elastic_utils,     only: calc_stiffness_GK, calc_K_from_G_nu
   use mod_stress_invariants, only: calc_sig_inv
   implicit none

   real(wp), parameter :: DEG_TO_RAD = acos(-1.0_wp) / 180.0_wp

   ! ---------------------------------------------------------------------------
   ! State type — separate so it can be copied cleanly for integrator rollback
   ! ---------------------------------------------------------------------------
   ! All index magic for STATEV lives only in load_state / save_state / snapshot / restore.
   ! Everywhere else state is accessed by name.

   type :: mcss_state_t
      real(wp) :: c, phi, psi     ! current softening values
      real(wp) :: eps_p(6)        ! accumulated plastic strain
   end type mcss_state_t

   ! ---------------------------------------------------------------------------
   ! Concrete model type
   ! ---------------------------------------------------------------------------

   type, extends(csm_model_t) :: mcss_model_t

      ! --- Parameters (fixed, populated from PROPS once at UMAT entry) ---
      real(wp) :: G, nu
      real(wp) :: c_peak,   c_res
      real(wp) :: phi_peak, phi_res     ! stored in radians
      real(wp) :: psi_peak, psi_res     ! stored in radians
      real(wp) :: shape_factor
      real(wp) :: yield_tol
      integer  :: max_iters
      real(wp) :: dt_min

      ! --- State (evolves during integration) ---
      type(mcss_state_t) :: state

   contains
      procedure :: yield_fn          => mcss_yield_fn
      procedure :: flow_rule         => mcss_flow_rule
      procedure :: plastic_potential => mcss_plastic_potential
      procedure :: update_hardening  => mcss_update_hardening
      procedure :: elastic_stiffness => mcss_elastic_stiffness
      procedure :: snapshot          => mcss_snapshot
      procedure :: restore           => mcss_restore
   end type mcss_model_t

contains

   ! ---------------------------------------------------------------------------
   ! UMAT boundary helpers — only place PROPS/STATEV index magic lives
   ! ---------------------------------------------------------------------------

   ! PROPS layout for MCSS:
   !   1: G          shear modulus
   !   2: nu         Poisson's ratio
   !   3: c_peak     peak cohesion
   !   4: c_res      residual cohesion
   !   5: phi_peak   peak friction angle (degrees)
   !   6: phi_res    residual friction angle (degrees)
   !   7: psi_peak   peak dilation angle (degrees)
   !   8: psi_res    residual dilation angle (degrees)
   !   9: shape_factor
   !  10: integration_flag   (0 = Euler, 1 = Ortiz-Simo)
   !  11: yield_tol
   !  12: max_iters
   !  13: dt_min

   function mcss_from_props(props) result(model)
      real(wp), intent(in) :: props(:)
      type(mcss_model_t) :: model

      model%G            = props(1)
      model%nu           = props(2)
      model%c_peak       = props(3)
      model%c_res        = props(4)
      model%phi_peak     = props(5) * DEG_TO_RAD
      model%phi_res      = props(6) * DEG_TO_RAD
      model%psi_peak     = props(7) * DEG_TO_RAD
      model%psi_res      = props(8) * DEG_TO_RAD
      model%shape_factor = props(9)
      model%yield_tol    = props(11)
      model%max_iters    = int(props(12))
      model%dt_min       = props(13)
   end function mcss_from_props

   ! STATEV layout for MCSS:
   !   1:   c      current cohesion
   !   2:   phi    current friction angle (radians)
   !   3:   psi    current dilation angle (radians)
   !   4-9: eps_p  accumulated plastic strain (Voigt, internal convention)
   !
   ! This is the ONLY place STATEV indices appear for this model.

   subroutine mcss_load_state(model, statev)
      type(mcss_model_t), intent(inout) :: model
      real(wp),           intent(in)    :: statev(:)

      model%state%c     = statev(1)
      model%state%phi   = statev(2)
      model%state%psi   = statev(3)
      model%state%eps_p = statev(4:9)
   end subroutine mcss_load_state

   subroutine mcss_save_state(model, statev)
      type(mcss_model_t), intent(in)    :: model
      real(wp),           intent(inout) :: statev(:)

      statev(1)   = model%state%c
      statev(2)   = model%state%phi
      statev(3)   = model%state%psi
      statev(4:9) = model%state%eps_p
   end subroutine mcss_save_state

   ! ---------------------------------------------------------------------------
   ! Implementations of the five deferred procedures
   ! ---------------------------------------------------------------------------

   function mcss_yield_fn(self, sig) result(F)
      class(mcss_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: F
      real(wp) :: p, q, lode_angle

      call calc_sig_inv(sig, p, q, lode_angle)
      ! Mohr-Coulomb: F = q - eta*p - c
      ! where eta = M(theta) derived from phi
      ! (simplified here — real version uses Lode-angle-dependent M)
      F = q - self%state%phi * p - self%state%c
   end function mcss_yield_fn

   function mcss_flow_rule(self, sig) result(dF_by_dsig)
      class(mcss_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: dF_by_dsig(6)

      ! ∂F/∂σ — gradient of yield surface w.r.t. stress
      ! Uses self%state%phi (current friction angle)
      ! ... invariant chain rule ...
      dF_by_dsig = 0.0_wp   ! placeholder
   end function mcss_flow_rule

   function mcss_plastic_potential(self, sig) result(dG_by_dsig)
      class(mcss_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: dG_by_dsig(6)

      ! ∂G/∂σ — flow direction
      ! Non-associated: same form as flow_rule but uses self%state%psi not phi
      ! ... invariant chain rule ...
      dG_by_dsig = 0.0_wp   ! placeholder
   end function mcss_plastic_potential

   subroutine mcss_update_hardening(self, deps_p)
      class(mcss_model_t), intent(inout) :: self
      real(wp), intent(in) :: deps_p(6)
      real(wp) :: eps_p_eq

      self%state%eps_p = self%state%eps_p + deps_p

      ! Equivalent plastic strain (scalar measure of accumulated plasticity)
      eps_p_eq = sqrt(2.0_wp/3.0_wp * dot_product(self%state%eps_p, self%state%eps_p))

      ! Softening: interpolate c, phi, psi between peak and residual values
      ! call calc_softening_params(eps_p_eq, self%c_peak, self%c_res, &
      !    self%phi_peak, self%phi_res, self%psi_peak, self%psi_res,  &
      !    self%shape_factor, self%state%c, self%state%phi, self%state%psi)
      ! TODO: implement mcss_functions.f90 and use it here
   end subroutine mcss_update_hardening

   function mcss_elastic_stiffness(self) result(stiff_e)
      class(mcss_model_t), intent(in) :: self
      real(wp) :: stiff_e(6,6)

      stiff_e = calc_stiffness_GK(self%G, calc_K_from_G_nu(self%G, self%nu))

   end function mcss_elastic_stiffness

   ! ---------------------------------------------------------------------------
   ! snapshot / restore — only place STATEV-style index layout appears internally
   ! The integrator calls these to roll back state during the substep loop.
   ! ---------------------------------------------------------------------------

   ! Serialise state to a flat array. Layout:
   !   1:   c
   !   2:   phi
   !   3:   psi
   !   4-9: eps_p

   subroutine mcss_snapshot(self, saved)
      class(mcss_model_t), intent(in)    :: self
      real(wp), allocatable, intent(out) :: saved(:)

      allocate(saved(9))
      saved(1)   = self%state%c
      saved(2)   = self%state%phi
      saved(3)   = self%state%psi
      saved(4:9) = self%state%eps_p
   end subroutine mcss_snapshot

   subroutine mcss_restore(self, saved)
      class(mcss_model_t), intent(inout) :: self
      real(wp),            intent(in)    :: saved(:)

      self%state%c     = saved(1)
      self%state%phi   = saved(2)
      self%state%psi   = saved(3)
      self%state%eps_p = saved(4:9)
   end subroutine mcss_restore

end module mod_mcss_model
