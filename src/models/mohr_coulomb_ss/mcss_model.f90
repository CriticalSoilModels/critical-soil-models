! =============================================================================
! PSEUDOCODE — does not compile
! Concrete implementation of csm_model_t for Mohr-Coulomb Strain Softening.
! Shows the pattern every model should follow:
!   1. Extend csm_model_t with named param + state fields
!   2. Implement the five deferred procedures
!   3. Provide from_props / load_state / save_state at the UMAT boundary
! =============================================================================

module mod_mcss_model
   use iso_fortran_env, only: dp => real64
   use mod_csm_model
   use mod_stress_invariants, only: get_invariants
   implicit none

   real(dp), parameter :: DEG_TO_RAD = acos(-1.0_dp) / 180.0_dp

   ! ---------------------------------------------------------------------------
   ! State type — separate so it can be copied cleanly for integrator rollback
   ! ---------------------------------------------------------------------------
   ! All index magic for STATEV lives only in load_state / save_state / snapshot / restore.
   ! Everywhere else state is accessed by name.

   type :: mcss_state_t
      real(dp) :: c, phi, psi     ! current softening values
      real(dp) :: eps_p(6)        ! accumulated plastic strain
   end type mcss_state_t

   ! ---------------------------------------------------------------------------
   ! Concrete model type
   ! ---------------------------------------------------------------------------

   type, extends(csm_model_t) :: mcss_model_t

      ! --- Parameters (fixed, populated from PROPS once at UMAT entry) ---
      real(dp) :: G, nu
      real(dp) :: c_peak,   c_res
      real(dp) :: phi_peak, phi_res     ! stored in radians
      real(dp) :: psi_peak, psi_res     ! stored in radians
      real(dp) :: shape_factor
      real(dp) :: yield_tol
      integer  :: max_iters
      real(dp) :: dt_min

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
      real(dp), intent(in) :: props(:)
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
      real(dp),           intent(in)    :: statev(:)

      model%state%c     = statev(1)
      model%state%phi   = statev(2)
      model%state%psi   = statev(3)
      model%state%eps_p = statev(4:9)
   end subroutine mcss_load_state

   subroutine mcss_save_state(model, statev)
      type(mcss_model_t), intent(in)    :: model
      real(dp),           intent(inout) :: statev(:)

      statev(1)   = model%state%c
      statev(2)   = model%state%phi
      statev(3)   = model%state%psi
      statev(4:9) = model%state%eps_p
   end subroutine mcss_save_state

contains

   ! ---------------------------------------------------------------------------
   ! Implementations of the five deferred procedures
   ! ---------------------------------------------------------------------------

   function mcss_yield_fn(self, sig) result(F)
      class(mcss_model_t), intent(in) :: self
      real(dp), intent(in) :: sig(6)
      real(dp) :: F
      real(dp) :: p, q, theta

      call get_invariants(sig, p, q, theta)
      ! Mohr-Coulomb: F = q - eta*p - c
      ! where eta = M(theta) derived from phi
      ! (simplified here — real version uses Lode-angle-dependent M)
      F = q - self%state%phi * p - self%state%c
   end function mcss_yield_fn

   function mcss_flow_rule(self, sig) result(dF_by_dsig)
      class(mcss_model_t), intent(in) :: self
      real(dp), intent(in) :: sig(6)
      real(dp) :: dF_by_dsig(6)

      ! ∂F/∂σ — gradient of yield surface w.r.t. stress
      ! Uses self%state%phi (current friction angle)
      ! ... invariant chain rule ...
      dF_by_dsig = 0.0_dp   ! placeholder
   end function mcss_flow_rule

   function mcss_plastic_potential(self, sig) result(dG_by_dsig)
      class(mcss_model_t), intent(in) :: self
      real(dp), intent(in) :: sig(6)
      real(dp) :: dG_by_dsig(6)

      ! ∂G/∂σ — flow direction
      ! Non-associated: same form as flow_rule but uses self%state%psi not phi
      ! ... invariant chain rule ...
      dG_by_dsig = 0.0_dp   ! placeholder
   end function mcss_plastic_potential

   subroutine mcss_update_hardening(self, deps_p)
      class(mcss_model_t), intent(inout) :: self
      real(dp), intent(in) :: deps_p(6)
      real(dp) :: eps_p_eq

      self%state%eps_p = self%state%eps_p + deps_p

      ! Equivalent plastic strain (scalar measure of accumulated plasticity)
      eps_p_eq = sqrt(2.0_dp/3.0_dp * dot_product(self%state%eps_p, self%state%eps_p))

      ! Softening: interpolate c, phi, psi between peak and residual values
      call calc_softening_params(eps_p_eq,            &
                                 self%c_peak,   self%c_res,   &
                                 self%phi_peak, self%phi_res, &
                                 self%psi_peak, self%psi_res, &
                                 self%shape_factor,           &
                                 self%state%c, self%state%phi, self%state%psi)
   end subroutine mcss_update_hardening

   function mcss_elastic_stiffness(self) result(stiff_e)
      class(mcss_model_t), intent(in) :: self
      real(dp) :: stiff_e(6,6)
      real(dp) :: lame_1, lame_2

      ! Params accessed flat — self%G, self%nu (no self%params% needed)
      lame_1 = 2*self%G*(1.0_dp - self%nu) / (1.0_dp - 2.0_dp*self%nu)
      lame_2 = 2*self%G*self%nu             / (1.0_dp - 2.0_dp*self%nu)

      stiff_e          = 0.0_dp
      stiff_e(1:3,1:3) = lame_2
      stiff_e(1,1) = lame_1;  stiff_e(2,2) = lame_1;  stiff_e(3,3) = lame_1
      stiff_e(4,4) = self%G;  stiff_e(5,5) = self%G;  stiff_e(6,6) = self%G
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
      real(dp), allocatable, intent(out) :: saved(:)

      allocate(saved(9))
      saved(1)   = self%state%c
      saved(2)   = self%state%phi
      saved(3)   = self%state%psi
      saved(4:9) = self%state%eps_p
   end subroutine mcss_snapshot

   subroutine mcss_restore(self, saved)
      class(mcss_model_t), intent(inout) :: self
      real(dp),            intent(in)    :: saved(:)

      self%state%c     = saved(1)
      self%state%phi   = saved(2)
      self%state%psi   = saved(3)
      self%state%eps_p = saved(4:9)
   end subroutine mcss_restore

end module mod_mcss_model
