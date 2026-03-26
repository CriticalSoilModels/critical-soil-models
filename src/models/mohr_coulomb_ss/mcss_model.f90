! Mohr-Coulomb Strain Softening — concrete implementation of csm_model_t.
!
! All constitutive math lives in mod_mcss_functions (pure, GPU-callable).
! The deferred procedure implementations here are thin wrappers only.
!
! PROPS layout:
!   Required (1–13):
!    1:  G          shear modulus [kPa]
!    2:  nu         Poisson's ratio [-]
!    3:  c_peak     peak cohesion [kPa]
!    4:  c_res      residual cohesion [kPa]
!    5:  phi_peak   peak friction angle [degrees]
!    6:  phi_res    residual friction angle [degrees]
!    7:  psi_peak   peak dilation angle [degrees]
!    8:  psi_res    residual dilation angle [degrees]
!    9:  factor     softening rate shape parameter [-]
!   10:  (integrator flag — consumed by umat_mcss, not stored here)
!   11:  yield_tol
!   12:  max_iters
!   13:  dt_min
!
!   Optional Abbo & Sloan smoothing (14–19):
!   If omitted or zero, DEFAULT_AS_PARAMS values are used (LodeT = 29.5 deg).
!   14:  lode_t       transition Lode angle [rad]
!   15:  A1           corner-rounding constant
!   16:  A2           corner-rounding constant
!   17:  B1           corner-rounding constant
!   18:  B2           corner-rounding constant
!   19:  smooth_coeff tip smoothing coefficient [-]
!
! STATEV layout (9 entries):
!   1:   c      current cohesion [kPa]
!   2:   phi    current friction angle [rad]
!   3:   psi    current dilation angle [rad]
!   4-9: eps_p  accumulated plastic strain (Voigt) [-]

module mod_mcss_model
   use mod_csm_kinds,      only: wp
   use mod_csm_model,      only: csm_model_t
   use mod_mcss_types,     only: mcss_params_t, mcss_state_t, abbo_sloan_params_t, DEFAULT_AS_PARAMS
   use mod_mcss_functions, only: mcss_yield_fn, mcss_flow_rule, mcss_plastic_potential, &
                                  mcss_update_hardening, mcss_elastic_stiffness
   implicit none
   private

   public :: mcss_model_t, mcss_from_props, mcss_load_state, mcss_save_state

   real(wp), parameter :: DEG_TO_RAD = acos(-1.0_wp) / 180.0_wp

   type, extends(csm_model_t) :: mcss_model_t
      type(mcss_params_t) :: params
      type(mcss_state_t)  :: state
      real(wp) :: yield_tol
      integer  :: max_iters
      real(wp) :: dt_min
   contains
      procedure :: yield_fn          => mcss_yield_fn_method
      procedure :: flow_rule         => mcss_flow_rule_method
      procedure :: plastic_potential => mcss_plastic_potential_method
      procedure :: update_hardening  => mcss_update_hardening_method
      procedure :: elastic_stiffness => mcss_elastic_stiffness_method
      procedure :: snapshot          => mcss_snapshot
      procedure :: restore           => mcss_restore
   end type mcss_model_t

contains

   ! ---------------------------------------------------------------------------
   ! UMAT boundary helpers — only place PROPS/STATEV index magic lives
   ! ---------------------------------------------------------------------------

   function mcss_from_props(props) result(model)
      real(wp), intent(in) :: props(:)
      type(mcss_model_t) :: model

      model%params%G        = props(1)
      model%params%nu       = props(2)
      model%params%c_peak   = props(3)
      model%params%c_res    = props(4)
      model%params%phi_peak = props(5) * DEG_TO_RAD
      model%params%phi_res  = props(6) * DEG_TO_RAD
      model%params%psi_peak = props(7) * DEG_TO_RAD
      model%params%psi_res  = props(8) * DEG_TO_RAD
      model%params%factor   = props(9)
      model%yield_tol       = props(11)
      model%max_iters       = int(props(12))
      model%dt_min          = props(13)

      model%params%as_params = mcss_as_params_from_props(props)
   end function mcss_from_props

   !> Build abbo_sloan_params_t from PROPS(14:19).
   !> Any entry that is missing (array too short) or zero uses the default value.
   pure function mcss_as_params_from_props(props) result(asp)
      real(wp), intent(in) :: props(:)
      type(abbo_sloan_params_t) :: asp

      asp = DEFAULT_AS_PARAMS   ! start from defaults

      if (size(props) >= 14) then
         if (props(14) /= 0.0_wp) asp%lode_t       = props(14)
      end if
      if (size(props) >= 15) then
         if (props(15) /= 0.0_wp) asp%A1           = props(15)
      end if
      if (size(props) >= 16) then
         if (props(16) /= 0.0_wp) asp%A2           = props(16)
      end if
      if (size(props) >= 17) then
         if (props(17) /= 0.0_wp) asp%B1           = props(17)
      end if
      if (size(props) >= 18) then
         if (props(18) /= 0.0_wp) asp%B2           = props(18)
      end if
      if (size(props) >= 19) then
         if (props(19) /= 0.0_wp) asp%smooth_coeff = props(19)
      end if
   end function mcss_as_params_from_props

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
   ! Thin wrappers around pure functions
   ! ---------------------------------------------------------------------------

   function mcss_yield_fn_method(self, sig) result(F)
      class(mcss_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: F
      F = mcss_yield_fn(self%params, self%state, sig)
   end function mcss_yield_fn_method

   function mcss_flow_rule_method(self, sig) result(dF_by_dsig)
      class(mcss_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: dF_by_dsig(6)
      dF_by_dsig = mcss_flow_rule(self%params, self%state, sig)
   end function mcss_flow_rule_method

   function mcss_plastic_potential_method(self, sig) result(dG_by_dsig)
      class(mcss_model_t), intent(in) :: self
      real(wp), intent(in) :: sig(6)
      real(wp) :: dG_by_dsig(6)
      dG_by_dsig = mcss_plastic_potential(self%params, self%state, sig)
   end function mcss_plastic_potential_method

   subroutine mcss_update_hardening_method(self, deps_p)
      class(mcss_model_t), intent(inout) :: self
      real(wp), intent(in) :: deps_p(6)
      call mcss_update_hardening(self%params, self%state, deps_p)
   end subroutine mcss_update_hardening_method

   function mcss_elastic_stiffness_method(self) result(stiff_e)
      class(mcss_model_t), intent(in) :: self
      real(wp) :: stiff_e(6,6)
      stiff_e = mcss_elastic_stiffness(self%params, self%state)
   end function mcss_elastic_stiffness_method

   ! ---------------------------------------------------------------------------
   ! snapshot / restore — for integrator rollback
   ! ---------------------------------------------------------------------------

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
