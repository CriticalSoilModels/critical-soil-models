!! Integration tests for the MCSS model stack.
!!
!! These tests drive `mcss_model_t` through `euler_substep` and verify
!! physically consistent outputs. They complement the unit tests in
!! `test_mcss_functions_suite.f90`, which test individual pure functions
!! in isolation. Here we verify the full cycle:
!!
!!   elastic predictor → yield check → return mapping → hardening update
!!
!! Full strain-path validation (stress-strain curves vs. published data)
!! lives in the incremental-driver repo.

module mod_test_mcss_integration_suite
   use mod_csm_kinds,     only: wp
   use mod_mcss_model,    only: mcss_model_t, mcss_from_props, &
                                 mcss_load_state, mcss_save_state
   use mod_euler_substep, only: euler_substep, integrator_params_t
   use testdrive,         only: new_unittest, unittest_type, error_type, check
   implicit none
   private
   public :: collect_mcss_integration_suite

   real(wp), parameter :: PI         = 4.0_wp * atan(1.0_wp)
   real(wp), parameter :: DEG_TO_RAD = PI / 180.0_wp

   real(wp), parameter :: G_TEST       = 10000.0_wp  !! Shear modulus [kPa]
   real(wp), parameter :: NU_TEST      = 0.3_wp       !! Poisson's ratio [-]
   real(wp), parameter :: C_PEAK       = 10.0_wp      !! Peak cohesion [kPa]
   real(wp), parameter :: C_RES        = 2.0_wp       !! Residual cohesion [kPa]
   real(wp), parameter :: PHI_PEAK_DEG = 30.0_wp      !! Peak friction angle [degrees]
   real(wp), parameter :: PHI_RES_DEG  = 20.0_wp      !! Residual friction angle [degrees]
   real(wp), parameter :: PSI_PEAK_DEG = 10.0_wp      !! Peak dilation angle [degrees]
   real(wp), parameter :: PSI_RES_DEG  = 0.0_wp       !! Residual dilation angle [degrees]
   real(wp), parameter :: FACTOR_STD   = 100.0_wp     !! Softening rate [-]
   real(wp), parameter :: YIELD_TOL    = 1.0e-8_wp    !! Yield surface tolerance [-]

   real(wp), parameter :: PHI_PEAK = PHI_PEAK_DEG * DEG_TO_RAD  !! [rad]
   real(wp), parameter :: PHI_RES  = PHI_RES_DEG  * DEG_TO_RAD  !! [rad]
   real(wp), parameter :: PSI_PEAK = PSI_PEAK_DEG * DEG_TO_RAD  !! [rad]

contains

   subroutine collect_mcss_integration_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("elastic step: no plastic strain",         test_elastic_no_plastic_strain),  &
         new_unittest("plastic step: stress on yield surface",   test_plastic_returns_to_surface), &
         new_unittest("large strain: approaches residual state", test_softening_to_residual)       &
      ]
   end subroutine collect_mcss_integration_suite

   ! ---------------------------------------------------------------------------
   ! Helper: build model from PROPS + STATEV at peak state
   ! ---------------------------------------------------------------------------

   subroutine make_model(model, factor)
      !! Construct an mcss_model_t at peak strength using the standard test parameters.
      type(mcss_model_t), intent(out) :: model
      real(wp),           intent(in)  :: factor  !! Softening rate parameter
      real(wp) :: props(13), statev(9)

      props = [G_TEST, NU_TEST, C_PEAK, C_RES, &
               PHI_PEAK_DEG, PHI_RES_DEG, PSI_PEAK_DEG, PSI_RES_DEG, &
               factor, 0.0_wp, YIELD_TOL, 100.0_wp, 1.0e-9_wp]

      statev      = 0.0_wp
      statev(1)   = C_PEAK   ! c   [kPa]
      statev(2)   = PHI_PEAK ! phi [rad] — mcss_load_state expects radians
      statev(3)   = PSI_PEAK ! psi [rad]

      model = mcss_from_props(props)
      call mcss_load_state(model, statev)
   end subroutine make_model

   ! ---------------------------------------------------------------------------
   ! Test 1: elastic step — integrator skips return mapping
   ! ---------------------------------------------------------------------------

   subroutine test_elastic_no_plastic_strain(error)
      !! Small volumetric increment from a hydrostatic state well inside the yield
      !! surface. The integrator should take a purely elastic step with no return
      !! mapping, leaving plastic strain zero and F < 0.
      type(error_type), allocatable, intent(out) :: error

      type(mcss_model_t) :: model
      real(wp) :: sig(6), deps(6)

      call make_model(model, FACTOR_STD)

      sig  = [-100.0_wp, -100.0_wp, -100.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
      deps = 1.0e-5_wp * [1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      call euler_substep(model, sig, deps, &
                         integrator_params_t(ftol=YIELD_TOL, stol=1.0e-4_wp, dt_min=1.0e-9_wp))

      call check(error, all(model%state%eps_p == 0.0_wp), .true., &
                 more="elastic step: plastic strain should remain zero")
      if (allocated(error)) return

      call check(error, model%yield_fn(sig) < 0.0_wp, .true., &
                 more="elastic step: stress should remain inside yield surface")
   end subroutine test_elastic_no_plastic_strain

   ! ---------------------------------------------------------------------------
   ! Test 2: plastic step — stress returned to yield surface
   ! ---------------------------------------------------------------------------

   subroutine test_plastic_returns_to_surface(error)
      !! Load from just inside the yield surface with a deviatoric increment whose
      !! elastic predictor overshoots. After integration: stress is on/inside the
      !! yield surface, plastic strain has accumulated, and cohesion has softened.
      type(error_type), allocatable, intent(out) :: error

      type(mcss_model_t) :: model
      real(wp) :: sig(6), deps(6)
      real(wp) :: p_val, J_yield

      call make_model(model, FACTOR_STD)

      ! 90% of the yield deviatoric stress at lode ≈ 0, p = -200 kPa
      p_val   = -200.0_wp
      J_yield = C_PEAK * cos(PHI_PEAK) - p_val * sin(PHI_PEAK)
      sig = [p_val + 0.9_wp * J_yield, p_val, p_val - 0.9_wp * J_yield, &
             0.0_wp, 0.0_wp, 0.0_wp]

      ! Pure deviatoric increment — elastic predictor crosses the yield surface
      deps = [1.0e-3_wp, 0.0_wp, -1.0e-3_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      call euler_substep(model, sig, deps, &
                         integrator_params_t(ftol=YIELD_TOL, stol=1.0e-4_wp, dt_min=1.0e-9_wp))

      ! Slightly looser than ftol to allow for drift-correction residual
      call check(error, model%yield_fn(sig) <= 1.0e-6_wp, .true., &
                 more="plastic step: stress should be on/inside yield surface after return mapping")
      if (allocated(error)) return

      call check(error, any(model%state%eps_p /= 0.0_wp), .true., &
                 more="plastic step: plastic strain should be non-zero")
      if (allocated(error)) return

      call check(error, model%state%c < C_PEAK, .true., &
                 more="plastic step: cohesion should have softened from peak")
   end subroutine test_plastic_returns_to_surface

   ! ---------------------------------------------------------------------------
   ! Test 3: large strain — state converges toward residual
   ! ---------------------------------------------------------------------------

   subroutine test_softening_to_residual(error)
      !! Thirty deviatoric increments from just inside the yield surface.
      !! Accumulated plastic strain drives exponential softening:
      !! c → c_res and φ → φ_res within engineering tolerance.
      type(error_type), allocatable, intent(out) :: error

      type(mcss_model_t) :: model
      real(wp) :: sig(6), deps(6), statev(9)
      real(wp) :: p_val, J_yield
      integer  :: i

      call make_model(model, FACTOR_STD)

      p_val   = -200.0_wp
      J_yield = C_PEAK * cos(PHI_PEAK) - p_val * sin(PHI_PEAK)
      sig  = [p_val + 0.9_wp * J_yield, p_val, p_val - 0.9_wp * J_yield, &
              0.0_wp, 0.0_wp, 0.0_wp]
      deps = [5.0e-3_wp, 0.0_wp, -5.0e-3_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      call mcss_save_state(model, statev)

      do i = 1, 30
         call mcss_load_state(model, statev)
         call euler_substep(model, sig, deps, &
                            integrator_params_t(ftol=YIELD_TOL, stol=1.0e-4_wp, dt_min=1.0e-9_wp))
         call mcss_save_state(model, statev)
      end do

      call check(error, abs(model%state%c - C_RES) < 0.5_wp, .true., &
                 more="residual: cohesion should be near c_res after large plastic deformation")
      if (allocated(error)) return

      call check(error, abs(model%state%phi - PHI_RES) < 1.0e-2_wp, .true., &
                 more="residual: friction angle should be near phi_res after large plastic deformation")
   end subroutine test_softening_to_residual

end module mod_test_mcss_integration_suite
