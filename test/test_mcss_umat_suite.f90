!! End-to-end UMAT tests for the Mohr-Coulomb Strain Softening model.
!!
!! These tests exercise the full call stack through the UMAT boundary:
!!
!!   umat_mcss → mcss_from_props → mcss_load_state → integrate_stress
!!             → mcss_save_state → DDSDDE
!!
!! Unlike `test_mcss_integration_suite`, which calls `euler_substep` directly,
!! the tests here call `umat_mcss` with raw PROPS and STATEV arrays and
!! thread STATEV output from one increment as input to the next.

module mod_test_mcss_umat_suite
   use mod_csm_kinds, only: wp
   use mod_umat_mcss, only: umat_mcss
   use testdrive,     only: new_unittest, unittest_type, error_type, check
   implicit none
   private
   public :: collect_mcss_umat_suite

   ! ---------------------------------------------------------------------------
   ! Shared test parameters
   ! ---------------------------------------------------------------------------
   real(wp), parameter :: PI         = 4.0_wp * atan(1.0_wp)
   real(wp), parameter :: DEG_TO_RAD = PI / 180.0_wp

   ! Material parameters
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
   real(wp), parameter :: MAX_ITERS    = 100.0_wp     !! Max substep iters (real for PROPS)
   real(wp), parameter :: DT_MIN       = 1.0e-9_wp    !! Min substep size

   ! STATEV initial values at peak strength
   real(wp), parameter :: PHI_PEAK = PHI_PEAK_DEG * DEG_TO_RAD  !! [rad]
   real(wp), parameter :: PHI_RES  = PHI_RES_DEG  * DEG_TO_RAD  !! [rad]
   real(wp), parameter :: PSI_PEAK = PSI_PEAK_DEG * DEG_TO_RAD  !! [rad]

   ! UMAT dimension arguments (full 3D)
   integer, parameter :: NDI_3D    = 3
   integer, parameter :: NSHR_3D   = 3
   integer, parameter :: NTENS_3D  = 6
   integer, parameter :: NSTATEV   = 9
   integer, parameter :: NPROPS    = 13

contains

   subroutine collect_mcss_umat_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("umat elastic step: STATEV unchanged",          test_umat_elastic_step),       &
         new_unittest("umat plastic step: state evolves",             test_umat_plastic_step),       &
         new_unittest("umat multi-increment: STATEV round-trip",      test_umat_multi_increment)     &
      ]
   end subroutine collect_mcss_umat_suite

   ! ---------------------------------------------------------------------------
   ! Helper: build a standard PROPS array (13 entries)
   ! ---------------------------------------------------------------------------

   subroutine make_props(props)
      !! Populate a PROPS(13) array with the standard test parameters.
      real(wp), intent(out) :: props(NPROPS)

      props(1)  = G_TEST
      props(2)  = NU_TEST
      props(3)  = C_PEAK
      props(4)  = C_RES
      props(5)  = PHI_PEAK_DEG   ! degrees — UMAT converts to radians internally
      props(6)  = PHI_RES_DEG
      props(7)  = PSI_PEAK_DEG
      props(8)  = PSI_RES_DEG
      props(9)  = FACTOR_STD
      props(10) = 0.0_wp         ! integrator switch (ignored)
      props(11) = YIELD_TOL
      props(12) = MAX_ITERS
      props(13) = DT_MIN
   end subroutine make_props

   ! ---------------------------------------------------------------------------
   ! Helper: populate STATEV at peak strength, zero plastic strain
   ! ---------------------------------------------------------------------------

   subroutine make_statev_peak(statev)
      !! Set STATEV(9) to the peak-strength, zero-plastic-strain initial state.
      real(wp), intent(out) :: statev(NSTATEV)

      statev      = 0.0_wp
      statev(1)   = C_PEAK    ! c   [kPa]
      statev(2)   = PHI_PEAK  ! phi [rad]
      statev(3)   = PSI_PEAK  ! psi [rad]
      ! statev(4:9) = 0 (eps_p — already zeroed above)
   end subroutine make_statev_peak

   ! ---------------------------------------------------------------------------
   ! Helper: build dummy UMAT arrays that are not under test
   ! ---------------------------------------------------------------------------

   subroutine make_dummy_umat_args(sse, spd, scd, rpl, ddsddt, drplde, drpldt, &
                                    stran, time, dtime, temp, dtemp,            &
                                    predef, dpred, coords, drot, pnewdt, celent, &
                                    dfgrd0, dfgrd1)
      real(wp), intent(out) :: sse, spd, scd, rpl, drpldt, dtime, temp, dtemp
      real(wp), intent(out) :: pnewdt, celent
      real(wp), intent(out) :: ddsddt(NTENS_3D), drplde(NTENS_3D), stran(NTENS_3D)
      real(wp), intent(out) :: time(2), predef(1), dpred(1), coords(3)
      real(wp), intent(out) :: drot(3,3), dfgrd0(3,3), dfgrd1(3,3)

      sse    = 0.0_wp;  spd  = 0.0_wp;  scd  = 0.0_wp
      rpl    = 0.0_wp;  drpldt = 0.0_wp
      dtime  = 1.0_wp;  temp = 0.0_wp;  dtemp = 0.0_wp
      pnewdt = 1.0_wp;  celent = 1.0_wp

      ddsddt = 0.0_wp;  drplde = 0.0_wp;  stran = 0.0_wp
      time   = 0.0_wp;  predef = 0.0_wp;  dpred = 0.0_wp
      coords = 0.0_wp

      drot   = 0.0_wp
      drot(1,1) = 1.0_wp;  drot(2,2) = 1.0_wp;  drot(3,3) = 1.0_wp

      dfgrd0 = 0.0_wp
      dfgrd0(1,1) = 1.0_wp;  dfgrd0(2,2) = 1.0_wp;  dfgrd0(3,3) = 1.0_wp

      dfgrd1 = dfgrd0
   end subroutine make_dummy_umat_args

   ! ---------------------------------------------------------------------------
   ! Test 1: elastic step — STATEV must remain at peak state
   ! ---------------------------------------------------------------------------

   subroutine test_umat_elastic_step(error)
      !! Start at hydrostatic stress well inside the yield surface.
      !! Apply a small purely volumetric strain increment.
      !!
      !! Expected after the call:
      !!   - STATEV(1) (c) equals C_PEAK  (no softening)
      !!   - STATEV(4:9) (eps_p) all zero (no plastic strain)
      !!   - DDSDDE diagonal entries positive (positive definite check)
      !!   - DDSDDE is symmetric
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: stress(NTENS_3D), statev(NSTATEV), ddsdde(NTENS_3D, NTENS_3D)
      real(wp) :: dstran(NTENS_3D), props(NPROPS)
      real(wp) :: sse, spd, scd, rpl, drpldt, dtime, temp, dtemp, pnewdt, celent
      real(wp) :: ddsddt(NTENS_3D), drplde(NTENS_3D), stran(NTENS_3D)
      real(wp) :: time(2), predef(1), dpred(1), coords(3)
      real(wp) :: drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
      character(80) :: cmname
      integer  :: noel, npt, layer, kspt, kstep, kinc
      integer  :: i, j

      call make_props(props)
      call make_statev_peak(statev)
      call make_dummy_umat_args(sse, spd, scd, rpl, ddsddt, drplde, drpldt, &
                                stran, time, dtime, temp, dtemp,             &
                                predef, dpred, coords, drot, pnewdt, celent, &
                                dfgrd0, dfgrd1)

      ! Hydrostatic stress well inside the yield surface
      stress    = 0.0_wp
      stress(1) = -100.0_wp   ! sig_11
      stress(2) = -100.0_wp   ! sig_22
      stress(3) = -100.0_wp   ! sig_33

      ! Small volumetric increment — stays elastic
      dstran    = 0.0_wp
      dstran(1) = 1.0e-5_wp
      dstran(2) = 1.0e-5_wp
      dstran(3) = 1.0e-5_wp

      ddsdde = 0.0_wp
      cmname = "MCSS"
      noel = 1;  npt = 1;  layer = 1;  kspt = 1;  kstep = 1;  kinc = 1

      call umat_mcss(stress, statev, ddsdde,                             &
                     sse, spd, scd,                                      &
                     rpl, ddsddt, drplde, drpldt,                        &
                     stran, dstran,                                      &
                     time, dtime, temp, dtemp,                           &
                     predef, dpred, cmname,                              &
                     NDI_3D, NSHR_3D, NTENS_3D, NSTATEV,                &
                     props, NPROPS,                                      &
                     coords, drot, pnewdt, celent,                       &
                     dfgrd0, dfgrd1,                                     &
                     noel, npt, layer, kspt, kstep, kinc)

      ! STATEV(1) = c must equal C_PEAK (no softening)
      call check(error, abs(statev(1) - C_PEAK) < 1.0e-10_wp, .true., &
                 more="elastic step: cohesion must remain at peak (no softening)")
      if (allocated(error)) return

      ! STATEV(4:9) = eps_p must all be zero
      call check(error, all(statev(4:9) == 0.0_wp), .true., &
                 more="elastic step: plastic strain components must remain zero")
      if (allocated(error)) return

      ! DDSDDE diagonal entries must be positive
      do i = 1, NTENS_3D
         call check(error, ddsdde(i,i) > 0.0_wp, .true., &
                    more="elastic step: DDSDDE diagonal must be positive")
         if (allocated(error)) return
      end do

      ! DDSDDE must be symmetric: |D(i,j) - D(j,i)| < tolerance
      do i = 1, NTENS_3D
         do j = i+1, NTENS_3D
            call check(error, abs(ddsdde(i,j) - ddsdde(j,i)) < 1.0e-6_wp, .true., &
                       more="elastic step: DDSDDE must be symmetric")
            if (allocated(error)) return
         end do
      end do

   end subroutine test_umat_elastic_step

   ! ---------------------------------------------------------------------------
   ! Test 2: plastic step — state must evolve, DDSDDE must be populated
   ! ---------------------------------------------------------------------------

   subroutine test_umat_plastic_step(error)
      !! Start near the yield surface and apply a deviatoric increment whose
      !! elastic predictor crosses the yield surface.
      !!
      !! Expected after the call:
      !!   - At least one of STATEV(4:9) is non-zero (plastic strain accumulated)
      !!   - STATEV(1) (c) is less than C_PEAK (softening occurred)
      !!   - DDSDDE is symmetric
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: stress(NTENS_3D), statev(NSTATEV), ddsdde(NTENS_3D, NTENS_3D)
      real(wp) :: dstran(NTENS_3D), props(NPROPS)
      real(wp) :: sse, spd, scd, rpl, drpldt, dtime, temp, dtemp, pnewdt, celent
      real(wp) :: ddsddt(NTENS_3D), drplde(NTENS_3D), stran(NTENS_3D)
      real(wp) :: time(2), predef(1), dpred(1), coords(3)
      real(wp) :: drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
      character(80) :: cmname
      integer  :: noel, npt, layer, kspt, kstep, kinc
      integer  :: i, j
      real(wp) :: p_val, j_yield

      call make_props(props)
      call make_statev_peak(statev)
      call make_dummy_umat_args(sse, spd, scd, rpl, ddsddt, drplde, drpldt, &
                                stran, time, dtime, temp, dtemp,             &
                                predef, dpred, coords, drot, pnewdt, celent, &
                                dfgrd0, dfgrd1)

      ! Start at 90% of the yield deviatoric stress at p = -200 kPa, Lode ~ 0
      p_val   = -200.0_wp
      j_yield = C_PEAK * cos(PHI_PEAK) - p_val * sin(PHI_PEAK)

      ! ANURA3D convention in 3D is identical to internal for pure deviatoric stress
      stress    = 0.0_wp
      stress(1) = p_val + 0.9_wp * j_yield   ! sig_11
      stress(2) = p_val                       ! sig_22
      stress(3) = p_val - 0.9_wp * j_yield   ! sig_33

      ! Deviatoric increment — elastic predictor crosses yield surface
      dstran    = 0.0_wp
      dstran(1) =  1.0e-3_wp
      dstran(3) = -1.0e-3_wp

      ddsdde = 0.0_wp
      cmname = "MCSS"
      noel = 1;  npt = 1;  layer = 1;  kspt = 1;  kstep = 1;  kinc = 1

      call umat_mcss(stress, statev, ddsdde,                             &
                     sse, spd, scd,                                      &
                     rpl, ddsddt, drplde, drpldt,                        &
                     stran, dstran,                                      &
                     time, dtime, temp, dtemp,                           &
                     predef, dpred, cmname,                              &
                     NDI_3D, NSHR_3D, NTENS_3D, NSTATEV,                &
                     props, NPROPS,                                      &
                     coords, drot, pnewdt, celent,                       &
                     dfgrd0, dfgrd1,                                     &
                     noel, npt, layer, kspt, kstep, kinc)

      ! Plastic strain must have accumulated
      call check(error, any(statev(4:9) /= 0.0_wp), .true., &
                 more="plastic step: at least one eps_p component must be non-zero")
      if (allocated(error)) return

      ! Cohesion must have softened
      call check(error, statev(1) < C_PEAK, .true., &
                 more="plastic step: cohesion must be less than C_PEAK after yielding")
      if (allocated(error)) return

      ! DDSDDE must be symmetric
      do i = 1, NTENS_3D
         do j = i+1, NTENS_3D
            call check(error, abs(ddsdde(i,j) - ddsdde(j,i)) < 1.0e-4_wp, .true., &
                       more="plastic step: DDSDDE must be symmetric")
            if (allocated(error)) return
         end do
      end do

   end subroutine test_umat_plastic_step

   ! ---------------------------------------------------------------------------
   ! Test 3: multi-increment STATEV round-trip
   ! ---------------------------------------------------------------------------

   subroutine test_umat_multi_increment(error)
      !! Run 5 consecutive UMAT calls, threading STATEV output as the next input.
      !!
      !! Expected:
      !!   - c decreases monotonically toward C_RES
      !!   - phi decreases monotonically toward PHI_RES
      !!   - norm of eps_p increases each step (plastic strain accumulates)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: stress(NTENS_3D), statev(NSTATEV), ddsdde(NTENS_3D, NTENS_3D)
      real(wp) :: dstran(NTENS_3D), props(NPROPS)
      real(wp) :: sse, spd, scd, rpl, drpldt, dtime, temp, dtemp, pnewdt, celent
      real(wp) :: ddsddt(NTENS_3D), drplde(NTENS_3D), stran(NTENS_3D)
      real(wp) :: time(2), predef(1), dpred(1), coords(3)
      real(wp) :: drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
      character(80) :: cmname
      integer  :: noel, npt, layer, kspt, kstep, kinc
      integer  :: step
      real(wp) :: c_prev, phi_prev, epsp_norm_prev, epsp_norm_curr
      real(wp) :: p_val, j_yield

      integer, parameter :: N_STEPS = 5

      call make_props(props)
      call make_statev_peak(statev)
      call make_dummy_umat_args(sse, spd, scd, rpl, ddsddt, drplde, drpldt, &
                                stran, time, dtime, temp, dtemp,             &
                                predef, dpred, coords, drot, pnewdt, celent, &
                                dfgrd0, dfgrd1)

      p_val   = -200.0_wp
      j_yield = C_PEAK * cos(PHI_PEAK) - p_val * sin(PHI_PEAK)

      stress    = 0.0_wp
      stress(1) = p_val + 0.9_wp * j_yield
      stress(2) = p_val
      stress(3) = p_val - 0.9_wp * j_yield

      ! Larger increment to drive more plastic flow per step
      dstran    = 0.0_wp
      dstran(1) =  5.0e-3_wp
      dstran(3) = -5.0e-3_wp

      cmname = "MCSS"
      noel = 1;  npt = 1;  layer = 1;  kspt = 1;  kstep = 1;  kinc = 1

      c_prev          = statev(1)
      phi_prev        = statev(2)
      epsp_norm_prev  = 0.0_wp

      do step = 1, N_STEPS

         ddsdde = 0.0_wp

         call umat_mcss(stress, statev, ddsdde,                             &
                        sse, spd, scd,                                      &
                        rpl, ddsddt, drplde, drpldt,                        &
                        stran, dstran,                                      &
                        time, dtime, temp, dtemp,                           &
                        predef, dpred, cmname,                              &
                        NDI_3D, NSHR_3D, NTENS_3D, NSTATEV,                &
                        props, NPROPS,                                      &
                        coords, drot, pnewdt, celent,                       &
                        dfgrd0, dfgrd1,                                     &
                        noel, npt, layer, kspt, kstep, kinc)

         ! c must decrease monotonically toward C_RES
         call check(error, statev(1) <= c_prev, .true., &
                    more="multi-increment: cohesion must decrease or stay equal each step")
         if (allocated(error)) return

         ! phi must decrease monotonically toward PHI_RES
         call check(error, statev(2) <= phi_prev, .true., &
                    more="multi-increment: friction angle must decrease or stay equal each step")
         if (allocated(error)) return

         ! eps_p norm must grow each step (plastic accumulation)
         epsp_norm_curr = sqrt(sum(statev(4:9)**2))
         call check(error, epsp_norm_curr > epsp_norm_prev, .true., &
                    more="multi-increment: eps_p norm must increase each step")
         if (allocated(error)) return

         c_prev         = statev(1)
         phi_prev       = statev(2)
         epsp_norm_prev = epsp_norm_curr

      end do

      ! After 5 large increments, c must be strictly less than C_PEAK
      call check(error, statev(1) < C_PEAK, .true., &
                 more="multi-increment: cohesion must have softened from peak after 5 steps")
      if (allocated(error)) return

      ! phi must be strictly less than PHI_PEAK
      call check(error, statev(2) < PHI_PEAK, .true., &
                 more="multi-increment: friction angle must have softened from peak after 5 steps")

   end subroutine test_umat_multi_increment

end module mod_test_mcss_umat_suite
