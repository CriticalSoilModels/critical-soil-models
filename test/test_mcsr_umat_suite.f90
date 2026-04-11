!! End-to-end UMAT tests for the MCSR model.
!!
!! These tests exercise the full call stack:
!!   umat_mcsr -> mcsr_from_props -> mcsr_load_state -> mcsr_update_rate_state
!!             -> integrate_stress -> mcsr_save_state -> DDSDDE
!!
!! The UMAT is called with raw PROPS/STATEV arrays, and STATEV output from
!! step N is fed back as input to step N+1, validating physical consistency
!! across multiple increments.

module mod_test_mcsr_umat_suite
   use mod_csm_kinds, only: wp
   use mod_umat_mcsr, only: umat_mcsr
   use testdrive,     only: new_unittest, unittest_type, error_type, check
   implicit none
   private
   public :: collect_mcsr_umat_suite

   ! Model parameters — same as test_mcsr_integration_suite
   real(wp), parameter :: G_0_TEST     = 10000.0_wp   !! [kPa]
   real(wp), parameter :: NU_TEST      = 0.3_wp
   real(wp), parameter :: M_TC_TEST    = 1.2_wp
   real(wp), parameter :: YIELD_TOL    = 1.0e-6_wp
   real(wp), parameter :: K_0_TEST     = &
      2.0_wp * G_0_TEST * (1.0_wp + NU_TEST) / (3.0_wp * (1.0_wp - 2.0_wp * NU_TEST))

   ! Initial stress: isotropic at p = -100 kPa
   real(wp), parameter :: SIG_INIT(6) = &
      [-100.0_wp, -100.0_wp, -100.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

   integer, parameter :: NPROPS_MCSR  = 18
   integer, parameter :: NSTATEV_MCSR = 14
   integer, parameter :: NTENS_MCSR   = 6
   integer, parameter :: NDI_MCSR     = 3
   integer, parameter :: NSHR_MCSR    = 3

contains

   subroutine collect_mcsr_umat_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("umat: elastic step STATEV unchanged",     test_elastic_statev),  &
         new_unittest("umat: plastic step state evolves",        test_plastic_evolves),  &
         new_unittest("umat: multi-increment STATEV round-trip", test_multiincrement)    &
      ]
   end subroutine collect_mcsr_umat_suite

   ! ---------------------------------------------------------------------------
   ! Helper: populate PROPS array
   ! ---------------------------------------------------------------------------
   subroutine make_props(props)
      real(wp), intent(out) :: props(NPROPS_MCSR)

      props(1)  = G_0_TEST    ! G_0      [kPa]
      props(2)  = NU_TEST     ! nu       [-]
      props(3)  = M_TC_TEST   ! M_tc     [-]
      props(4)  = 0.0_wp      ! N        [-]  Nova coupling
      props(5)  = 0.3_wp      ! D_min    [-]
      props(6)  = 2.0_wp      ! h        [-]  hardening parameter
      props(7)  = 0.2_wp      ! alpha_G  [-]
      props(8)  = 0.5_wp      ! alpha_K  [-]
      props(9)  = 0.2_wp      ! alpha_D  [-]
      props(10) = 1.0e-3_wp   ! D_part   [mm]
      props(11) = 2.65_wp     ! G_s      [-]
      props(12) = 1.0e-3_wp   ! ref_e_rate [1/s]
      props(13) = 0.0_wp      ! switch_smooth (0 = no smoothing)
      props(14) = 5.0_wp      ! N_S      [-]
      props(15) = 1.0_wp      ! switch_original (1 = Wang)
      props(16) = YIELD_TOL   ! yield_tol
      props(17) = 100.0_wp    ! max_iters
      props(18) = 0.0_wp      ! integrator switch (ignored)
   end subroutine make_props

   ! ---------------------------------------------------------------------------
   ! Helper: populate initial STATEV array
   ! ---------------------------------------------------------------------------
   subroutine make_statev(statev)
      real(wp), intent(out) :: statev(NSTATEV_MCSR)

      statev       = 0.0_wp
      statev(1)    = G_0_TEST   ! G   [kPa]
      statev(2)    = K_0_TEST   ! K   [kPa]
      statev(3)    = M_TC_TEST  ! eta_y = M_tc at isotropic state
      statev(4)    = 0.0_wp     ! dilation
      statev(5)    = 1.0e-3_wp  ! I_coeff (at reference rate)
      ! statev(6:14) = 0 already set above
   end subroutine make_statev

   ! ---------------------------------------------------------------------------
   ! Helper: build dummy UMAT scalar/array arguments
   ! ---------------------------------------------------------------------------
   subroutine make_dummy_args(SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
                               STRAN, TIME, TEMP, DTEMP, PREDEF, DPRED,    &
                               COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1)
      real(wp), intent(out) :: SSE, SPD, SCD, RPL, DRPLDT, TEMP, DTEMP, PNEWDT, CELENT
      real(wp), intent(out) :: DDSDDT(NTENS_MCSR), DRPLDE(NTENS_MCSR)
      real(wp), intent(out) :: STRAN(NTENS_MCSR), TIME(2)
      real(wp), intent(out) :: PREDEF(1), DPRED(1), COORDS(3)
      real(wp), intent(out) :: DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)

      integer :: i

      SSE    = 0.0_wp;  SPD    = 0.0_wp;  SCD   = 0.0_wp
      RPL    = 0.0_wp;  DRPLDT = 0.0_wp;  TEMP  = 0.0_wp
      DTEMP  = 0.0_wp;  PNEWDT = 1.0_wp;  CELENT = 1.0_wp

      DDSDDT  = 0.0_wp;  DRPLDE = 0.0_wp
      STRAN   = 0.0_wp;  TIME   = 0.0_wp
      PREDEF  = 0.0_wp;  DPRED  = 0.0_wp;  COORDS = 0.0_wp

      DROT   = 0.0_wp
      DFGRD0 = 0.0_wp;  DFGRD1 = 0.0_wp
      do i = 1, 3
         DROT(i,i)   = 1.0_wp
         DFGRD0(i,i) = 1.0_wp
         DFGRD1(i,i) = 1.0_wp
      end do
   end subroutine make_dummy_args

   ! ---------------------------------------------------------------------------
   ! Test 1: elastic step — eps_p must remain zero, DDSDDE diagonal > 0
   ! ---------------------------------------------------------------------------
   subroutine test_elastic_statev(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: STRESS(NTENS_MCSR), STATEV(NSTATEV_MCSR), DDSDDE(NTENS_MCSR, NTENS_MCSR)
      real(wp) :: PROPS(NPROPS_MCSR), DSTRAN(NTENS_MCSR)
      real(wp) :: SSE, SPD, SCD, RPL, DRPLDT, TEMP, DTEMP, PNEWDT, CELENT
      real(wp) :: DDSDDT(NTENS_MCSR), DRPLDE(NTENS_MCSR)
      real(wp) :: STRAN(NTENS_MCSR), TIME(2), PREDEF(1), DPRED(1), COORDS(3)
      real(wp) :: DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)
      character(80) :: CMNAME
      integer :: NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      real(wp), parameter :: DTIME = 0.01_wp
      integer  :: i

      call make_props(PROPS)
      call make_statev(STATEV)
      call make_dummy_args(SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
                           STRAN, TIME, TEMP, DTEMP, PREDEF, DPRED,    &
                           COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1)

      STRESS = SIG_INIT
      CMNAME = '                                                                                '

      ! Small volumetric increment (elastic — well inside yield surface)
      DSTRAN = 1.0e-5_wp * [1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      NOEL = 1;  NPT = 1;  LAYER = 1;  KSPT = 1;  KSTEP = 1;  KINC = 1
      DDSDDE = 0.0_wp

      call umat_mcsr(STRESS, STATEV, DDSDDE,           &
                     SSE, SPD, SCD,                     &
                     RPL, DDSDDT, DRPLDE, DRPLDT,      &
                     STRAN, DSTRAN,                     &
                     TIME, DTIME, TEMP, DTEMP,          &
                     PREDEF, DPRED, CMNAME,             &
                     NDI_MCSR, NSHR_MCSR, NTENS_MCSR,  &
                     NSTATEV_MCSR,                      &
                     PROPS, NPROPS_MCSR,                &
                     COORDS, DROT, PNEWDT, CELENT,      &
                     DFGRD0, DFGRD1,                    &
                     NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

      ! eps_p (STATEV 7:12) must remain zero for a purely elastic step
      call check(error, all(STATEV(7:12) == 0.0_wp), .true., &
                 more="elastic step: plastic strains must remain zero")
      if (allocated(error)) return

      ! DDSDDE diagonal entries must be positive
      do i = 1, NTENS_MCSR
         call check(error, DDSDDE(i,i) > 0.0_wp, .true., &
                    more="elastic step: DDSDDE diagonal must be positive")
         if (allocated(error)) return
      end do

   end subroutine test_elastic_statev

   ! ---------------------------------------------------------------------------
   ! Test 2: plastic step — at least one eps_p component non-zero,
   !         DDSDDE must be symmetric
   ! ---------------------------------------------------------------------------
   subroutine test_plastic_evolves(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: STRESS(NTENS_MCSR), STATEV(NSTATEV_MCSR), DDSDDE(NTENS_MCSR, NTENS_MCSR)
      real(wp) :: PROPS(NPROPS_MCSR), DSTRAN(NTENS_MCSR)
      real(wp) :: SSE, SPD, SCD, RPL, DRPLDT, TEMP, DTEMP, PNEWDT, CELENT
      real(wp) :: DDSDDT(NTENS_MCSR), DRPLDE(NTENS_MCSR)
      real(wp) :: STRAN(NTENS_MCSR), TIME(2), PREDEF(1), DPRED(1), COORDS(3)
      real(wp) :: DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)
      character(80) :: CMNAME
      integer :: NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      real(wp), parameter :: DTIME = 0.01_wp
      integer  :: i, j
      real(wp) :: asym

      call make_props(PROPS)
      call make_statev(STATEV)
      call make_dummy_args(SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
                           STRAN, TIME, TEMP, DTEMP, PREDEF, DPRED,    &
                           COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1)

      STRESS = SIG_INIT
      CMNAME = '                                                                                '

      ! Deviatoric increment large enough to cause yielding (same as integration tests)
      DSTRAN = [6.0e-3_wp, 0.0_wp, -6.0e-3_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      NOEL = 1;  NPT = 1;  LAYER = 1;  KSPT = 1;  KSTEP = 1;  KINC = 1
      DDSDDE = 0.0_wp

      call umat_mcsr(STRESS, STATEV, DDSDDE,           &
                     SSE, SPD, SCD,                     &
                     RPL, DDSDDT, DRPLDE, DRPLDT,      &
                     STRAN, DSTRAN,                     &
                     TIME, DTIME, TEMP, DTEMP,          &
                     PREDEF, DPRED, CMNAME,             &
                     NDI_MCSR, NSHR_MCSR, NTENS_MCSR,  &
                     NSTATEV_MCSR,                      &
                     PROPS, NPROPS_MCSR,                &
                     COORDS, DROT, PNEWDT, CELENT,      &
                     DFGRD0, DFGRD1,                    &
                     NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

      ! At least one plastic strain component must be non-zero
      call check(error, any(STATEV(7:12) /= 0.0_wp), .true., &
                 more="plastic step: at least one eps_p component must be non-zero")
      if (allocated(error)) return

      ! DDSDDE must be symmetric: |D(i,j) - D(j,i)| < small tolerance
      asym = 0.0_wp
      do i = 1, NTENS_MCSR
         do j = i+1, NTENS_MCSR
            asym = max(asym, abs(DDSDDE(i,j) - DDSDDE(j,i)))
         end do
      end do
      call check(error, asym < 1.0e-6_wp * maxval(abs(DDSDDE)) + 1.0e-10_wp, .true., &
                 more="plastic step: DDSDDE must be symmetric")

   end subroutine test_plastic_evolves

   ! ---------------------------------------------------------------------------
   ! Test 3: multi-increment STATEV round-trip
   ! Run 5 consecutive increments, feeding output STATEV back as input.
   ! Verify: eps_p norm grows each step, G and K stay positive.
   ! ---------------------------------------------------------------------------
   subroutine test_multiincrement(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: STRESS(NTENS_MCSR), STATEV(NSTATEV_MCSR), DDSDDE(NTENS_MCSR, NTENS_MCSR)
      real(wp) :: PROPS(NPROPS_MCSR), DSTRAN(NTENS_MCSR)
      real(wp) :: SSE, SPD, SCD, RPL, DRPLDT, TEMP, DTEMP, PNEWDT, CELENT
      real(wp) :: DDSDDT(NTENS_MCSR), DRPLDE(NTENS_MCSR)
      real(wp) :: STRAN(NTENS_MCSR), TIME(2), PREDEF(1), DPRED(1), COORDS(3)
      real(wp) :: DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)
      character(80) :: CMNAME
      integer :: NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      real(wp), parameter :: DTIME = 0.01_wp
      integer  :: step
      real(wp) :: eps_p_norm_prev, eps_p_norm

      call make_props(PROPS)
      call make_statev(STATEV)
      call make_dummy_args(SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
                           STRAN, TIME, TEMP, DTEMP, PREDEF, DPRED,    &
                           COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1)

      STRESS = SIG_INIT
      CMNAME = '                                                                                '

      ! Deviatoric increment causing plastic flow each step
      DSTRAN = [6.0e-3_wp, 0.0_wp, -6.0e-3_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      NOEL = 1;  NPT = 1;  LAYER = 1;  KSPT = 1;  KSTEP = 1

      eps_p_norm_prev = 0.0_wp

      do step = 1, 5
         KINC   = step
         DDSDDE = 0.0_wp

         call umat_mcsr(STRESS, STATEV, DDSDDE,           &
                        SSE, SPD, SCD,                     &
                        RPL, DDSDDT, DRPLDE, DRPLDT,      &
                        STRAN, DSTRAN,                     &
                        TIME, DTIME, TEMP, DTEMP,          &
                        PREDEF, DPRED, CMNAME,             &
                        NDI_MCSR, NSHR_MCSR, NTENS_MCSR,  &
                        NSTATEV_MCSR,                      &
                        PROPS, NPROPS_MCSR,                &
                        COORDS, DROT, PNEWDT, CELENT,      &
                        DFGRD0, DFGRD1,                    &
                        NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

         ! Shear modulus G (STATEV(1)) must stay positive
         call check(error, STATEV(1) > 0.0_wp, .true., &
                    more="multi-increment: G must remain positive")
         if (allocated(error)) return

         ! Bulk modulus K (STATEV(2)) must stay positive
         call check(error, STATEV(2) > 0.0_wp, .true., &
                    more="multi-increment: K must remain positive")
         if (allocated(error)) return

         eps_p_norm = norm2(STATEV(7:12))

         ! eps_p norm must increase (plastic flow accumulating)
         call check(error, eps_p_norm > eps_p_norm_prev, .true., &
                    more="multi-increment: eps_p norm must grow each step")
         if (allocated(error)) return

         eps_p_norm_prev = eps_p_norm
      end do

   end subroutine test_multiincrement

end module mod_test_mcsr_umat_suite
