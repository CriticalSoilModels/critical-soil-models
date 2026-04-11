!! End-to-end UMAT tests for the NorSand model.
!!
!! These tests exercise the full call stack:
!!   umat_norsand -> norsand_from_props -> norsand_load_state -> pre_step
!!                -> integrate_stress -> norsand_save_state -> DDSDDE
!!
!! The UMAT is called with raw PROPS(13)/STATEV(15) arrays, and STATEV output
!! from step N is fed back as input to step N+1, validating physical consistency
!! across multiple increments.

module mod_test_norsand_umat_suite
   use mod_csm_kinds,   only: wp
   use mod_umat_norsand, only: umat_norsand
   use testdrive,        only: new_unittest, unittest_type, error_type, check
   implicit none
   private
   public :: collect_norsand_umat_suite

   ! Model parameters
   real(wp), parameter :: G_0_TEST     = 10000.0_wp   !! Reference shear modulus [kPa]
   real(wp), parameter :: P_REF_TEST   = 100.0_wp     !! Reference mean stress [kPa]
   real(wp), parameter :: NG_TEST      = 0.5_wp       !! Shear modulus exponent [-]
   real(wp), parameter :: NU_TEST      = 0.3_wp       !! Poisson's ratio [-]
   real(wp), parameter :: M_TC_TEST    = 1.2_wp       !! Critical friction ratio [-]

   ! Computed elastic constants at p_ref
   real(wp), parameter :: G_INIT = G_0_TEST   ! G_0 * (100/100)^0.5 = G_0
   real(wp), parameter :: K_INIT = &
      G_INIT * 2.0_wp * (1.0_wp + NU_TEST) / (3.0_wp * (1.0_wp - 2.0_wp * NU_TEST))

   ! Initial stress: isotropic at p = -100 kPa
   real(wp), parameter :: SIG_INIT(6) = &
      [-100.0_wp, -100.0_wp, -100.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

   integer, parameter :: NPROPS_NS  = 13
   integer, parameter :: NSTATEV_NS = 15
   integer, parameter :: NTENS_NS   = 6
   integer, parameter :: NDI_NS     = 3
   integer, parameter :: NSHR_NS    = 3

contains

   subroutine collect_norsand_umat_suite(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
         new_unittest("umat: elastic step STATEV unchanged",     test_elastic_statev),  &
         new_unittest("umat: plastic step state evolves",        test_plastic_evolves),  &
         new_unittest("umat: multi-increment STATEV round-trip", test_multiincrement)    &
      ]
   end subroutine collect_norsand_umat_suite

   ! ---------------------------------------------------------------------------
   ! Helper: populate PROPS array
   ! ---------------------------------------------------------------------------
   subroutine make_props(props)
      real(wp), intent(out) :: props(NPROPS_NS)

      props(1)  = G_0_TEST    ! G_0      [kPa]
      props(2)  = P_REF_TEST  ! p_ref    [kPa]
      props(3)  = NG_TEST     ! nG       [-]
      props(4)  = NU_TEST     ! nu       [-]
      props(5)  = 0.7_wp      ! e_o      [-]
      props(6)  = 0.9_wp      ! Gamma    [-]
      props(7)  = 0.05_wp     ! lambda_c [-]
      props(8)  = 1.0_wp      ! R        [-]
      props(9)  = M_TC_TEST   ! M_tc     [-]
      props(10) = 0.0_wp      ! N        [-]
      props(11) = 0.0_wp      ! chi_tc   [-]
      props(12) = 100.0_wp    ! H_0      [-]
      props(13) = 0.0_wp      ! H_y      [-]
   end subroutine make_props

   ! ---------------------------------------------------------------------------
   ! Helper: populate initial STATEV array
   ! p_i = p/e (on yield surface at q=0), psi=0, e = Gamma - lambda_c*ln(-p_i)
   ! ---------------------------------------------------------------------------
   subroutine make_statev(statev)
      real(wp), intent(out) :: statev(NSTATEV_NS)

      real(wp) :: p_i_init, e_init, M_i_init
      real(wp), parameter :: GAMMA    = 0.9_wp
      real(wp), parameter :: LAMBDA_C = 0.05_wp

      p_i_init = SIG_INIT(1) / exp(1.0_wp)   ! p_i = p/e => on yield surface
      e_init   = GAMMA - LAMBDA_C * log(-p_i_init)
      M_i_init = M_TC_TEST * 3.0_wp / (3.0_wp + M_TC_TEST)

      statev       = 0.0_wp
      statev(1)    = G_INIT    ! G   [kPa]
      statev(2)    = K_INIT    ! K   [kPa]
      statev(3)    = SIG_INIT(1)  ! p = -100 [kPa]
      statev(4)    = e_init    ! e   [-]
      statev(5)    = 0.0_wp    ! psi [-]
      statev(6)    = 0.0_wp    ! chi_tce [-]
      statev(7)    = p_i_init  ! p_i [kPa]
      statev(8)    = M_i_init  ! M_i [-]
      statev(9)    = 0.0_wp    ! switch_yield
      ! statev(10:15) = 0 already set above

   end subroutine make_statev

   ! ---------------------------------------------------------------------------
   ! Helper: build dummy UMAT scalar/array arguments
   ! ---------------------------------------------------------------------------
   subroutine make_dummy_args(SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
                               STRAN, TIME, TEMP, DTEMP, PREDEF, DPRED,    &
                               COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1)
      real(wp), intent(out) :: SSE, SPD, SCD, RPL, DRPLDT, TEMP, DTEMP, PNEWDT, CELENT
      real(wp), intent(out) :: DDSDDT(NTENS_NS), DRPLDE(NTENS_NS)
      real(wp), intent(out) :: STRAN(NTENS_NS), TIME(2)
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
   ! Test 1: elastic step — eps_p (STATEV 10:15) must stay zero,
   !         DDSDDE diagonal entries must be positive
   ! ---------------------------------------------------------------------------
   subroutine test_elastic_statev(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: STRESS(NTENS_NS), STATEV(NSTATEV_NS), DDSDDE(NTENS_NS, NTENS_NS)
      real(wp) :: PROPS(NPROPS_NS), DSTRAN(NTENS_NS)
      real(wp) :: SSE, SPD, SCD, RPL, DRPLDT, TEMP, DTEMP, PNEWDT, CELENT
      real(wp) :: DDSDDT(NTENS_NS), DRPLDE(NTENS_NS)
      real(wp) :: STRAN(NTENS_NS), TIME(2), PREDEF(1), DPRED(1), COORDS(3)
      real(wp) :: DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)
      character(80) :: CMNAME
      integer :: NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      real(wp), parameter :: DTIME = 0.01_wp
      integer :: i

      call make_props(PROPS)
      call make_statev(STATEV)
      call make_dummy_args(SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
                           STRAN, TIME, TEMP, DTEMP, PREDEF, DPRED,    &
                           COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1)

      STRESS = SIG_INIT
      CMNAME = 'NORSAND' // repeat(' ', 73)

      ! Small volumetric increment (elastic — well inside yield surface)
      DSTRAN = 1.0e-5_wp * [1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      NOEL = 1;  NPT = 1;  LAYER = 1;  KSPT = 1;  KSTEP = 1;  KINC = 1
      DDSDDE = 0.0_wp

      call umat_norsand(STRESS, STATEV, DDSDDE,          &
                        SSE, SPD, SCD,                    &
                        RPL, DDSDDT, DRPLDE, DRPLDT,      &
                        STRAN, DSTRAN,                    &
                        TIME, DTIME, TEMP, DTEMP,         &
                        PREDEF, DPRED, CMNAME,            &
                        NDI_NS, NSHR_NS, NTENS_NS,        &
                        NSTATEV_NS,                       &
                        PROPS, NPROPS_NS,                 &
                        COORDS, DROT, PNEWDT, CELENT,     &
                        DFGRD0, DFGRD1,                   &
                        NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

      ! eps_p (STATEV 10:15) must remain zero for a purely elastic step
      call check(error, all(STATEV(10:15) == 0.0_wp), .true., &
                 more="elastic step: plastic strains must remain zero")
      if (allocated(error)) return

      ! DDSDDE diagonal entries must be positive
      do i = 1, NTENS_NS
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

      real(wp) :: STRESS(NTENS_NS), STATEV(NSTATEV_NS), DDSDDE(NTENS_NS, NTENS_NS)
      real(wp) :: PROPS(NPROPS_NS), DSTRAN(NTENS_NS)
      real(wp) :: SSE, SPD, SCD, RPL, DRPLDT, TEMP, DTEMP, PNEWDT, CELENT
      real(wp) :: DDSDDT(NTENS_NS), DRPLDE(NTENS_NS)
      real(wp) :: STRAN(NTENS_NS), TIME(2), PREDEF(1), DPRED(1), COORDS(3)
      real(wp) :: DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)
      character(80) :: CMNAME
      integer :: NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      real(wp), parameter :: DTIME = 0.01_wp
      integer :: i, j
      real(wp) :: asym

      call make_props(PROPS)
      call make_statev(STATEV)
      call make_dummy_args(SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
                           STRAN, TIME, TEMP, DTEMP, PREDEF, DPRED,    &
                           COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1)

      STRESS = SIG_INIT
      CMNAME = 'NORSAND' // repeat(' ', 73)

      ! Deviatoric increment large enough to cause yielding
      ! q_trial ≈ 2*10000*5e-3*sqrt(3) ≈ 173 kPa >> yield radius ≈ 85.7 kPa
      DSTRAN = [5.0e-3_wp, 0.0_wp, -5.0e-3_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      NOEL = 1;  NPT = 1;  LAYER = 1;  KSPT = 1;  KSTEP = 1;  KINC = 1
      DDSDDE = 0.0_wp

      call umat_norsand(STRESS, STATEV, DDSDDE,          &
                        SSE, SPD, SCD,                    &
                        RPL, DDSDDT, DRPLDE, DRPLDT,      &
                        STRAN, DSTRAN,                    &
                        TIME, DTIME, TEMP, DTEMP,         &
                        PREDEF, DPRED, CMNAME,            &
                        NDI_NS, NSHR_NS, NTENS_NS,        &
                        NSTATEV_NS,                       &
                        PROPS, NPROPS_NS,                 &
                        COORDS, DROT, PNEWDT, CELENT,     &
                        DFGRD0, DFGRD1,                   &
                        NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

      ! At least one plastic strain component must be non-zero
      call check(error, any(STATEV(10:15) /= 0.0_wp), .true., &
                 more="plastic step: at least one eps_p component must be non-zero")
      if (allocated(error)) return

      ! DDSDDE must be symmetric: |D(i,j) - D(j,i)| < small tolerance
      asym = 0.0_wp
      do i = 1, NTENS_NS
         do j = i+1, NTENS_NS
            asym = max(asym, abs(DDSDDE(i,j) - DDSDDE(j,i)))
         end do
      end do
      call check(error, asym < 1.0e-6_wp * maxval(abs(DDSDDE)) + 1.0e-10_wp, .true., &
                 more="plastic step: DDSDDE must be symmetric")

   end subroutine test_plastic_evolves

   ! ---------------------------------------------------------------------------
   ! Test 3: multi-increment STATEV round-trip
   ! Run 5 consecutive increments, feeding output STATEV back as input.
   ! Verify: G (STATEV(1)) > 0, K (STATEV(2)) > 0 each step,
   !         eps_p norm grows each step.
   ! ---------------------------------------------------------------------------
   subroutine test_multiincrement(error)
      type(error_type), allocatable, intent(out) :: error

      real(wp) :: STRESS(NTENS_NS), STATEV(NSTATEV_NS), DDSDDE(NTENS_NS, NTENS_NS)
      real(wp) :: PROPS(NPROPS_NS), DSTRAN(NTENS_NS)
      real(wp) :: SSE, SPD, SCD, RPL, DRPLDT, TEMP, DTEMP, PNEWDT, CELENT
      real(wp) :: DDSDDT(NTENS_NS), DRPLDE(NTENS_NS)
      real(wp) :: STRAN(NTENS_NS), TIME(2), PREDEF(1), DPRED(1), COORDS(3)
      real(wp) :: DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)
      character(80) :: CMNAME
      integer :: NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      real(wp), parameter :: DTIME = 0.01_wp
      integer :: step
      real(wp) :: eps_p_norm_prev, eps_p_norm

      call make_props(PROPS)
      call make_statev(STATEV)
      call make_dummy_args(SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
                           STRAN, TIME, TEMP, DTEMP, PREDEF, DPRED,    &
                           COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1)

      STRESS = SIG_INIT
      CMNAME = 'NORSAND' // repeat(' ', 73)

      ! Deviatoric increment causing plastic flow each step
      DSTRAN = [5.0e-3_wp, 0.0_wp, -5.0e-3_wp, 0.0_wp, 0.0_wp, 0.0_wp]

      NOEL = 1;  NPT = 1;  LAYER = 1;  KSPT = 1;  KSTEP = 1

      eps_p_norm_prev = 0.0_wp

      do step = 1, 5
         KINC   = step
         DDSDDE = 0.0_wp

         call umat_norsand(STRESS, STATEV, DDSDDE,          &
                           SSE, SPD, SCD,                    &
                           RPL, DDSDDT, DRPLDE, DRPLDT,      &
                           STRAN, DSTRAN,                    &
                           TIME, DTIME, TEMP, DTEMP,         &
                           PREDEF, DPRED, CMNAME,            &
                           NDI_NS, NSHR_NS, NTENS_NS,        &
                           NSTATEV_NS,                       &
                           PROPS, NPROPS_NS,                 &
                           COORDS, DROT, PNEWDT, CELENT,     &
                           DFGRD0, DFGRD1,                   &
                           NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

         ! Shear modulus G (STATEV(1)) must stay positive
         call check(error, STATEV(1) > 0.0_wp, .true., &
                    more="multi-increment: G must remain positive")
         if (allocated(error)) return

         ! Bulk modulus K (STATEV(2)) must stay positive
         call check(error, STATEV(2) > 0.0_wp, .true., &
                    more="multi-increment: K must remain positive")
         if (allocated(error)) return

         eps_p_norm = norm2(STATEV(10:15))

         ! eps_p norm must increase (plastic flow accumulating)
         call check(error, eps_p_norm > eps_p_norm_prev, .true., &
                    more="multi-increment: eps_p norm must grow each step")
         if (allocated(error)) return

         eps_p_norm_prev = eps_p_norm
      end do

   end subroutine test_multiincrement

end module mod_test_norsand_umat_suite
