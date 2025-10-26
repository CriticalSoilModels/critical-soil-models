module MOD_NAMC_ESM
   !**********************************************************************
   !
   ! Module: Contains all functions and subroutines required for the Non-Associative Mohr-Coulomb constitutive model
   !
   !   Info: Constitutive model was presented in Constitutive modelling of non-cohesive soils under high-strain rates DOI: 10.1680/jgeot.21.00192
   !         Current implementation was modified to use Ortiz-Simo yield surface correction
   !
   ! Note: Integer type and real type is not specified in this module. This was done because the ESMs are usually compiled into external .dlls
   !        and the type wouldn't be specfied there
   ! TODO: Add integer and real type. Also determine if doubles are really needed for this computation. As they are cast into reals does it matter?
   !
   !     $Revision: ????? $
   !     $Date: 2024-01-19 10:34 +0500 (WaveHello, 28 Dec 2023) $
   !
   !**********************************************************************

   use kind_precision_module, only: Real_Type => dp
   use mod_bool_helper, only: dbltobool, logic2dbl
   use mod_SRMC, only: SRMC_HSR
   use mod_array_helper, only: reorder_real_array

   ! implicit none

   private ! Makes all function private to this module (No other modules can get access)
   public UMAT_NAMC ! Overides private status for specific subroutine
contains


   Subroutine UMAT_NAMC(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, &
      DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED,&
      CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, COORDS, &
      drot, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, &
      LAYER, KSPT, KSTEP, KINC)

      implicit none

      !Defining inputs
      ! Arguments:
      !          I/O  Type
      !  PROPS    I   R()  : List with model parameters
      !  DSTRAN   I   R()  : Strain increment
      !  DTIME    I   R()  : Time increment
      !  DDSDDE   O   R(,) : Material stiffness matrix
      !  STRESS  I/O  R()  : stresses
      !  STATEV  I/O  R()  : state variables
      !

      integer, intent(in) :: NSTATEV, NPROPS, NPT
      integer :: NTENS
      integer, intent(in) :: NOEL
      real (Real_Type), dimension(NTENS), intent(inout) :: STRESS
      real(Real_Type), dimension(NSTATEV), intent(inout) :: STATEV
      real(Real_Type), dimension(NTENS, NTENS), intent(inout) :: DDSDDE
      real(Real_Type), intent(in) :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
      real(Real_Type), intent(in) :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
      real(Real_Type), dimension(NTENS), intent(in) :: DDSDDT
      real(Real_Type), dimension(NTENS), intent(in) :: DRPLDE
      real(Real_Type)                               :: DRPLDT
      real(Real_Type), dimension(NTENS), intent(in) :: STRAN
      real(Real_Type), dimension(NTENS), intent(inout) :: DSTRAN ! TODO: Change this back to intent(in) - 
                                                                 ! Making it intent(inout) to reorder it for incremental driver
      real(Real_Type), dimension(2), intent(in) :: TIME
      real(Real_Type), dimension(1), intent(in) :: PREDEF
      real(Real_Type), dimension(1), intent(in) :: DPRED
      real(Real_Type), dimension(NPROPS), intent(in) :: PROPS
      real(Real_Type), dimension(3), intent(in) :: COORDS
      real(Real_Type), dimension(3,3), intent(in) :: DFGRD0
      real(Real_Type), dimension(3,3), intent(in) :: DFGRD1, drot
      REAL(Real_Type), intent(in) :: PNEWDT,  TEMP, DTEMP, CELENT
      double precision, intent(in) :: DTIME
      character(len = 80), intent(in):: CMNAME
      integer, intent(in) :: NDI, NSHR, LAYER, KSPT, KSTEP, KINC

      ! Local variables:
      !
      !  DE        : Linear Elastic constitutive matrix
      !  dSig	     : Stress increment vector
      !  Sig	     : Stress vector
      !  dEpsE     : Elastic strain increment vector
      !  dEpsP     : Plastic strain increment vector
      !  dEps      : Total strain increment vector
      !  EpsE      : Elastic strain vector
      !  EpsP      : Plastic strain vector
      !  Eps       : Total strain vector
      !  EpsRate	 : Total strain rate tensor
      !
      integer     			:: N_S, N_i
      real(Real_Type), dimension(6,6) :: DE
      real(Real_Type), dimension(6)   :: Sig
      real(Real_Type), dimension(6)   :: dEpsP
      real(Real_Type), dimension(6)   :: EpsP, ERate
      real(Real_Type), dimension(6) :: dEpsE, EpsE, dEps, Eps
      logical :: switch_smooth, switch_original, switch_yield
      real(Real_Type) :: G_0, enu, eM_tc, eN, D_min, eh, alpha_G, alpha_K, alpha_D, D_part, G_s, RefERate !SSMC props local variables, (props)
      real(Real_Type) :: G, bk, eta_y, DP, eI_coeff, Sum_rate  !SSMC state variables (statv)
      double precision :: Error_yield_1, Error_yield_2, Error_Euler_max, Error_Yield_last, Error_Yield_max
      double precision :: F1, F2, bK_0, FTOL
      integer :: i, max_stress_iters
      integer :: switch_plastic_integration
      integer, parameter :: A3D_voigt_order(6) = [1, 2, 3, 4, 6, 5]

      ! Viscoplastic NAMC model with smoothed outer surface in pi plane
      !
      ! Contents of PROPS(10) NAMC with HSR
      !  1 : G_0				shear modulus
      !  2 : enu				Poisson's ratio
      !  3 : eM_tc            Critical stress ratio for triaxial compression
      !  4 : eN               Nova's vol coupling coefficient
      !  5 : D_min			Minimum dilation
      !  6 : eh				hardening parameter
      !  7 : alpha_G			Shear modulus viscosity factor
      !  8 : alpha_K			Bulk modulus viscosity factor
      !  9 : alpha_D			dilation viscosity
      ! 10 : D_part			Particle diameter
      ! 11 : G_s				Specific gravity
      ! 12 : RefERate			Reference strain rate
      ! 13 : Switch_smooth	Boolean switch for activating strain rate smoothing
      ! 14 : N_S				Degree of smoothening
      ! 15 : switch_original	Changes from Wang's to Olzak&Perzyna consistency
      ! 16 : FTOL             Yield surface tolerance
      ! 17 : max_stress_iters Maximum stress integration iterations
      ! 18 : switch_plastic_integration: Switch that controls which integration scheme is used to split the elastic and plastic portions for the strain increment

      G_0         = PROPS(1)         ! shear modulus
      enu         = PROPS(2)         ! Poisson's ratio
      eM_tc       = PROPS(3)         ! Critical stress ratio
      eN          = PROPS(4)         ! Nova's vol coupling coefficient
      D_min       = PROPS(5)         ! Minimum dilation
      eh          = PROPS(6)         ! hardening parameter
      alpha_G     = PROPS(7)         ! Shear modulus viscosity factor
      alpha_K     = PROPS(8)         ! Bulk modulus viscosity factor

      alpha_D     = PROPS(9)         ! dilation angle viscosity
      
      D_part      = PROPS(10)        ! Associated particle diameter, value in mm
      G_s         = PROPS(11)        ! Specific gravity
      RefERate    = PROPS(12)        ! reference strain rate
      
      call dbltobool(PROPS(13), switch_smooth)  ! switch for activating strain rate smoothening
      
      N_S        = int(PROPS(14))    ! Degree of smoothening
      
      call dbltobool(PROPS(15), switch_original)! (1 for Wang, 0 for Olzak&Perzyna)
      
      FTOL = PROPS(16)
      max_stress_iters = int(PROPS(17))
      switch_plastic_integration = int(PROPS(18))

      if (RefERate==0.0d0) then
         RefERate=2.5e-5
      endif
      
      if ((alpha_D==0.0d0).and.(alpha_G>0.0d0)) then ! in case user uses 0 it will be computed from nu
         alpha_D= alpha_G
      endif
      
      if ((alpha_K==0.0d0).and.(alpha_G>0.0d0)) then ! in case user uses 0 it will be computed from nu
         alpha_K= 2.5*alpha_G
      endif
      
      !TODO max this an error
      if (FTOL<= 1e-10) then
         print *, "A FTOL of 1e-10 is too small"
      end if

      if(max_stress_iters ==0) then
         print *, "max stress iters is zero"
      end if

      G            = STATEV(1)               ! shear modulus
      bK           = STATEV(2)               ! bulk modulus
      eta_y        = STATEV(3)               ! friction ratio
      Dp           = STATEV(4)               ! Dilation
      eI_coeff     = STATEV(5)               ! Inertial coefficient
      call dbltobool(STATEV(6),switch_yield) ! Point is yielding
      do i=1,6
         EpsP(i)   = STATEV(6+i)            ! Plastic strain
      end do
      N_i          = STATEV(13)            ! Current number of strain rate sums
      SUM_rate     = STATEV(14)            ! Current sum of strain rates

      !_____Error control state parameters__________________________________________________
      Error_yield_1=0.0d0                                                                   !
      Error_yield_2=0.0d0                                                                   !
      Error_Euler_max=0.0d0															      !
      Error_Yield_last=0.0d0															      !
      Error_Yield_max=0.0d0																  !
      !							                                                          !
      !_____________________________________________________________________________________!
      if (DTIME==0.0d0) then
         ERate= 0.0d0    ! Current strain rate
      else
         ERate= (1/DTIME)*DSTRAN ! Current strain rate
      end if

      bK_0= 2*G_0*(1+ENU)/(3*(1-2*ENU))

      !***********************************************************************************
      !Call the refined modified Euler algorithm

      ! Change the order of the Strain and Stress vectors to follow the conventiion of Anura3D
      ! I want to make sure that there's something fishy going on there
      ! I'm pretty sure the order program was written in was for incremental driver but I want to make sure

      ! Convert to A3D order
      ! Sig    = reorder_real_array(Sig, A3D_voigt_order)
      ! Stress = reorder_real_array(Stress, A3D_voigt_order)
      ! DSTRAN = reorder_real_array(DSTRAN, A3D_voigt_order)
      ! EpsP   = reorder_real_array(EpsP, A3D_voigt_order)
      ! Erate  = reorder_real_array(Erate, A3D_voigt_order)

      call SRMC_HSR(NOEL, G_0, enu, eM_tc, eN, D_min, eh, alpha_G, alpha_K, alpha_D, D_part, G_s,&
         switch_smooth, RefERate, N_S, switch_original,&
         G, bK, eta_y, Dp, EpsP, eI_coeff, switch_yield, N_i, SUM_rate,&
         DSTRAN, STRESS, Sig, Erate, DTIME,&
         Error_yield_1, Error_yield_2, Error_Euler_max, Error_Yield_last, &
         Error_Yield_max, FTOL, max_stress_iters, switch_plastic_integration)
      !************************************************************************************
      !************************************************************************************
      
      ! Convert back to incremental driver order
      ! Sig    = reorder_real_array(Sig, A3D_voigt_order)
      ! Stress = reorder_real_array(Stress, A3D_voigt_order)
      ! DSTRAN = reorder_real_array(DSTRAN, A3D_voigt_order)
      ! EpsP   = reorder_real_array(EpsP, A3D_voigt_order)
      ! Erate  = reorder_real_array(Erate, A3D_voigt_order)

      !Stress and state variables updating
      Do i=1,NTENS
         STRESS(i) = Sig(i)
      End Do
      STATEV(1) = G
      STATEV(2) = bK
      STATEV(3) = eta_y
      STATEV(4) = Dp
      STATEV(5) =	eI_coeff
      STATEV(6) = logic2dbl(switch_yield)
      do i=1,6
         STATEV(6+i) = EpsP(i)
      end do
      STATEV(13)=N_i
      STATEV(14)=SUM_rate

      !_____Error control state parameters__________________________________________________
      !     Comment if not wanted                                                           !
      !_____________________________________________________________________________________!
      ! STATEV(28)=Error_yield_1
      ! STATEV(29)=Error_yield_2
      ! STATEV(30)=Error_Euler_max
      ! STATEV(31)=Error_Yield_last
      ! STATEV(32)=Error_Yield_max

      !************************************************************************************
      !************************************************************************************
      !Tangent stiffness matrix to be returned done by elastic stiffness
      F1  = bK+(4*G/3)
      F2  = bK-(2*G/3)
      DDSDDE = 0.0
      DDSDDE(1:3,1:3) = F2
      DDSDDE(1,1) = F1
      DDSDDE(2,2) = F1
      DDSDDE(3,3) = F1
      DDSDDE(4,4) = G
      DDSDDE(5,5) = G
      DDSDDE(6,6) = G
      !*************************************************************************************
      !End of UMAT
      Return
   end subroutine UMAT_NAMC

end module MOD_NAMC_ESM
