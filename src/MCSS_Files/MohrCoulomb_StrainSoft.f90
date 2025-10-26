!*****************************************************************************
!                                       ____  _____
!           /\                         |___ \|  __ \
!          /  \   _ __  _   _ _ __ __ _  __) | |  | |
!         / /\ \ | '_ \| | | | '__/ _` ||__ <| |  | |
!        / ____ \| | | | |_| | | | (_| |___) | |__| |
!       /_/    \_\_| |_|\__,_|_|  \__,_|____/|_____/
!
!
!	Anura3D - Numerical modelling and simulation of large deformations
!   and soil�water�structure interaction using the material point method (MPM)
!
!	Copyright (C) 2023  Members of the Anura3D MPM Research Community
!   (See Contributors file "Contributors.txt")
!
!	This program is free software: you can redistribute it and/or modify
!	it under the terms of the GNU Lesser General Public License as published by
!	the Free Software Foundation, either version 3 of the License, or
!	(at your option) any later version.
!
!	This program is distributed in the hope that it will be useful,
!	but WITHOUT ANY WARRANTY; without even the implied warranty of
!	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!	GNU Lesser General Public License for more details.
!
!	You should have received a copy of the GNU Lesser General Public License
!	along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!*****************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !███╗░░░███╗░█████╗░░██████╗░██████╗
 !████╗░████║██╔══██╗██╔════╝██╔════╝
 !██╔████╔██║██║░░╚═╝╚█████╗░╚█████╗░
 !██║╚██╔╝██║██║░░██╗░╚═══██╗░╚═══██╗
 !██║░╚═╝░██║╚█████╔╝██████╔╝██████╔╝
 !╚═╝░░░░░╚═╝░╚════╝░╚═════╝░╚═════╝░
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


module MOD_MCSS_ESM
   !**********************************************************************
   !
   ! Module: contains all functions and subroutines required for the Mohr-Coloumb Strain Softening constitutive model
   !
   ! Note: Integer type and real type is not specified in this module. This was done because the ESMs are usually compiled into external .dlls
   !        and the type wouldn't be specfied there
   ! TODO: Add integer and real type. Also determine if doubles are really needed for this computation. As they are cast into reals does it matter?
   !
   !     $Revision: ????? $
   !     $Date: 2023-12-28 11:41 +0500 (WaveHello, 28 Dec 2023) $
   !
   !**********************************************************************
   !TODO: Figure out why only this module created a .mod file
   use mod_MCSS_funcs, only: EndOfStepCorrection, DetermineElasticProportionPegasusMethod, &
      DetermineYieldFunctionValue, &
      CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain, &
      CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain , &
      CalculateEpsPEq, CalculateSofteningParameters, DetermineDSigAndDEpsP, &
      MCSS_Ortiz_Simo_Integration
   
   use mod_array_helper, only: reorder_real_array

   implicit none
   private ! Makes all function private to this module (No other modules can get access)
   public ESM_MohrCoulombStrainSoftening, UMAT_MohrCoulombStrainSoftening ! Overides private status for specific subroutine

contains
   SUBROUTINE ESM_MohrCoulombStrainSoftening(NPT,NOEL,IDSET,STRESS,EUNLOADING,PLASTICMULTIPLIER, &
      DSTRAN,NSTATEV,STATEV,NADDVAR,ADDITIONALVAR,CMNAME, &
      NPROPS,PROPS,NUMBEROFPHASES,NTENS)

      CHARACTER*80 CMNAME
      integer :: NPT, NOEL, IDSET, NSTATEV, NADDVAR, NPROPS, NUMBEROFPHASES, NTENS
      integer :: IStep, TimeStep
      double precision :: Eunloading, PLASTICMULTIPLIER
      double precision :: STRESS(NTENS), DSTRAN(NTENS), STATEV(NSTATEV), ADDITIONALVAR(NADDVAR), &
         PROPS(NPROPS)

      !---Local variables required in standard UMAT
      double precision, dimension(:), allocatable :: ddsddt ! only for fully coupled thermal analysis: variation of stress increment due to temperature
      double precision, dimension(:), allocatable :: drplde ! only for fully coupled thermal analysis: variation of volumetric heat generation due to strain increment
      double precision, dimension(:), allocatable :: stran
      double precision, dimension(:), allocatable :: time
      double precision, dimension(:), allocatable :: predef
      double precision, dimension(:), allocatable :: dpred
      double precision, dimension(:), allocatable :: coords
      double precision, dimension(:,:), allocatable :: ddsdde ! Jacobian matrix of the constitutive model (tangent stiffness matrix in case of MC)
      double precision, dimension(:,:), allocatable :: drot
      double precision, dimension(:,:), allocatable :: dfgrd0
      double precision, dimension(:,:), allocatable :: dfgrd1
      double precision :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
      double precision :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
      double precision :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
      double precision :: pnewdt, dtime, temp, dtemp, celent
      double precision :: Value ! auxiliary variable holding any real valued number
      double precision :: Porosity, WaterPressure, WaterPressure0, GasPressure, GasPressure0, DegreeSaturation


      integer :: ndi, nshr, layer, kspt, kstep, kinc, IDTask

      allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
         coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )

      ! Initialization
      Eunloading = 0.0
      PlasticMultiplier = 0.0

      !Rename additional variables
      Porosity = AdditionalVar(1)
      WaterPressure = AdditionalVar(2)
      WaterPressure0 = AdditionalVar(3)
      GasPressure = AdditionalVar(4)
      GasPressure0 = AdditionalVar(5)
      DegreeSaturation = AdditionalVar(6)
      time(1) = AdditionalVar(7)   !TotalRealTime
      time(2) = AdditionalVar(8)   !OverallTotalTime
      dtime = AdditionalVar(9)     !TimeIncrement
      IStep = AdditionalVar(10)
      TimeStep = AdditionalVar(11)   !Note: Very first time and load step: Istep=1 and TimeStep=1

      IDTask = 0

      IF((IStep==1).and.(TimeStep==1)) IDTask = 1

      IF (IDTask == 1) then ! initialisation of state variables
         STATEV(1)=PROPS(3)
         STATEV(2)=PROPS(5)
         STATEV(3)=PROPS(7)
      END IF ! IDTask = 1

      !---Call the UMAT
      call umat_MohrCoulombStrainSoftening(stress, statev, ddsdde, sse, spd, scd, &
         rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
         dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props,&
         nprops, coords, drot, pnewdt, celent, dfgrd0, &
         dfgrd1, noel, npt, layer, kspt, kstep, kinc)

      !---Definition of Eunloading -> required to define the max time step
      Eunloading = max(ddsdde(1,1),ddsdde(2,2),ddsdde(3,3))
      !---Always define this value to run the simulation

      ! PlasticMultiplier can be given as an output because plastic points can be plotted as a result
      return
   end subroutine ESM_MohrCoulombStrainSoftening

   SUBROUTINE UMAT_MohrCoulombStrainSoftening(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, &
      RPL, DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP,&
      DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS,&
      NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT,&
      LAYER, KSPT, KSTEP, KINC)

      ! Extra variables that are being passed in but not used
      CHARACTER*80 CMNAME
      double precision :: SSE, SPD, SCD, rpl, DDSDDT(NTENS), DRPLDE(NTENS), DRPLDT, STRAN(NTENS), &
         TIME(2), DTIME, TEMP, DTEMP, PREDEF(1), DPRED(1), COORDS(3), DROT(3,3), pnewdt,  &
         celent, DFGRD0(3,3), DFGRD1(3,3)
      integer :: ndi, nshr, noel, npt, layer, kspt, kstep, kinc

      ! Variables that are being used
      integer :: ntens, nstatev, nprops
      double precision :: STRESS(NTENS), STATEV(NSTATEV), DDSDDE(NTENS,NTENS), DSTRAN(NTENS), PROPS(NPROPS)

      !---  Local variables --------
      double precision :: Rad ! Store the value of a radian

      !Soil properties
      double precision :: G, ENU, cp, cr, phip, phir, psip, psir, factor, c, phi, psi, &
         F1, F2
      ! Stress variables
      double precision :: YTOL, Euler_DT_min

      ! Define elastic matrix (DE), stress increment, stress, plastic strain incremenent, and total plastic strain
      double precision :: DE(6,6), dSig(6), Sig(6), dEpsP(6), EpsP(6)

      integer :: i, integration_flag, num_integration_iters, ipl, intGlo
      
      ! Define an array to store the conversion from Anura3D to incremental driver voigt notation
      integer, parameter :: inc_driver_voigt_reorder(6) = [1, 2, 3, 4, 6, 5]

      ! Put the unused variables in an if statement so the compiler doesn't show a warning
      if (.False.) then
         print *, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, STRAN, TIME, DTIME,TEMP, DTEMP,     &
            PREDEF, DPRED, CMNAME, NDI, NSHR, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, &
            NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      end if

      ! Incremental driver uses a different voigt vector notation than Anura3D.
      ! Incremental driver uses: [11, 22, 33, 12, 13, 23] (Incremental driver stores engineering shear strains)
      ! Anura3D uses           : [11, 22, 33, 12, 23, 31] (Anura3D stores engineering shear strain [e11, e22, e33, 2e12, 2e23, 2e31])

      ! Change the order of the stress and strain vectors from the incremental driver input to the Anura3D input
      STRESS = reorder_real_array(STRESS, inc_driver_voigt_reorder)
      STRAN = reorder_real_array(STRAN, inc_driver_voigt_reorder)
      DSTRAN = reorder_real_array(DSTRAN, inc_driver_voigt_reorder) 


      ! Mohr-Coulomb Strain Softening model
      !
      ! Contents of PROPS(9) MCSS
      !  1 : G       shear modulus
      !  2 : ENU     Poisson's ratio
      !  3 : cp      peak cohesion
      !  4 : cr      residual cohesion
      !  5 : phip    peak friction angle
      !  6 : phir    residual friction angle
      !  7 : psip    peak dilation angle
      !  8 : psir    residual dilation angle
      !  9 : factor  shape factor
      ! 10 : integration_flag - Controls the stress integration scheme 0 for Euler and 1 for Ortiz-Simo
      Rad  = 45d0 / datan(1d0)
      !*
      !* ... start correction routine
      !*
      G      = PROPS(1)         ! shear modulus
      ENU    = PROPS(2)         ! Poisson's ratio
      cp     = PROPS(3)         ! peak cohesion
      cr     = PROPS(4)         ! residual cohesion
      phip   = PROPS(5)/Rad     ! peak friction angle (rad)
      phir   = PROPS(6)/Rad     ! residual friction angle (rad)

      if (.False.) then
         print *, "Peak Friction angle", phip
         print *, "Residual Friciton angle", phir
      end if

      psip   = PROPS(7)/Rad     ! peak dilation angle (rad)
      psir   = PROPS(8)/Rad     ! residual dilation angle (rad)
      factor = PROPS(9)         ! shape factor
      integration_flag = int(PROPS(10)) ! Stress Integration flag
      YTOL = PROPS(11) ! Yield surface tolerance
      num_integration_iters = int(PROPS(12)) ! Number of stress integration iterations allowed
      Euler_DT_min = PROPS(13) ! Minimum pseudo-time step for Sloan integration

      c    = STATEV(1)          ! cohesion
      phi  = STATEV(2)          ! friction angle
      psi  = STATEV(3)          ! dilatancy angle
      Do i = 1,NTENS
         EpsP(i) = STATEV(3+i)
      end do

      ipl     =   0
      !*
      ! Fill elastic material matrix
      F1  = 2*G*(1-ENU)/(1-2*ENU)
      F2  = 2*G*( ENU )/(1-2*ENU)
      DE  = 0.0
      DE(1:3,1:3) = F2
      DE(1,1) = F1
      DE(2,2) = F1
      DE(3,3) = F1
      DE(4,4) = G
      DE(5,5) = G
      DE(6,6) = G
      
      !*
      ! elastic stress increment
      dSig = matmul(DE, DSTRAN)

      ! ! elastic stress
      Sig = STRESS + dSig

      if (.False.) then
         print *, "Stress after first elastic update:", Sig
      end if
      call MOHRStrainSoftening(IntGlo,F1,F2,G,cp,cr,phip,phir,psip,psir,factor,c,phi,psi,stress,EpsP,&
         DSTRAN,dEpsP,Sig,IPL, integration_flag, YTOL, num_integration_iters, Euler_DT_min)

      !*
      !* ... stress state parameters update
      !*
      Do i=1,NTENS
         STRESS(i) = Sig(i)
      End Do
      
      if (.False.) then
         print *, "printing the stress:", stress
      end if
      STATEV(1) = c
      STATEV(2) = phi
      STATEV(3) = psi
      Do i = 1,NTENS
         STATEV(3+i) = EpsP(i)
      end do

      !*
      !* ... Tangent stiffness matrix to be returned (done by elastic stiffness)
      !*
      G       =   PROPS(1)       ! G
      ENU     =   PROPS(2)       ! nu
      F1  = 2*G*(1-ENU)/(1-2*ENU)
      F2  = 2*G*( ENU )/(1-2*ENU)
      DDSDDE = 0.0
      DDSDDE(1:3,1:3) = F2
      DDSDDE(1,1) = F1
      DDSDDE(2,2) = F1
      DDSDDE(3,3) = F1
      DDSDDE(4,4) = G
      DDSDDE(5,5) = G
      DDSDDE(6,6) = G
      
      ! Change the order back from incremental driver order to Anura3D order
      STRESS = reorder_real_array(STRESS, inc_driver_voigt_reorder)
      STRAN = reorder_real_array(STRAN, inc_driver_voigt_reorder)
      DSTRAN = reorder_real_array(DSTRAN, inc_driver_voigt_reorder) 

   end SUBROUTINE UMAT_MohrCoulombStrainSoftening

   Subroutine MOHRStrainSoftening(IntGlo,D1,D2, GG,cp,cr,phip,phir, psip,psir,factor,c, &
      phi,psi,Sig0,EpsP,DEps,DEpsP,SigC,IPL, integration_flag, &
      YTOL, num_integration_iters, Euler_DT_min)

      !**********************************************************************
      !
      ! Elastoplastic constitutive model with STRAIN SOFTENING, based on the
      ! MOHR-COULOMB criterion (considering modifications of Abbo & Sloan (1995))
      ! Following Ortiz and Simo (1986) to determine stress update
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: i,n,m,it
      integer, parameter :: ONE =1, ZERO = 0
      double precision :: F,F0,F2 !Evaluation of the Yield function
      double precision :: alpha !Elastic Strain proportion
      double precision :: SSTOL !Tolerance Relative Error
      double precision :: SPTOL !Tolerance Softening parameters
      double precision :: Rn !Relative error
      double precision :: T,DT,T1,beta,DTmin !Substepping parameters
      double precision :: c1,phi1,psi1,c2,phi2,psi2
      double precision :: ctol,phitol,psitol !c,phi,psi tolerances
      double precision :: Dcr,Dphir,Dpsir !Diference between current and residial values
      double precision :: moduleEr,moduleSigDSig
      double precision :: EpsPEq,EpsPEq1,EpsPEq2 !Equivalent Plastic Deformation
      double precision :: DEpsPEq !Derivative Strain in function of Equivalent Plastic Deformation
      double precision, dimension(6) :: SigYield, SigYield2
      double precision, dimension(6) :: DSigPP,DSigP1,DSigP2, DSIGE
      double precision, dimension(6) :: DEpsPP,DEpsPP1,DEpsPP2
      double precision, dimension(6) :: DEpsS,DEpsSS
      double precision, dimension(6) :: EpsP1,EpsP2
      double precision, dimension(6) :: DEpsPEqDPS,DEpsPEqDPS1
      double precision, dimension(6) :: sumSg,Er
      double precision, dimension(3) :: DSPDPEq,DSPDPEq1 !Variation of softening parameters (c,phi,psi) in function of plastic strain
      double precision, dimension(6, 6) :: DE

      !In variables
      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
      integer, intent(in) :: integration_flag, num_integration_iters
      double precision, intent(in) :: D1,D2,GG !Elastic Parameters
      double precision, intent(in) :: cp,cr,phip,phir,psip,psir,factor !Softening parameter
      double precision, intent(in) :: YTOL !Tolerance Error on the Yield surface (10-6 to 10-9)
      double precision, intent(in) :: Euler_DT_min ! Minimum psuedo-time step size

      !Inout variables
      double precision, intent(inout), dimension(6) :: DEps !Incremental total strain
      double precision, intent(inout):: c,phi,psi !cohesion,friction angle and dilatancy angle
      double precision, intent(inout), dimension(6) :: EpsP !Accumulated Plastic Strain
      double precision, intent(inout), dimension(6) :: Sig0 !Initial Stress
      double precision, intent(inout), dimension(6) :: SigC !Final Stress
      double precision, intent(inout), dimension(6) :: DEpsP !Incremental plastic strain

      !Out variables
      integer, intent(out) :: IPL

      if (num_integration_iters == 0) then
         print *, "Number of MCSS stress iterations set to zero (Param 12 in .GOM file)"
      end if

      !Initialization
      DEpsPEq = 0.0d0
      EpsPEq = 0.0d0
      SigYield = 0.0d0
      DEpsP = 0.0d0
      F = 0.0d0
      it = 0

      if (c > cp.or.phi > phip.or.psi > psip) then
         c = min(c,cp)
         phi = min(phi,phip)
         psi = min(psi,psip)
      end if

      if (c < cr.or.phi < phir.or.psi < psir) then
         c = max(c,cr)
         phi = max(phi,phir)
         psi = max(psi,psir)
      end if

      !Tolerances
      SSTOL = 1e-5 !Tolerance Relative Error (10-3 to 10-5)
      SPTOL = 0.01d0 !Tolerance Softening Parameters (0.0001d0)
      ctol = abs(cp-cr)*SPTOL
      phitol = abs(phip-phir)*SPTOL
      psitol = abs(psip-psir)*SPTOL

      if (Euler_DT_min <= 1e-10) then
         print *, "Input DTmin, MCSS param 13 is less then 1e-10, that's too small, DTmin set to 1e-8"
         DTmin = 1e-8
      else
         DTmin = Euler_DT_min
      endif

      !Check the yield function value
      call DetermineYieldFunctionValue(IntGlo,SigC,c,phi,F)

      !If F<0 then the behaviour is elastic --> Return
      if (F <= YTOL) then
         IPL = 0
         return
      end if

      !If F>0, the behaviour is elastoplastic --> Continue
      Dcr = abs(c - cr)
      Dphir = abs(phi - phir)
      Dpsir = abs(psi - psir)

      !Check if we are in residual conditions or in softening conditions
      if (Dcr <= ctol.and.Dphir <= phitol.and.Dpsir <= psitol) then
         IPL = 1 !IPL=1 Residual conditions --> no changes of the strength parameters
         c = cr
         phi = phir
         psi = psir
      else
         IPL = 2 !IPL=2 Softening Conditions --> changes of the strength parameters
      end if

      select case(integration_flag)

       case(ZERO)
         !Determine the proportion (alpha) of the stress increment that lies within the yield function.
         !The PEGASUS ALGORITHM SCHEME FOR CONVENTIONAL ELASTOPLASTIC MODELS has been used
         call DetermineYieldFunctionValue(IntGlo,Sig0,c,phi,F0)

         ! Form the stiffness matix
         DE = 0.0
         DE(1:3,1:3) = D2
         DE(1,1) = D1
         DE(2,2) = D1
         DE(3,3) = D1
         DE(4,4) = GG
         DE(5,5) = GG
         DE(6,6) = GG

         ! Calc the stress increment
         DSigE = matmul(DE, DEps)
         
         if (F0 < -YTOL) then !In this Time increment there is part of elastic behavior
            call DetermineElasticProportionPegasusMethod(IntGlo,Sig0,DSigE,DEps,c,phi,YTOL,alpha)
         else
            alpha = 0.0d0 !In this time increment all is plastic behavior
         end if

         !Calculate the direction of the stress path--> missing
         !It is assumed that the direction is always outside the yield surface.

         !Determine the elastic portion of the stress increment
         DSigE = alpha * DSigE !Correct Incremental Elastic Stress

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Determine the plastic portion of the stress increment.
         !The method used is the MODIFIED EULER INTEGRATION SCHEME with error control
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !Initialise parameters
         SigYield = Sig0 + DSigE !Sigma on the Yield surface
         DEpsS = (1.0d0-alpha) * DEps !Incremental Plastic Strain

         T = 0.0d0
         DT = 1.0d0

         !Start the plastification
         Do while (T <= 1.0d0)

            m = 0 !Counter
            Rn = 100

            call CalculateEpsPEq(EpsP,EpsPEq) !Determine Equivalent plastic Strain (EpsPEq)

            Do while (Rn > SSTOL.and.m < num_integration_iters)
               !1)Calculation of the portion of the plastic strain increment (DEpsPP)
               DEpsSS = DT * DEpsS !Portion of the plastic strain increment

               !Calculate a first estimate of the associated stress
               !hardening/softening parameter changes
               call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,phip,phir,psip,psir,&
                  EpsPEq,DSPDPEq)
               call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP,EpsPEq,DEpsPEqDPS)
               call DetermineDSigAndDEpsP(IntGlo,D1,D2,GG,c,phi,psi,SigYield,DEpsPEqDPS,DSPDPEq,DEpsSS,DSigP1,DEpsPP1)
               EpsP1 = EpsP + DEpsPP1

               call CalculateEpsPEq(EpsP1,EpsPEq1) !Determine Equivalent plastic Strain (EpsPEq)

               !if (IPL == 1) then !Residual conditions --> no changes of the strength parameters
               !    c1 = c
               !    phi1 = phi
               !    psi1 = psi
               !else !IPL=2 Softening Conditions --> changes of the strength parameters
               call CalculateSofteningParameters(EpsPEq1,factor,cp,cr,phip,phir,psip,psir,c1,phi1,psi1)
               !end if

               !2)Calculate a second estimate of the associated stress
               !hardening/softening parameter changes
               SigYield2 = SigYield + DSigP1

               call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,phip,phir,psip,psir,&
                  EpsPEq1,DSPDPEq1)
               call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP1,EpsPEq1,DEpsPEqDPS1)
               call DetermineDSigAndDEpsP(IntGlo,D1,D2,GG,c1,phi1,psi1,SigYield2,DEpsPEqDPS1,DSPDPEq1,DEpsSS,DSigP2,DEpsPP2)
               EpsP2 = EpsP + DEpsPP2

               call CalculateEpsPEq(EpsP2,EpsPEq2) !Determine Equivalent plastic Strain (EpsPEq)

               !if (IPL == 1) then !Residual conditions --> no changes of the strength parameters
               !    c2 = c
               !    phi2 = phi
               !    psi2 = psi
               !else  !IPL=2 Softening Conditions --> changes of the strength parameters
               call CalculateSofteningParameters(EpsPEq2,factor,cp,cr,phip,phir,psip,psir,c2,phi2,psi2)
               !end if

               !3)Obtain a more accurate modified Euler estimate of the changes in stress,
               !plastic strain and hardening/softening parameters
               DSigPP = 0.5d0 * (DSigP1 + DSigP2)

               !Calculation of the relative error
               Er = 0.5d0 * (DSigP1 - DSigP2)
               moduleEr = sqrt(Er(1)*Er(1)+Er(2)*Er(2)+ Er(3)*Er(3)+ Er(4)*Er(4)+Er(5)*Er(5)+Er(6)*Er(6))

               sumSg = SigYield + DSigPP
               moduleSigDSig = sqrt(sumSg(1)*sumSg(1) + sumSg(2)*sumSg(2) + sumSg(3)*sumSg(3)+ &
                  sumSg(4)*sumSg(4) + sumSg(5)*sumSg(5) + sumSg(6)*sumSg(6))

               !Check the relative error (Rn) of the new stresses, with the defined tolerance (SSTOL)
               Rn = (moduleEr/moduleSigDSig)

               ! check whether decreasing of DT is possible, if not exit loop
               if (DT == DTmin) then
                  exit
               end if

               !4)If Rn>SSTOL, the loop is not finished and the substep is recalculated smaller
               if (Rn > SSTOL) then
                  beta = max (0.9d0*(sqrt(SSTOL/Rn)), 0.1d0)
                  beta = min (beta, 1.1d0)
                  DT = max (DT*beta, DTmin)
                  m = m + 1 ! Update counter
               end if

            end do

            !Update the accumulated stresses, plastic strain and softening parameters
            SigYield = SigYield + DSigPP
            DEpsPP = 0.5d0 * (DEpsPP1 + DEpsPP2)
            DEpsP = DEpsP + DEpsPP
            EpsP = EpsP + DEpsPP

            call CalculateEpsPEq(EpsP,EpsPEq) !Determine Equivalent plastic Strain (EpsPEq)

            call CalculateSofteningParameters(EpsPEq,factor,cp,cr,phip,phir,psip,psir,c,phi,psi)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!! END OF STEP CORRECTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Check if we are on/in the yield surface, otherwise we are still outside (F>0)
            !and a correction is needed.
            call DetermineYieldFunctionValue(IntGlo,SigYield,c,phi,F)
            n=0 !Counter
            do while (abs(F) >= YTOL .and. n < 10000) !The correction is needed
               n = n + 1
               call CalculateEpsPEq(EpsP,EpsPEq)             !Determine Equivalent plastic Strain (EpsPEq)
               call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,phip,phir,psip,psir,&
                  EpsPEq,DSPDPEq)
               call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP,EpsPEq,DEpsPEqDPS)
               call EndOfStepCorrection(IntGlo,D1,D2,GG,IPL,F,SigYield,DSPDPEq,DEpsPEqDPS,EpsP,c,phi,psi)
            end do
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !The substep is updated
            T1 = T + DT

            !If T1>1 the calculation is finished
            If (T1 >= 1d0) then
               SigC = SigYield   !Determine Final stresses
               return
            end if

            !If T1<1, calculation of the next substep DT
            beta = min (0.9d0*(sqrt(SSTOL/Rn)), 1.1d0)
            if (m > 1) then ! the previous step failed
               beta = min (beta, 1.0d0)
               DT = beta * DT
               it = it+1
            else
               DT = beta * DT
               it = 0
            end if
            DT = max (DT, DTmin)
            DT = min (DT, 1.0d0-T1)
            T = T1


         end do  !If T=1 the loop is finishedv   !Determine the proportion (alpha) of the stress increment that lies within the yield function.

       case(ONE)

         call MCSS_Ortiz_Simo_Integration(GG, D1, D2, IntGlo, Sig0, c, phi, psi, factor, DEps, EpsP, dEpsP, &
            cr, phir, psir, cp, phip, psip, ctol, phitol, psitol, YTOL      , &
            num_integration_iters)

         ! State parameters {phi, psi, c} updated inside ortiz-simo
         ! EpsP updated inside of the integration
         ! dEpsP updated inside of ortiz-Simo

         ! Update the Stress
         SigC = Sig0

         ! Increment of elastic stress not updated
      end select
   end subroutine MOHRStrainSoftening


end module MOD_MCSS_ESM
