module mod_MCSS_funcs
! Functions that are used for MCSS. Moved them to a new module to make testing easier
   implicit none

contains

   Subroutine EndOfStepCorrection(IntGlo,D1,D2,GG,IPL,F,Sig,DSPDPEq,DEpsPEqDPS,EpsP,c,phi,psi)
      !**********************************************************************
      !
      ! Final correction of the yield surface drift (END OF STEP CORRECTION).
      ! The stresses, the plastic strain and the strength parameters are corrected.
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: i
      double precision :: p,J,Lode,S3TA !Invariants
      double precision :: Lambda,param,c2,phi2,psi2,F2
      double precision :: Denom,A,Ai
      double precision, dimension(2) :: DFDSP
      double precision, dimension(6) :: DPPDSig,DFDSig,Sig2,DEpsP,EpsP2
      double precision, dimension(6) :: Denom1
      double precision, dimension(3) :: Dh
      !In Variables
      integer, intent(in) :: IntGlo,IPL !Global ID of Gauss point or particle
      double precision, intent(in):: D1,D2,GG
      double precision, intent(in), dimension(3) :: DSPDPEq !Derivatives respect Equivalent Plastic Strain
      double precision, intent(in), dimension(6) :: DEpsPEqDPS !Derivatives respect Equivalent Plastic Strain
      !InOut Variables
      double precision, intent(inout):: c,phi,psi
      double precision, intent(inout), dimension(6) :: Sig
      double precision, intent(inout), dimension(6) :: EpsP
      double precision, intent(inout):: F

      call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
      call CalculateDerivativesYieldFunctAndPlasticPotential(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,DPPDSig)
      call CalculateDerivativesYieldFunctSofteningParameters(p,J,Lode,S3TA,c,phi,DFDSP)

      !Parameter A (hardening/softening parameter)
      A = 0.0d0
      Ai = (DFDSP(1)*DSPDPEq(1) + DFDSP(2)*DSPDPEq(2))
      do i=1,6
         A = A + Ai * DEpsPEqDPS(i) * DPPDSig(i)
      end do

      Denom1(1) = DPPDSig(1)*D1 + DPPDSig(2)*D2 + DPPDSig(3)*D2
      Denom1(2) = DPPDSig(1)*D2 + DPPDSig(2)*D1 + DPPDSig(3)*D2
      Denom1(3) = DPPDSig(1)*D2 + DPPDSig(2)*D2 + DPPDSig(3)*D1
      Denom1(4) = DPPDSig(4)*GG
      Denom1(5) = DPPDSig(5)*GG
      Denom1(6) = DPPDSig(6)*GG

      Denom = Denom1(1)*DFDSig(1) + Denom1(2)*DFDSig(2) + &
         Denom1(3)*DFDSig(3) + Denom1(4)*DFDSig(4) + &
         Denom1(5)*DFDSig(5) + Denom1(6)*DFDSig(6) - A

      Lambda = F/Denom !factor correction

      Sig2 = Sig - Lambda * Denom1 ! Sig2 = Sig + fact * Denom1 Stress corrected
      DEpsP = Lambda * DPPDSig
      EpsP2 = EpsP + DEpsP

      if (IPL == 1)then
         Dh = 0.0d0
      else
         param = DEpsPEqDPS(1) * DEpsP(1) + DEpsPEqDPS(2) * DEpsP(2) + DEpsPEqDPS(3) * DEpsP(3) + &
            DEpsPEqDPS(4) * DEpsP(4) + DEpsPEqDPS(5) * DEpsP(5) + DEpsPEqDPS(6) * DEpsP(6)
         Dh(1) = min (DSPDPEq(1)*param, 0.0d0)
         Dh(2) = min (DSPDPEq(2)*param, 0.0d0)
         Dh(3) = min (DSPDPEq(3)*param, 0.0d0)
      end if

      c2 = c + Dh(1)
      phi2 = phi + Dh(2)
      psi2 = psi + Dh(3)

      call DetermineYieldFunctionValue(IntGlo,Sig2,c2,phi2,F2)

      if ((abs(F2) > abs(F)).or.(Denom == 0.0d0)) then !NormalCorrectionScheme
         Denom = 0.0d0
         Denom = DFDSig(1)*DFDSig(1) + DFDSig(2)*DFDSig(2) + &
            DFDSig(3)*DFDSig(3) + DFDSig(4)*DFDSig(4) + &
            DFDSig(5)*DFDSig(5) + DFDSig(6)*DFDSig(6)

         Lambda = F/Denom
         Sig = Sig - Lambda * DFDSig
         DEpsP = Lambda * DPPDSig
         EpsP = EpsP + DEpsP
         call DetermineYieldFunctionValue(IntGlo,Sig,c,phi,F)
      else
         Sig = Sig2
         EpsP = EpsP2
         c = c2
         phi = phi2
         psi = psi2
         F = F2
      end if

   end subroutine EndOfStepCorrection

   Subroutine DetermineElasticProportionPegasusMethod(IntGlo,Sig,DSig,DEps,c,phi,YTOL,alpha)
      !**********************************************************************
      !
      ! The PEGASUS METHOD method is used
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: Its, max_iterations
      double precision :: alpha0,alpha1,F0,F1,F
      double precision, dimension(6) :: Sig0,Sig1,SigNew
      !In variables
      double precision, intent(in), dimension(6) :: Sig, DSig
      double precision, intent(in), dimension(6) :: DEps
      double precision, intent(in) :: c,phi,YTOL
      integer, intent(in) :: IntGlo       !Global ID of Gauss point or particle
      !Out variables
      double precision, intent(out) :: alpha

      alpha0 = 0.0d0
      alpha1 = 1.0d0

      Sig0 = Sig + alpha0*DSig ! = Sig0
      Sig1 = Sig + alpha1*DSig ! = SigE

      call DetermineYieldFunctionValue(IntGlo,Sig0,c,phi,F0)
      call DetermineYieldFunctionValue(IntGlo,Sig1,c,phi,F1)

      F=YTOL+1000
      Its = 0 !Counter
      max_iterations = 1000 ! Maximum newton iterations

      do while (abs(F) > YTOL.and.Its < max_iterations)
         alpha = alpha1 - F1*(alpha1-alpha0)/(F1-F0)
         SigNew = Sig + alpha*DSig

         call DetermineYieldFunctionValue(IntGlo,SigNew,c,phi,F)

         if ((F*F1) < 0.0d0) then
            alpha0 = alpha1
            F0 = F1
         else
            F0 = F1*F0/(F1+F)
         end if

         alpha1 = alpha
         F1 = F
         Its = Its + 1

      end do
      if (Its >= max_iterations) then
         alpha = 0.0d0
      end if
   end subroutine DetermineElasticProportionPegasusMethod

   Subroutine MCSS_Ortiz_Simo_Integration(G, D1, D2, IntGlo, Sig, c, phi, psi, factor, dEps, EpsP, dEpsP,&
      cr, phir, psir, cp, phip, &
      psip, ctol, phitol, psitol, FTOL, max_iterations)
      !**********************************************************************
      ! Function:  To update the stress (Sig) and plastic strain (EpsP)
      ! Follows Ortiz and Simo (1986) https://doi.org/10.1002/nme.1620230303
      !
      ! Last Modified: 11/10/2023
      ! Author: Jonathan J. Moore
      !**********************************************************************
      ! Subroutine: Determines the change in stress (Sig), increment of plastic strain dEpsP

      ! List the input variables
      ! G: Shear modulus
      ! D1, D2: Components of the elastic stiffness  matrix
      ! IntGlo: Global iD of the Gauss point or particle
      ! Sig: Current stress state
      ! c: Cohesion
      ! phi: Friction angle
      ! psi: Dilatancy angle
      ! factor: Softening parameter
      ! dEps: Total strain increment
      ! EpsP:  Accumulated plastic strain
      ! cr: residual cohesion value
      ! phir: residual friction angle
      ! psir: residual dilatancy angle
      ! cp: peak cohesion value
      ! phip: peak friction angle
      ! psip: peak dilatancy angle
      ! ctol: cohesion tolerance around the residual value (cr)
      ! phitol: friction angle tolerance around the residual value (phir)
      ! psitol: dilatnacy angle tolerance around the residual value (psir)
      ! max_iterations: The maximum number of time the gradient descent  method should be used

      ! -------------------------- Variable Definitions --------------------------
      ! ------------- Scalar Values -------------
      ! In
      double precision, intent(in) :: G, D1, D2, factor, cr, phir, psir, cp, phip, psip, ctol, phitol, psitol, FTOL
      integer, intent(in) :: IntGlo, max_iterations

      ! In/Out
      double precision, intent(inout) :: c, phi, psi
      ! Out
      ! double precision, intent(out) ::

      ! ------------- Vector Values -------------
      ! In
      !double precision, intent(inout), dimension(6) ::
      ! In/Out
      double precision, intent(inout), dimension(6) :: Sig, dEps, EpsP, dEpsP

      ! ------------- Local Variables -------------
      ! Variable definitions
      ! F: Yield surface value
      ! cu:  Updated cohesion value
      ! phiu: Updated friction angle value
      ! psiu: updated dilatancy value
      ! p: Mean stress
      ! J: Deviatoric stress
      ! Lode: Lode angle
      ! S3TA: ??
      ! dummyVal_1, dummyVal_2, dummyVal_3: Free variables
      ! H: Hardening parameter (dF/dLambda)
      ! D1, D2: Dummy values  to store stfiffness matrix components
      ! epsPEq: Equivalent plastic  strain (constant scaled norm of the plastic strain)
      ! dLambda: Increment of the plastic multiplier
      ! counter: Track current number of iterations

      ! dummyVec_6d: length 6 free  vector
      ! dEpsPu: Updated increment of plastic strain
      ! EpsPu: Updated value of the total plastic strain
      ! Sigu: Updated stress value
      ! dSigu: Updated increment of stress
      ! m_vec: Normal to the plastic potential (dP/dSig)
      ! n_vec: Normal to the Yield surface (dF/dSig)
      ! DE_m: Elastic stiffness matrix times the plastic potential normal
      ! DEpsPEqDPS: Derivative of the Equivalent plastic strain  wrt to the plastic strain (dEpsPEq/dEpsP)
      ! DSPDPEq: Derivative of the state  parameters wrt the equivalent strain (dXs/dEpsEq)

      ! dummyVec_3d: length 3 free vector

      ! dFdSP: Derivative of the yield surface wrt  to the  state parameters (dF/dXs)

      ! DE: Stiffness matricx

      ! Local scalar values

      double precision :: F, cu, Phiu, Psiu , J, Lode, S3TA, dummyVal_1, dummyVal_2, dummyVal_3, H, epsPEq, dLambda
      double precision :: p
      integer:: counter

      ! Local vector values
      double precision, dimension(6):: dummyVec_6d, dEpsPu, EpsPu, Sigu, dSigu, &
         m_vec, n_vec, DE_m, DEpsPEqDPS

      double precision, dimension(3):: dummyVec_3d, DSPDPEq

      double precision, dimension(2):: dFdSP ! Derivative of the yield function with respect to the softening parameters (phi, c)

      ! Local matrix values
      double precision, dimension(6,6):: DE

      ! -------------------------- Begin Calculations --------------------------


      ! Store variables for updating
      Sigu = Sig
      EpsPu = EpsP

      cu = c ! Updated cohesion
      Phiu = Phi! Updated friction angle
      Psiu = Psi ! Updated dilatancy

      ! Form the stiffness matix
      DE = 0.0
      DE(1:3,1:3) = D2
      DE(1,1) = D1
      DE(2,2) = D1
      DE(3,3) = D1
      DE(4,4) = G
      DE(5,5) = G
      DE(6,6) = G

      ! Keep the State varaibles constant

      ! Calc the elastic predictor for the stresses
      ! (Assumes that all of strain increment is elastic therfore there is no change in the equivalent plastic strain)
      dSigu = matmul(DE, dEps)

      ! Update the stresses
      Sigu = Sigu + dSigu

      ! Evalue the yield surface
      call DetermineYieldFunctionValue(IntGlo, Sigu, cu, phiu, F)

      if (F <= FTOL) then
         ! Prediction of the stress and strain values are correct and the values can be updated and returned

         ! Update Sig, EpsP, dEpsP
         Sig = Sigu
         EpsP(:) = 0
         dEpsP(:) = 0

         ! Update yield surface values
         c = cu
         phi = phiu
         psi = psiu

         ! Exit the subroutine
         return
      end if

      ! Max number of plastic descent iterations
      counter = 0

      do while(abs(F) >= FTOL .and. counter <= max_iterations)
         call CalculateInvariants(IntGlo, Sigu, p, J, Lode, S3TA)

         ! Calc the equivalent plastic strain
         call CalculateEpsPEq(EpsPu, epsPEq)

         ! Calc n_vec, m_vec
         call CalculateDerivativesYieldFunctAndPlasticPotential(Sigu, p, J, Lode, S3TA, cu, phiu, psiu, n_vec, m_vec)

         ! dF/dXs
         call CalculateDerivativesYieldFunctSofteningParameters(p, J, Lode, S3TA, cu, phiu, dFdSP)

         ! dXs/dEpsPEq
         call CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,&
            phip,phir,psip,psir,EpsPEq,DSPDPEq)

         ! Calc dEpsPEq/dEpsP
         call CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsPu, EpsPEq, DEpsPEqDPS)

         ! Calc D * m
         DE_m = matmul(DE, m_vec)

         ! Compute n_vec.DE.m_vec
         dummyVal_1 = dot_product(n_vec, DE_m)

         ! Make a 1x3 vector to store dF/dXs
         dummyVec_3d(:) = 0
         dummyVec_3d(1) = dFdSP(1)
         dummyVec_3d(2) = dFdSP(2)

         ! Calc the dot product between dF/dXs.dXs/dEpsPEq
         dummyVal_2 = dot_product(dummyVec_3d, DSPDPEq)

         ! Calc the dot product between dEpsPEq/dEpsP.dP/dSig
         dummyVal_3 = dot_product(DEpsPEqDPS, m_vec)

         ! Need to calc the hardening/softening parameter (H)
         ! H = dF/dXs.dXs/dEpsEq * dEpsEq/dEpsP.dP/dSig
         H = dummyVal_2 * dummyVal_3

         ! calc dLambda (Increment of the plastic multiplier)
         dLambda = F/(dummyVal_1 - H)
         !dLambda = F/(dummyVal_1)

         ! Compute the stress update
         Sigu = Sigu - dLambda * DE_m

         ! Accumulate plastic strain
         EpsPu = EpsPu + dLambda * m_vec

         ! Calc the equivalent plastic strain
         call CalculateEpsPEq(EpsPu, epsPEq)

         ! Update the state parameters (c, phi, psi)
         call CalculateSofteningParameters(epsPEq,factor,cp,cr,phip,phir,psip,psir,cu,phiu,psiu)

         ! Calc the yield function value
         call DetermineYieldFunctionValue(IntGlo,Sigu,cu,phiu,F)

         ! Update the counter
         counter = counter + 1
      end do

      ! Retun the integrated parameters
      Sig = Sigu
      dEpsP = EpsPu-EpsP
      EpsP = EpsPu

      c = cu
      Phi = Phiu
      Psi = Psiu
      print *, "Number of iterations", counter
   end Subroutine MCSS_Ortiz_Simo_Integration

   Subroutine CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
      !**********************************************************************
      !
      ! Calcuation of the invariants (defined as Abbo & Sloan (1995))
      !
      !**********************************************************************

      implicit none

      !Local variables
      double precision :: Sx,Sy,Sz,SqTxy,SqTyz,SqTxz,suma,h1,h2,J2,J3
      double precision, parameter :: C00000 = 0.0D0
      double precision, parameter :: C00001 = 1.0D0
      double precision, parameter :: C00P16 = 0.166666666666666D0
      double precision, parameter :: C00002 = 2.0D0
      double precision, parameter :: C00003 = 3.0D0
      double precision, parameter :: CP3333 = 0.333333333333333D0
      double precision, parameter :: C00IR3 = 0.577350269189626D0 ! 1/sqrt(3)
      double precision, parameter :: TINY = 0.000000000000001D0
      !In variables
      double precision, intent(in), dimension(6) :: Sig
      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
      !Out variables
      double precision, intent(out) :: p,J,Lode,S3TA !Invariants

      p = C00000
      J = C00000
      Lode = C00000

      !Mean stress (p)
      p = CP3333 * (Sig(1) + Sig(2) + Sig(3))

      !Deviatoric stress (J)
      Sx = Sig(1) - p
      Sy = Sig(2) - p
      Sz = Sig(3) - p

      suma = (Sig(1)-Sig(2))**2 + (Sig(1)-Sig(3))**2 + (Sig(2)-Sig(3))**2

      SqTxy =  Sig(4) * Sig(4)
      SqTyz =  Sig(5) * Sig(5)
      SqTxz =  Sig(6) * Sig(6)

      J2 = C00P16 * suma + SqTxy + SqTyz + SqTxz
      J = SQRT(J2)
      
      J3 = Sx*Sy*Sz + C00002 * Sig(4)*Sig(5)*Sig(6) - Sx*SqTyz - Sy*SqTxz - Sz*SqTxy

      !Lode's angle (Lode)
      if (J2 > C00000) then
         
         ! Calculate S3TA (This is the term inside the arcsin() part of the lode angle )
         h1 = -C00003/(C00002*C00IR3)
         h2 = J3/(J**3)
         S3TA = h1*h2
         
         ! Check if S3TA is outside the allowed bounds of asin ([-1, 1])
         if (S3TA < -C00001) then
            S3TA = -C00001
         else if (S3TA > C00001) then
            S3TA = C00001
         end if

         Lode = CP3333*asin(S3TA)

      else  !Special case of zero deviatoric stress
         Lode = C00000
         S3TA = C00000
      end if

   end subroutine CalculateInvariants

   Subroutine DetermineYieldFunctionValue(IntGlo,Sig,c,phi,F)
      !**********************************************************************
      !
      ! In this subroutine the yield function evaluated is a smooth hyperbolic approximation to the
      ! Mohr-Coulomb yield criterion (Abbo and Sloan, 1995).
      !
      ! The edges of the hexagonal pyramid and the tip have been smoothed.
      ! There are two parameters aSmooth (smoothes the tip) and ATTRAN(smoothes the edges)
      ! In this case aSmooth=0.0005*c*cot(phi) and LodeT=25�.
      ! If aSmooth=0 and LodeT=30� the original Mohr-Coulomb is obtained.
      ! TODO: Store rounded params in module array and pass those values to the
      ! associated functions
      ! There's a couple functions that use these values
      !**********************************************************************

      implicit none

      !Local variables
      double precision ::  p,J,Lode,S3TA !Invariants
      double precision ::  COH, SPHI, CPHI, COTPHI, STA, CTA, K, aSmooth, ASPHI2, SGN, A, B
      double precision, parameter :: C00001 = 1.0d0 !Parameters
      double precision, parameter :: C00003 = 3.0d0
      double precision, parameter :: C00P50 = 0.0005d0
      double precision, parameter :: C00000 = 0.0d0
      double precision, parameter :: C00IR3 = 0.577350269189626d0
      double precision, parameter :: C000P1 = 0.00000000001D0
      !Constants for rounded K function (for LodeT=25)
      !double precision, parameter :: A1 = 1.432052062044227d0
      !double precision, parameter :: A2 = 0.406941858374615d0
      !double precision, parameter :: B1 = 0.544290524902313d0
      !double precision, parameter :: B2 = 0.673903324498392d0
      !double precision, parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=29.5)
      double precision, parameter :: A1 = 7.138654723242414d0
      double precision, parameter :: A2 = 6.112267270920612d0
      double precision, parameter :: B1 = 6.270447753139589d0
      double precision, parameter :: B2 = 6.398760841429403d0
      double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=30)
      !double precision, parameter :: A1 = -138300705.446275
      !double precision, parameter :: A2 = -138300706.472675
      !double precision, parameter :: B1 = -138300706.3123
      !double precision, parameter :: B2 = 0.192450089729875
      !double precision, parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians
      !In variables
      double precision, intent(in), dimension(6) :: Sig
      double precision, intent(in) :: c,phi
      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle

      !Out variables
      double precision, intent(out) :: F

      F = C00000

      !Calculation of the invariants (p',J,Lode)
      call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!Evaluation of the yield function with Smoothing!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Material parameters
      COH = c     !Cohesion
      SPHI = sin(phi)
      CPHI = cos(phi)

      COTPHI = CPHI/SPHI

      aSmooth = C00P50*COH*COTPHI !Smoothing parameter
      
      ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
      
      if (abs(phi) == C00000) then
         ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
      end if

      !Calculate K function
      if (abs(Lode) < ATTRAN) then
         STA = sin(Lode)
         CTA = cos(Lode)
         K = CTA - STA*SPHI*C00IR3
      else
         SGN = SIGN(C00001,Lode)
         A = A1 + A2*SGN*SPHI
         B = B1*SGN + B2*SPHI
         K = A - B*S3TA
      end if

      !Calculate value of Hyperbolic Yield function
      F = p*SPHI + sqrt(J*J*K*K+ASPHI2) - COH*CPHI

   end subroutine DetermineYieldFunctionValue

   Subroutine CalculateDerivativesYieldFunctAndPlasticPotential(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,DPPDSig)
      !**********************************************************************
      !
      ! Calculation of the derivatives of the yield function (F) and the plastic potencial punction (P).
      ! Based on Abbo & Sloan (1995)
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: i
      double precision :: COH, SPHI, CPHI, TPHI, COTPHI, STA, CTA, A, B,&
         D, aSmooth, ASPHI2, SGN, T3TA, C3TA, J2, psi2
      double precision ::   K, dKdLode
      double precision :: SPSI, CPSI, TPSI, COTPSI, ASPSI2
      double precision :: i1, i2, Sx, Sy, Sz
      double precision :: DFDp,DFDJ,DFDLode !Derivatives F respect Invariants
      double precision :: DPDp,DPDJ,DPDLode !Derivatives P respect Invariants
      double precision :: C1, C2, C3
      double precision, dimension(6):: DpDSig,DJDSig,DJ3DSig !Derivatives Invariants

      double precision, parameter :: C00001 = 1.0D0 !Parameters
      double precision, parameter :: C000P5 = 0.5D0
      double precision, parameter :: C00P50 = 0.0005D0
      double precision, parameter :: C00000 = 0.0D0
      double precision, parameter :: C00003 = 3.0D0
      double precision, parameter :: C00004 = 4.0D0
      double precision, parameter :: C00002 = 2.0D0
      double precision, parameter :: CP3333 = 0.333333333333333D0
      double precision, parameter :: C00IR3 = 0.577350269189626D0
      double precision, parameter :: C0R3I2 = 0.866025403784439D0
      double precision, parameter :: C000P1 = 0.000000000000001D0
      double precision, parameter :: J0 = 0.001D0
      !Constants for rounded K function (for LodeT=25)
      !double precision, parameter :: A1 = 1.432052062044227d0
      !double precision, parameter :: A2 = 0.406941858374615d0
      !double precision, parameter :: B1 = 0.544290524902313d0
      !double precision, parameter :: B2 = 0.673903324498392d0
      !double precision, parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=29.5)
      double precision, parameter :: A1 = 7.138654723242414d0
      double precision, parameter :: A2 = 6.112267270920612d0
      double precision, parameter :: B1 = 6.270447753139589d0
      double precision, parameter :: B2 = 6.398760841429403d0
      double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=30)
      !double precision, parameter :: A1 = -138300705.446275
      !double precision, parameter :: A2 = -138300706.472675
      !double precision, parameter :: B1 = -138300706.3123
      !double precision, parameter :: B2 = 0.192450089729875
      !double precision, parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians
      !In variables
      double precision, intent(in) ::  c,phi,psi !Soft Parameters
      double precision, intent(in), dimension(6) :: Sig
      !Out variables
      double precision, intent(out), dimension(6) :: DFDSig, DPPDSig !Derivatives respect Sigma
      !Inout variables
      double precision, intent(inout) :: p,J,Lode,S3TA !Invariants

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!! DFDSig = C1*DPDSig + C2*DJDSig + C3*DJ3DSig  !!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Material parameters
      COH = c !Cohesion
      SPHI = sin(phi)
      CPHI = cos(phi)
      COTPHI = CPHI/SPHI
      aSmooth = C00P50*COH*COTPHI !Smoothing parameter
      ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
      if (abs(phi) == C00000) then
         ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
      end if

      if (J == C00000) then
         J2 = C000P1
         J = sqrt(J2)
      else
         J2 = J*J
      end if

      CTA = cos(Lode)
      C3TA = CTA*(C00004*CTA*CTA-C00003)
      T3TA = S3TA/C3TA

      !Calculate K function and its derivative
      if (abs(Lode) < ATTRAN) then
         STA = S3TA/(C00004*CTA*CTA-C00001)
         K = CTA - STA*SPHI*C00IR3
         dKdLode =  - STA - C00IR3*SPHI*CTA
      else
         SGN = SIGN(C00001,Lode) ! It puts the Lode's sign to the number 1
         A = A1 + A2*SGN*SPHI
         B = B1*SGN + B2*SPHI
         K = A - B*S3TA
         dKdLode = - C00003*B*C3TA
      end if

      !Derivative Dp/DSig
      DpDSig(1) = CP3333
      DpDSig(2) = CP3333
      DpDSig(3) = CP3333
      DpDSig(4) = C00000
      DpDSig(5) = C00000
      DpDSig(6) = C00000

      !Derivative DJ/DSig
      i1 = C000P5/J
      if (J < 0.0001) then
         i1 = 0.0d0
      end if
      Sx = Sig(1)-p
      Sy = Sig(2)-p
      Sz = Sig(3)-p

      DJDSig(1) = i1 * Sx
      DJDSig(2) = i1 * Sy
      DJDSig(3) = i1 * Sz
      DJDSig(4) = i1 * C00002 * Sig(4)
      DJDSig(5) = i1 * C00002 * Sig(5)
      DJDSig(6) = i1 * C00002 * Sig(6)

      !Derivative DJ3/DSig
      i2 = CP3333*J*J
      DJ3DSig(1) = (Sy*Sz - Sig(5)*Sig(5) + i2)
      DJ3DSig(2) = (Sx*Sz - Sig(6)*Sig(6) + i2)
      DJ3DSig(3) = (Sx*Sy - Sig(4)*Sig(4) + i2)
      DJ3DSig(4) = C00002*(Sig(5)*Sig(6) - Sz*Sig(4))
      DJ3DSig(5) = C00002*(Sig(6)*Sig(4) - Sx*Sig(5))
      DJ3DSig(6) = C00002*(Sig(4)*Sig(5) - Sy*Sig(6))

      D = J*K/(sqrt(J2*K*K + ASPHI2))

      !C1F
      C1 = SPHI
      !C2F
      C2 = D*K - T3TA*D*dKdLode
      !C3F
      C3 = -C0R3I2*dKdLode*D/(J2*C3TA)

      !DFDSig!
      do i=1,6
         DFDSig(i) = C1*DpDSig(i) + C2*DJDSig(i) + C3*DJ3DSig(i)
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!! DPPDSig = DFDSig (if associated Flow Rule)  !!!!!!!!!!!!!!!!!!!!!!
      !!!!! or
      !!!!! DPPDSig = C1*DPDSig + C2*DJDSig + C3*DJ3DSig  !!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (abs(J) < J0) then
         psi2 = phi - abs(J)*(phi - psi)/J0
      else
         psi2 = psi
      end if

      if (phi == psi2) then !If Associated Flow Rule, then DPPDSig = DFDSig
         DPPDSig = DFDSig

      else !If Non-Associated Flow Rule, then calculate...
         !Material parameters
         SPSI = sin(psi2)
         CPSI = cos(psi2)
         if (SPSI<0.0001) then
            COTPSI=0
         else
            COTPSI = CPSI/SPSI
         end if
         aSmooth = C00P50*COH*COTPSI !Smoothing parameter
         ASPSI2 = aSmooth*aSmooth*SPSI*SPSI
         if (abs(psi2) == C00000)then
            ASPSI2 = C00000
         end if

         !Calculate K function and its derivative
         if (abs(Lode) <= ATTRAN) then
            K = CTA - STA*SPSI*C00IR3
            dKdLode = - STA - C00IR3*SPSI*CTA
         else
            A = A1 + A2*SGN*SPSI
            B = B1*SGN + B2*SPSI
            K = A - B*S3TA
            dKdLode = - C00003*B*C3TA
         end if

         D = J*K/(sqrt(J*J*K*K + ASPSI2))

         !C1F
         C1 = SPSI
         !C2F
         C2 = D*K - T3TA*D*dKdLode
         !C3F
         C3 = -C0R3I2*dKdLode*D/(J2*C3TA)

         !DPPDSig
         do i=1,6
            DPPDSig(i) = C1*DpDSig(i) + C2*DJDSig(i) + C3*DJ3DSig(i)
         end do

      end if

   end subroutine CalculateDerivativesYieldFunctAndPlasticPotential

   Subroutine CalculateDerivativesYieldFunctSofteningParameters(p,J,Lode,S3TA,c,phi,DFDSP)
      !**********************************************************************
      !
      ! Calculation of the derivatives of the yield function (F) with respect the strength parameters
      ! The strength parameters are: cohesion (COH) and friction angle (PHI)
      !
      !**********************************************************************

      implicit none

      !Local variables
      double precision :: COH, SPHI, CPHI, TPHI, COTPHI, STA, CTA, A, B,&
         Denom, Num, aSmooth, ASPHI2, SGN
      double precision :: K, dKdPhi, dadc, dadPhi
      double precision, parameter :: C00001 = 1.0D0 !Parameters
      double precision, parameter :: C00P50 = 0.0005D0
      double precision, parameter :: C00000 = 0.0D0
      double precision, parameter :: C00003 = 3.0D0
      double precision, parameter :: C00002 = 2.0D0
      double precision, parameter :: C00IR3 = 0.577350269189626D0
      double precision, parameter :: C000P1 = 0.00000000001D0
      !Constants for rounded K function (for LodeT=25)
      !double precision, parameter :: A1 = 1.432052062044227d0
      !double precision, parameter :: A2 = 0.406941858374615d0
      !double precision, parameter :: B1 = 0.544290524902313d0
      !double precision, parameter :: B2 = 0.673903324498392d0
      !double precision, parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=29.5)
      double precision, parameter :: A1 = 7.138654723242414d0
      double precision, parameter :: A2 = 6.112267270920612d0
      double precision, parameter :: B1 = 6.270447753139589d0
      double precision, parameter :: B2 = 6.398760841429403d0
      double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=30)
      !double precision, parameter :: A1 = -138300705.446275
      !double precision, parameter :: A2 = -138300706.472675
      !double precision, parameter :: B1 = -138300706.3123
      !double precision, parameter :: B2 = 0.192450089729875
      !double precision, parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians

      !In variables
      double precision, intent(in) :: p,J,Lode,S3TA !Invariants
      double precision, intent(in) :: c,phi !Soft Parameters
      !Out variables
      double precision, intent(out), dimension(2) :: DFDSP !Derivatives respect Soft Parameters


      !Material parameters
      COH = c !Cohesion
      SPHI = sin(phi)
      CPHI = cos(phi)
      COTPHI = CPHI/SPHI

      !Calculate aSmooth and its derivatives
      if (abs(phi) == C00000) then
         COTPHI = C00000
         dadc = C00000
         dadPhi = C00000
      else
         dadc = C00P50*CPHI/SPHI
         dadPhi = - C00P50*COH/(SPHI*SPHI)
      end if
      aSmooth = C00P50*COH*COTPHI !Smoothing parameter
      ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
      if (abs(phi) == C00000) then
         ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
      end if

      !Calculate K function and its derivatives
      if (abs(Lode) <= ATTRAN) then
         STA = sin(Lode)
         CTA = cos(Lode)
         K = CTA - STA*SPHI*C00IR3
         dKdPhi = - C00IR3*CPHI*STA
      else
         SGN = SIGN(C00001,Lode) !It puts the Lode's sign to the number 1
         A = A1 + A2*SGN*SPHI
         B = B1*SGN + B2*SPHI
         K = A - B*S3TA
         dKdPhi = A2*SGN*CPHI - B2*CPHI*S3TA
      end if

      !Operating..
      Denom = (sqrt(J*J*K*K + ASPHI2))
      Num =  J*J*K*dKdPhi + aSmooth*SPHI*SPHI*dadPhi + aSmooth*aSmooth*SPHI*CPHI

      !Derivative DF/Dc
      DFDSP(1) = aSmooth*SPHI*SPHI*dadc/Denom - CPHI

      !Derivative DF/Dphi
      DFDSP(2) = p*CPHI + Num/Denom + COH*SPHI

      if (J <= C00000) then
         DFDSP(1) = - CPHI
         DFDSP(2) = p*CPHI + COH*SPHI
      end if

   end subroutine CalculateDerivativesYieldFunctSofteningParameters

   subroutine CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain(factor,cp,cr,&
      phip,phir,psip,psir,EpsPEq,DSPDPEq)
      !**********************************************************************
      !
      ! Calculation of the derivatives of the strength parameters with respect
      ! the equivalent plastic shear strain
      ! TODO: shorten code length store-> -factor * (exp(-factor*EpsPEq)) as
      !       softening factor
      !**********************************************************************

      implicit none

      ! In Variables
      double precision, intent(in) :: EpsPEq
      double precision, intent(in) :: factor,cp,cr,phip,phir,psip,psir
      ! Out Variables
      double precision, intent(out), dimension(3):: DSPDPEq

      ! Derivative Cohesion respect Equivalent Plastic Strain (Dc/DPEq)
      DSPDPEq(1) = -factor * (cp - cr) * (exp(-factor*EpsPEq))
      ! Derivative Friction angle respect Equivalent Plastic Strain (Dphi/DPEq)
      DSPDPEq(2) = -factor * (phip - phir) * (exp(-factor*EpsPEq))
      ! Derivative Dilatancy angle respect Equivalent Plastic Strain (Dpsi/DPEq)
      DSPDPEq(3) = -factor * (psip - psir) * (exp(-factor*EpsPEq))

   end subroutine CalculateDerivativesStrSoftParamRespectEquivalentPlasticStrain

   Subroutine CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain(EpsP,EpsPEq,DEpsPEqDPS)
      !**********************************************************************
      !
      ! Calculation of the derivatives of the equivalent plastic shear strain
      ! with respect the plastic strain
      !
      !**********************************************************************

      implicit none

      !Local Variables
      double precision :: k1, k2, k3
      double precision :: EpsPM
      double precision, dimension(3) :: EpsDev
      !In Variables
      double precision, intent(in), dimension(6) :: EpsP
      double precision, intent(in) :: EpsPEq
      !Out Variables
      double precision, intent(out), dimension(6):: DEpsPEqDPS

      if (EpsPEq < 0.00000000001d0) then
         k1 = 0.0d0
      else
         k1 = 2.0d0/(3.0d0*EpsPEq)
      end if

      k2 = k1 * 1.0d0/3.0d0
      k3 = k1 * 2.0d0

      EpsPM = k2 * (EpsP(1) + EpsP(2) + EpsP(3))
      EpsDev(1) = EpsP(1)-EpsPM
      EpsDev(2) = EpsP(2)-EpsPM
      EpsDev(3) = EpsP(3)-EpsPM

      DEpsPEqDPS(1) = k2 * ( 2.0d0*EpsDev(1) - EpsDev(2) - EpsDev(3))
      DEpsPEqDPS(2) = k2 * (-EpsDev(1) + 2.0d0*EpsDev(2) - EpsDev(3))
      DEpsPEqDPS(3) = k2 * (-EpsDev(1) - EpsDev(2) + 2.0d0*EpsDev(3))
      DEpsPEqDPS(4) = k3 * EpsP(4)
      DEpsPEqDPS(5) = k3 * EpsP(5)
      DEpsPEqDPS(6) = k3 * EpsP(6)

   end subroutine CalculateDerivativesEquivalentPlasticStrainRespectPlasticStrain

   subroutine CalculateEpsPEq(EpsP,EpsPEq)
      !**********************************************************************
      !
      ! Calculation of the equivalent plastic shear strain
      !
      !**********************************************************************

      implicit none

      !Local variables
      double precision:: EpsPM, C1, C2
      double precision, dimension(3) :: EpsDev
      !In variables
      double precision, intent(in), dimension(6) :: EpsP
      !Out variables
      double precision, intent(out) :: EpsPEq

      !EpsPEq = ((2/3)ep:ep)^(1/2), ep is the deviatoric plastic strain

      EpsPM = (1.0d0/3.0d0) * (EpsP(1) + EpsP(2) + EpsP(3))
      EpsDev(1) = EpsP(1)-EpsPM
      EpsDev(2) = EpsP(2)-EpsPM
      EpsDev(3) = EpsP(3)-EpsPM
      C1 = 2.0d0/3.0d0
      C2 = C1 * 2.0d0

      EpsPEq = sqrt(C1*EpsDev(1)*EpsDev(1) + C1*EpsDev(2)*EpsDev(2) +  C1*EpsDev(3)*EpsDev(3) +&
         C2*EpsP(4)*EpsP(4) + C2*EpsP(5)*EpsP(5) + C2*EpsP(6)*EpsP(6))

   end subroutine CalculateEpsPEq

   Subroutine CalculateSofteningParameters(EpsPEq,factor,cp,cr,phip,phir,psip,psir,c,phi,psi)
      !**********************************************************************
      !
      ! Calculation of strenght parameters (c, phi, psi)
      !
      !**********************************************************************

      implicit none

      !In variables
      double precision, intent(in) :: EpsPEq,factor,cp,cr,phip,phir,psip,psir
      !Out variables
      double precision, intent(out) :: c,phi,psi

      c   = cr + (cp-cr)*exp(-factor*EpsPEq)
      phi = phir + (phip-phir)*exp(-factor*EpsPEq)
      psi = psir + (psip-psir)*exp(-factor*EpsPEq)

   end subroutine CalculateSofteningParameters

   Subroutine DetermineDSigAndDEpsP(IntGlo,D1,D2,GG,c,phi,psi,Sig,DEpsPEqDPS,DSPDPEq,DEps,DSig,DEpsP)
      !**********************************************************************
      !
      ! Calculation of the stress increment and plastic strain increment
      !
      !         dSig = Dep * dEps
      !         dEpsP = Lambda * DPDSig
      ! TODO: Turn the long calculation into a matrix vector product
      !**********************************************************************

      implicit none

      !Local variables
      integer :: i,k
      double precision :: A,Ai,Denom,Fact,LambdaNum,Lambda
      double precision :: p,J,Lode,S3TA !Invariants
      double precision, dimension(6,6) :: Num,Num1,Prod
      double precision, dimension(6) :: Denom1
      double precision, dimension(6) :: DPPDSig !Derivatives Plastic potential respect net stress
      double precision, dimension(6) :: DFDSig !Derivatives Yield function respect net stress
      double precision, dimension(2) :: DFDSP !Derivatives Yield function respect Soft Parameters
      double precision, dimension(6,6) :: Dep !Elastoplastic Constitutive Matrix
      !In Variables
      double precision, intent(in) :: c,phi,psi !Softening parameters
      double precision, intent(in) :: D1,D2,GG !Elastic parameters
      double precision, intent(in), dimension(6):: DEpsPEqDPS
      double precision, intent(in), dimension(6) :: Sig
      double precision, intent(in), dimension(3) :: DSPDPEq !Derivatives respect Equivalent Plastic Strain
      double precision, intent(in), dimension(6) :: DEps
      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
      !Out Variables
      double precision, intent(out), dimension(6) :: DSig
      double precision, intent(out), dimension(6) :: DEpsP

      call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
      call CalculateDerivativesYieldFunctAndPlasticPotential(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,DPPDSig)
      call CalculateDerivativesYieldFunctSofteningParameters(p,J,Lode,S3TA,c,phi,DFDSP)

      !Parameter A (H = -A --> A>0 softening / A<0 hardening)
      A = 0.0d0
      Ai = (DFDSP(1)*DSPDPEq(1) + DFDSP(2)*DSPDPEq(2))
      do i=1,6
         A = A + Ai * DEpsPEqDPS(i) * DPPDSig(i)
      end do

      !Elastoplastic Constitutive Matrix (Dep)
      do i=1,6
         do k=1,6
            Prod(i,k) =  DPPDSig(i) * DFDSig(k)
         end do
      end do

      Num1(1,1) = D1*Prod(1,1) + D2*Prod(2,1) + D2*Prod(3,1)
      Num1(1,2) = D1*Prod(1,2) + D2*Prod(2,2) + D2*Prod(3,2)
      Num1(1,3) = D1*Prod(1,3) + D2*Prod(2,3) + D2*Prod(3,3)
      Num1(1,4) = D1*Prod(1,4) + D2*Prod(2,4) + D2*Prod(3,4)
      Num1(1,5) = D1*Prod(1,5) + D2*Prod(2,5) + D2*Prod(3,5)
      Num1(1,6) = D1*Prod(1,6) + D2*Prod(2,6) + D2*Prod(3,6)

      Num1(2,1) = D2*Prod(1,1) + D1*Prod(2,1) + D2*Prod(3,1)
      Num1(2,2) = D2*Prod(1,2) + D1*Prod(2,2) + D2*Prod(3,2)
      Num1(2,3) = D2*Prod(1,3) + D1*Prod(2,3) + D2*Prod(3,3)
      Num1(2,4) = D2*Prod(1,4) + D1*Prod(2,4) + D2*Prod(3,4)
      Num1(2,5) = D2*Prod(1,5) + D1*Prod(2,5) + D2*Prod(3,5)
      Num1(2,6) = D2*Prod(1,6) + D1*Prod(2,6) + D2*Prod(3,6)

      Num1(3,1) = D2*Prod(1,1) + D2*Prod(2,1) + D1*Prod(3,1)
      Num1(3,2) = D2*Prod(1,2) + D2*Prod(2,2) + D1*Prod(3,2)
      Num1(3,3) = D2*Prod(1,3) + D2*Prod(2,3) + D1*Prod(3,3)
      Num1(3,4) = D2*Prod(1,4) + D2*Prod(2,4) + D1*Prod(3,4)
      Num1(3,5) = D2*Prod(1,5) + D2*Prod(2,5) + D1*Prod(3,5)
      Num1(3,6) = D2*Prod(1,6) + D2*Prod(2,6) + D1*Prod(3,6)

      Num1(4,1) = GG*Prod(4,1)
      Num1(4,2) = GG*Prod(4,2)
      Num1(4,3) = GG*Prod(4,3)
      Num1(4,4) = GG*Prod(4,4)
      Num1(4,5) = GG*Prod(4,5)
      Num1(4,6) = GG*Prod(4,6)

      Num1(5,1) = GG*Prod(5,1)
      Num1(5,2) = GG*Prod(5,2)
      Num1(5,3) = GG*Prod(5,3)
      Num1(5,4) = GG*Prod(5,4)
      Num1(5,5) = GG*Prod(5,5)
      Num1(5,6) = GG*Prod(5,6)

      Num1(6,1) = GG*Prod(6,1)
      Num1(6,2) = GG*Prod(6,2)
      Num1(6,3) = GG*Prod(6,3)
      Num1(6,4) = GG*Prod(6,4)
      Num1(6,5) = GG*Prod(6,5)
      Num1(6,6) = GG*Prod(6,6)



      Num(1,1) = D1*Num1(1,1) + D2*Num1(1,2) + D2*Num1(1,3)
      Num(1,2) = D2*Num1(1,1) + D1*Num1(1,2) + D2*Num1(1,3)
      Num(1,3) = D2*Num1(1,1) + D2*Num1(1,2) + D1*Num1(1,3)
      Num(1,4) = GG*Num1(1,4)
      Num(1,5) = GG*Num1(1,5)
      Num(1,6) = GG*Num1(1,6)

      Num(2,1) = D1*Num1(2,1) + D2*Num1(2,2) + D2*Num1(2,3)
      Num(2,2) = D2*Num1(2,1) + D1*Num1(2,2) + D2*Num1(2,3)
      Num(2,3) = D2*Num1(2,1) + D2*Num1(2,2) + D1*Num1(2,3)
      Num(2,4) = GG*Num1(2,4)
      Num(2,5) = GG*Num1(2,5)
      Num(2,6) = GG*Num1(2,6)

      Num(3,1) = D1*Num1(3,1) + D2*Num1(3,2) + D2*Num1(3,3)
      Num(3,2) = D2*Num1(3,1) + D1*Num1(3,2) + D2*Num1(3,3)
      Num(3,3) = D2*Num1(3,1) + D2*Num1(3,2) + D1*Num1(3,3)
      Num(3,4) = GG*Num1(3,4)
      Num(3,5) = GG*Num1(3,5)
      Num(3,6) = GG*Num1(3,6)

      Num(4,1) = D1*Num1(4,1) + D2*Num1(4,2) + D2*Num1(4,3)
      Num(4,2) = D2*Num1(4,1) + D1*Num1(4,2) + D2*Num1(4,3)
      Num(4,3) = D2*Num1(4,1) + D2*Num1(4,2) + D1*Num1(4,3)
      Num(4,4) = GG*Num1(4,4)
      Num(4,5) = GG*Num1(4,5)
      Num(4,6) = GG*Num1(4,6)

      Num(5,1) = D1*Num1(5,1) + D2*Num1(5,2) + D2*Num1(5,3)
      Num(5,2) = D2*Num1(5,1) + D1*Num1(5,2) + D2*Num1(5,3)
      Num(5,3) = D2*Num1(5,1) + D2*Num1(5,2) + D1*Num1(5,3)
      Num(5,4) = GG*Num1(5,4)
      Num(5,5) = GG*Num1(5,5)
      Num(5,6) = GG*Num1(5,6)

      Num(6,1) = D1*Num1(6,1) + D2*Num1(6,2) + D2*Num1(6,3)
      Num(6,2) = D2*Num1(6,1) + D1*Num1(6,2) + D2*Num1(6,3)
      Num(6,3) = D2*Num1(6,1) + D2*Num1(6,2) + D1*Num1(6,3)
      Num(6,4) = GG*Num1(6,4)
      Num(6,5) = GG*Num1(6,5)
      Num(6,6) = GG*Num1(6,6)



      Denom1(1) = DFDSig(1)*D1 + DFDSig(2)*D2 + DFDSig(3)*D2
      Denom1(2) = DFDSig(1)*D2 + DFDSig(2)*D1 + DFDSig(3)*D2
      Denom1(3) = DFDSig(1)*D2 + DFDSig(2)*D2 + DFDSig(3)*D1
      Denom1(4) = DFDSig(4)*GG
      Denom1(5) = DFDSig(5)*GG
      Denom1(6) = DFDSig(6)*GG

      Denom =   Denom1(1)*DPPDSig(1) + Denom1(2)*DPPDSig(2) + &
         Denom1(3)*DPPDSig(3) + Denom1(4)*DPPDSig(4) + &
         Denom1(5)*DPPDSig(5) + Denom1(6)*DPPDSig(6) - A

      Fact = 1d0/Denom

      !Dep
      Dep(1,1) = D1 - Fact*Num(1,1)
      Dep(1,2) = D2 - Fact*Num(1,2)
      Dep(1,3) = D2 - Fact*Num(1,3)
      Dep(1,4) = -Fact*Num(1,4)
      Dep(1,5) = -Fact*Num(1,5)
      Dep(1,6) = -Fact*Num(1,6)

      Dep(2,1) = D2 - Fact*Num(2,1)
      Dep(2,2) = D1 - Fact*Num(2,2)
      Dep(2,3) = D2 - Fact*Num(2,3)
      Dep(2,4) = -Fact*Num(2,4)
      Dep(2,5) = -Fact*Num(2,5)
      Dep(2,6) = -Fact*Num(2,6)

      Dep(3,1) = D2 - Fact*Num(3,1)
      Dep(3,2) = D2 - Fact*Num(3,2)
      Dep(3,3) = D1 - Fact*Num(3,3)
      Dep(3,4) = -Fact*Num(3,4)
      Dep(3,5) = -Fact*Num(3,5)
      Dep(3,6) = -Fact*Num(3,6)

      Dep(4,1) = -Fact*Num(4,1)
      Dep(4,2) = -Fact*Num(4,2)
      Dep(4,3) = -Fact*Num(4,3)
      Dep(4,4) = GG - Fact*Num(4,4)
      Dep(4,5) = -Fact*Num(4,5)
      Dep(4,6) = -Fact*Num(4,6)

      Dep(5,1) = -Fact*Num(5,1)
      Dep(5,2) = -Fact*Num(5,2)
      Dep(5,3) = -Fact*Num(5,3)
      Dep(5,4) = -Fact*Num(5,4)
      Dep(5,5) = GG - Fact*Num(5,5)
      Dep(5,6) = -Fact*Num(5,6)

      Dep(6,1) = -Fact*Num(6,1)
      Dep(6,2) = -Fact*Num(6,2)
      Dep(6,3) = -Fact*Num(6,3)
      Dep(6,4) = -Fact*Num(6,4)
      Dep(6,5) = -Fact*Num(6,5)
      Dep(6,6) = GG - Fact*Num(6,6)

      !!!!!!!!! Calculate Plastic multipliler(Lambda)!!!!!!!!!!!!!!!!!
      LambdaNum =   Denom1(1)*DEps(1) + Denom1(2)*DEps(2) + &
         Denom1(3)*DEps(3) + Denom1(4)*DEps(4) + &
         Denom1(5)*DEps(5) + Denom1(6)*DEps(6)
      Lambda =  LambdaNum/Denom

      !!!!!!!!! Determine DSig --> (DSig = Dep*dEps) !!!!!!!!!!!
      do i=1,6
         DSig(i) = 0.0d0
         do k=1,6
            DSig(i) =  DSig(i) + Dep(i,k) * DEps(k)
         end do
      end do

      !!!!!!!!! Determine DEpsP --> (DEpsP = Lambda*DPDSig) !!!!!!!!!!!!
      do i=1,6
         DEpsP(i) = Lambda * DPPDSig(i)
      end do

   end subroutine DetermineDSigAndDEpsP

   Subroutine CalculatePrincipalStresses(IntGlo,Sig,SigPrin)
      !**********************************************************************
      !
      ! Implemented in the frame of the MPM project.
      !
      !**********************************************************************

      implicit none

      !Local variables
      double precision, dimension(3) :: xN1,xN2,xN3
      double precision :: Sig1,Sig2,Sig3,p,q
      !In Variables
      integer, intent(in) :: IntGlo ! Global ID of Gauss point or particle
      double precision, intent(in), dimension(6) :: Sig
      !Out Variables
      double precision, intent(out), dimension(6) :: SigPrin

      call PrincipalSig(1,Sig,xN1,xN2,xN3,Sig1,Sig2,Sig3,P,Q)

      If (Sig1 >= Sig2.and.Sig2 >= Sig3) then
         SigPrin(1) = Sig1
         SigPrin(2) = Sig2
         SigPrin(3) = Sig3
      else if (Sig1 >= Sig3.and.Sig3 >= Sig2) then
         SigPrin(1) = Sig1
         SigPrin(2) = Sig3
         SigPrin(3) = Sig2
      else if (Sig3 >= Sig1.and.Sig1 >= Sig2) then
         SigPrin(1) = Sig3
         SigPrin(2) = Sig1
         SigPrin(3) = Sig2
      else if (Sig3 >= Sig2.and.Sig2 >= Sig1) then
         SigPrin(1) = Sig3
         SigPrin(2) = Sig2
         SigPrin(3) = Sig1
      else if (Sig2 >= Sig1.and.Sig1 >= Sig3) then
         SigPrin(1) = Sig2
         SigPrin(2) = Sig1
         SigPrin(3) = Sig3
      else if (Sig2 >= Sig3.and.Sig3 >= Sig1) then
         SigPrin(1) = Sig2
         SigPrin(2) = Sig3
         SigPrin(3) = Sig1
      end if

      SigPrin(4) = 0.0d0
      SigPrin(5) = 0.0d0
      SigPrin(6) = 0.0d0

   end subroutine CalculatePrincipalStresses

   Subroutine PrincipalSig(IOpt, S, xN1, xN2, xN3, S1, S2, S3, P, Q)
      implicit none

      !! TODO: Need to assign intents for each variable

      integer :: iOPt
      ! Pricipal stress values
      double precision :: s1, s2, s3
      ! Stress invariants
      double precision :: P, Q
      double precision ::  S(:), xN1(:), xN2(:), xN3(:)

      If (iOpt.Eq.1) Then
         ! Calculates Eigenvalues and eigenvectors
         Call Eig_3_MohrCoulombStrainSoftening(0, S, xN1, xN2, xN3, S1, S2, S3, P, Q)
      Else

         Call Eig_3a_MohrCoulombStrainSoftening(0, S, S1, S2, S3, P, Q) ! no Eigenvectors
      End If
      Return
   end subroutine PrincipalSig

   Subroutine Eig_3_MohrCoulombStrainSoftening(iOpt,St,xN1,xN2,xN3,S1,S2,S3,P,Q)
      implicit none
      !TODO: Remove iO

      integer :: iOpt, i, j, k, it, itmax, is3, is2, is1, iq, ip, sign_tau
      ! ^ Optional Integer value (Not used)
      double precision :: v3p, v2p, v1p, t, tol, tau, s, s1, s2, s3, c, ap1, ap2, ap3, abs_max_s, a1p, a2p, a3p
      ! Stress invariants p, q
      double precision :: P, Q
      double precision :: St(6),A(3,3),V(3,3),xN1(3),xN2(3),xN3(3)
      !     *          xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues/Eigenvectors for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      ! PGB : adaption to Principal stress calculation
      !
      ! Applied on principal stresses, directions
      ! Stress vector St(): XX, YY, ZZ, XY, YZ, ZX
      !
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      ! Set V to unity matrix
      V(1,1) = 1
      V(2,1) = 0
      V(3,1) = 0

      V(1,2) = 0
      V(2,2) = 1
      V(3,2) = 0

      V(1,3) = 0
      V(2,3) = 0
      V(3,3) = 1


      abs_max_s=0.0
      Do i=1,3
         Do j=1,3
            if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
         End Do
      End Do
      Tol = 1d-20 * abs_max_s
      it = 0
      itmax = 50
      Do While ( it.Lt.itMax .And. abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )
         it=it+1
         Do k=1,3
            If (k .Eq. 1) Then
               ip=1
               iq=2
            Else If (k .Eq.2) Then
               ip=2
               iq=3
            Else
               ip=1
               iq=3
            End If
            If (a(ip,iq) .Ne. 0.0) Then
               tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
               If (tau .Ge.0.0) Then
                  sign_tau=1.0
               Else
                  sign_tau=-1.0
               End If
               t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
               c=1.0/sqrt(1.0+t*t)
               s=t*c
               a1p=c*a(1,ip)-s*a(1,iq)
               a2p=c*a(2,ip)-s*a(2,iq)
               a3p=c*a(3,ip)-s*a(3,iq)
               a(1,iq)=s*a(1,ip)+c*a(1,iq)
               a(2,iq)=s*a(2,ip)+c*a(2,iq)
               a(3,iq)=s*a(3,ip)+c*a(3,iq)
               a(1,ip)=a1p
               a(2,ip)=a2p
               a(3,ip)=a3p

               v1p=c*v(1,ip)-s*v(1,iq)
               v2p=c*v(2,ip)-s*v(2,iq)
               v3p=c*v(3,ip)-s*v(3,iq)
               v(1,iq)=s*v(1,ip)+c*v(1,iq)
               v(2,iq)=s*v(2,ip)+c*v(2,iq)
               v(3,iq)=s*v(3,ip)+c*v(3,iq)
               v(1,ip)=v1p
               v(2,ip)=v2p
               v(3,ip)=v3p

               ap1=c*a(ip,1)-s*a(iq,1)
               ap2=c*a(ip,2)-s*a(iq,2)
               ap3=c*a(ip,3)-s*a(iq,3)
               a(iq,1)=s*a(ip,1)+c*a(iq,1)
               a(iq,2)=s*a(ip,2)+c*a(iq,2)
               a(iq,3)=s*a(ip,3)+c*a(iq,3)
               a(ip,1)=ap1
               a(ip,2)=ap2
               a(ip,3)=ap3
            End If ! a(ip,iq)<>0
         End Do ! k
      End Do ! While
      ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)
      ! Derived invariants
      P = (S1+S2+S3)/3
      Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

      ! Sort eigenvalues S1 <= S2 <= S3
      is1 = 1
      is2 = 2
      is3 = 3
      if (s1.Gt.s2) Then
         t   = s2
         s2  = s1
         s1  = t
         it  = is2
         is2 = is1
         is1 = it
      End If
      if (s2.Gt.s3) Then
         t   = s3
         s3  = s2
         s2  = t
         it  = is3
         is3 = is2
         is2 = it
      End If
      if (s1.Gt.s2) Then
         t   = s2
         s2  = s1
         s1  = t
         it  = is2
         is2 = is1
         is1 = it
      End If
      Do i=1,3
         xN1(i) = v(i,is1) ! first  column
         xN2(i) = v(i,is2) ! second column
         xN3(i) = v(i,is3) ! third  column
      End Do
      Return
   end Subroutine Eig_3_MohrCoulombStrainSoftening ! Eig_3

   Subroutine Eig_3a_MohrCoulombStrainSoftening(iOpt,St,S1,S2,S3,P,Q) ! xN1,xN2,xN3,
      implicit none

      integer :: iOpt

      ! Principal stress values
      double precision :: s1, s2, s3

      ! Stress invariants
      double precision :: P, Q

      double precision :: St(6),A(3,3)   !  V(3,3),xN1(3),xN2(3),xN3(3)

      ! Local Variables
      integer :: i, j, k, it, itmax, is3, is2, is1, iq, ip, sign_tau
      double precision :: v3p, v2p, v1p, t, tol, tau, s, c, ap1, ap2, ap3, abs_max_s, a1p, a2p, a3p

      ! Get Eigenvalues ( no Eigenvectors) for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      !
      ! Applied on principal stresses, directions
      ! Stress vector XX, YY, ZZ, XY, YZ, ZX

      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      abs_max_s=0.0
      Do i=1,3
         Do j=1,3
            if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
         End Do
      End Do
      Tol = 1d-20 * abs_max_s
      If (iOpt.Eq.1) Tol = 1d-50*abs_max_s
      it=0
      itmax = 50
      Do While ( it.lt.itmax .And.&
         abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )

         it=it+1
         Do k=1,3
            If (k .Eq. 1) Then
               ip=1
               iq=2
            Else If (k .Eq.2) Then
               ip=2
               iq=3
            Else
               ip=1
               iq=3
            End If
            If (a(ip,iq) .Ne. 0.0) Then         ! ongelijk nul ?
               tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
               If (tau .Ge.0.0) Then
                  sign_tau=1.0
               Else
                  sign_tau=-1.0
               End If
               t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
               c=1.0/sqrt(1.0+t*t)
               s=t*c
               a1p=c*a(1,ip)-s*a(1,iq)
               a2p=c*a(2,ip)-s*a(2,iq)
               a3p=c*a(3,ip)-s*a(3,iq)
               a(1,iq)=s*a(1,ip)+c*a(1,iq)
               a(2,iq)=s*a(2,ip)+c*a(2,iq)
               a(3,iq)=s*a(3,ip)+c*a(3,iq)
               a(1,ip)=a1p
               a(2,ip)=a2p
               a(3,ip)=a3p

               ap1=c*a(ip,1)-s*a(iq,1)
               ap2=c*a(ip,2)-s*a(iq,2)
               ap3=c*a(ip,3)-s*a(iq,3)
               a(iq,1)=s*a(ip,1)+c*a(iq,1)
               a(iq,2)=s*a(ip,2)+c*a(iq,2)
               a(iq,3)=s*a(ip,3)+c*a(iq,3)
               a(ip,1)=ap1
               a(ip,2)=ap2
               a(ip,3)=ap3
            End If ! a(ip,iq)<>0
         End Do ! k
      End Do ! While

      ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)

      ! Derived invariants
      P = (S1+S2+S3)/3
      Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

      if (s1.Gt.s2) Then
         t   = s2
         s2  = s1
         s1  = t
      End If
      if (s2.Gt.s3) Then
         t   = s3
         s3  = s2
         s2  = t
      End If
      if (s1.Gt.s2) Then
         t   = s2
         s2  = s1
         s1  = t
      End If
      Return
   end Subroutine Eig_3a_MohrCoulombStrainSoftening
end module mod_MCSS_funcs
