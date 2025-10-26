module mod_SRMC_Ortiz_Simo

   use mod_SRMC_funcs        , only: MatVec, DotProduct_2 ! Subroutines that are required for calculations
   use mod_state_params      , only: check4crossing, Update_GK, Get_Dp, Get_M
   use mod_state_params_deriv, only: Get_dD_to_dEpsP
   use mod_stress_invariants , only: Get_invariants
   use mod_strain_invariants , only: Get_strain_invariants
   use mod_yield_function    , only: YieldFunction, Get_dF_to_dSigma
   use mod_plastic_potential , only: Get_dP_to_dSigma
   use mod_stress_invariants , only : Get_invariants

   implicit none

contains
   subroutine Ortiz_Simo_Integration(G_0, nu, M_tc, M, No, D_min, h, G, K, eta_y, Dp, &
      I_0, I, dI, k_G, k_K, k_D, Sig, EpsP, dEps, &
      FTOL, NOEL, max_stress_iters)

      !------------------------Function Information------------------------
      ! Ortiz and Simo (1986) integration scheme (Algorithm on pg. 6)
      ! Returns:
      ! Increment of plastic strain, plastic strain, and change of state parameters
      !----------------------End Function Information----------------------

      implicit none
      ! Variable Info
      !!!------Begin External variables----!!!
      ! G_0: Initial Shear Modulus
      ! nu: Poisson ratio
      ! M_tc: Critical stress ratio for triaxial loading
      ! M: Current critical stress ratio
      ! No: Nova's volumetric coupling coefficient
      ! D_min: Minimum dilatancy
      ! h: Dilatancy hardening parameter
      ! G: Current Shear modulus
      ! K: Current Bulk Modulus
      ! eta_y: Current Stress ratio
      ! Dp: Current dilatancy
      ! I_0: Reference interial coefficient
      ! I: Current inertial coefficient
      ! dI: Increment of the inertial coeff
      ! k_G: Shear modulus viscosity coeff
      ! k_K: Bulk modulus viscosity coeff
      ! k_D: Dilatancy viscosity coeff
      ! Sig: Input stress
      ! EpsP: Plastic strain
      ! dEps: Increment of strain
      ! dEpsp: Increment of Plastic Strain
      ! FTOL: Tolerance for initial predictor
      !!!------End External variables-----!!!

      !!!------Begin new output variables------!!!
      !!!---------End Returned calculated variables-----!!!

      !--------------Input variables--------------!

      ! Input scalar values
      double precision, intent(in):: G_0, nu, M_tc, No, D_min, h, &
         k_G, k_K, k_D, FTOL
      integer, intent(in) :: NOEL, max_stress_iters
      ! Input vector values
      double precision, dimension (6), intent(in):: dEps

      !-------------End Input Variables-----------!

      !--------------Output Variables-------------!
      ! In/Out scalar values
      double precision, intent(inout):: G, K, eta_y, Dp, I_0, I, dI, M

      ! In/Out Vector values
      double precision, dimension(6), intent(inout):: Sig, EpsP !, dEpsP

      ! Out scalar values
      !double precision, intent(out)::

      ! Out vector values
      ! double precision, dimension(6), intent(out):: dSig

      !-------------End Output Variables--------!

      !-------------local Variables-------------!
      ! Local scalar values
      double precision:: I_f, F, p, q, epsv_p, epsq_p, eta_yu, Du, Mu, dummyVal, a_Dot_m, L, H_term, &
         D1, D2, b, dLambda, dD

      logical:: ApplyStrainRateUpdate = .false.
      integer:: counter

      ! Local vector values
      double precision, dimension(6):: dEpsE, dummyVec, dEpsPu, EpsPu, Sigu, &
         m_vec, n_vec, DE_m, a, dSig

      ! Local matrix values
      double precision, dimension(6,6):: DE

      !------------End local Variables----------!

      ! Store variables for updating
      Sigu = Sig
      EpsPu = EpsP
      ! dEpsPu = dEpsP
      eta_yu = eta_y
      Du = Dp
      Mu = M

      ! Apply strain rate updates
      I_f=I+dI
      call check4crossing(I,  I_f, dI, I_0, ApplyStrainRateUpdate)
      call Get_strain_invariants(EpsPu, epsv_p, epsq_p)
      call Update_GK(G_0, nu, I_f, I_0, k_G, k_K, G, K)
      call Get_Dp(h, D_min, I_f, I_0, epsq_p, k_D, ApplyStrainRateUpdate, Du) ! Need this for the strain rate but also to actually calculate Du
      eta_yu = Mu-du*(1.0 * No)

      ! Turn off strain rate affects
      ApplyStrainRateUpdate = .false.

      !--------------------Compute elastic predictor---------------------------!
      dEpsE = dEps ! This is Luis's option for predicting the elastic strain

      !-----------Ensemble elastic matrix------------------!
      D1  = K+(4*G/3)
      D2  = K-(2*G/3)
      DE  = 0.0
      DE(1:3,1:3) = D2
      DE(1,1) = D1
      DE(2,2) = D1
      DE(3,3) = D1
      DE(4,4) = G
      DE(5,5) = G
      DE(6,6) = G
      !----------End Ensemble elastic matrix--------------!

      ! Calc the stress predictor
      call MatVec(DE, 6, dEpsE, 6, dSig)

      !Update the stresses
      Sigu = Sigu + dSig
      !-----------------End Compute elastic predictor------------------------!
      ! Since all strain is assumed to be elastic no plastic strain increment

      !-------------------Begin Yielding Check--------------------------!

      ! Compute stress invariants
      call Get_invariants(Sigu, p, q, dummyVal)

      ! M = M_tc*(1 + 0.25(cos(1.5 * theta + 0.25 *pi))
      call Get_M(M_tc, dummyVal, Mu)

      eta_yu = Mu - Du * (1.0 -No)

      ! Compute the value of the yield Function
      call YieldFunction(q, p, eta_yu, F)

      if (F < FTOL) then
         ! Prediction is correct, stress and strain values can be updated and returned

         ! Update Sig, EpsP, dEpsP, eta_y, Dp
         Sig = Sigu
         EpsP = EpsPu
         ! dEpsP = dEpsPu
         eta_y = eta_yu
         Dp = Du
         ! Exit out of the Subroutine and return values
         return
      endif
      ! else
      ! Predicted stress too large, iterations are needed
      !-------------------End Yielding Check--------------------------!


      !-----------------------Begin Plastic Descent-----------------------!


      counter = 0

      call Get_dD_to_dEpsP(D_min, h, I_0, k_D, epsq_p, epsv_p, &
         dEpsPu, I, ApplyStrainRateUpdate, a) !a=dD/dEpsP

      ! Compute stress invariants
      call Get_invariants(Sigu, p, q, dummyVal)

      do while (abs(F) >= FTOL .and. counter <= max_stress_iters)
         !---------------------Begin Compute derivatives--------------------------!
         call Get_dF_to_dSigma(Mu, eta_yu, Sigu, n_vec) !n=dF/dSig
         call Get_dP_to_dSigma(Du, Sigu, m_vec) !m=dP/dSig
         L = -p * (1-No) !dF/Xs = dF/dDp = Xi in Ortiz & Simo

         !-----------------------End Compute derivatives--------------------------!

         !----------------Compute Denominator of dLambda--------------------------!
         ! Denominator = v:D:r - Xi.h == n_vec:DE:m_vec - H_term

         ! Compute DE.m_vec
         call MatVec(DE, 6, m_vec, 6,  DE_m)

         ! Compute n_vec.DE.m_vec
         call DotProduct_2(n_vec, DE_m, 6, dummyVal)

         ! compute Xi.h = H_term = L.a.m
         call DotProduct_2(a, m_vec, 6, a_Dot_m)
         H_term = L * a_Dot_m

         !--------------End Compute Denominator of dLambda-----------------------!

         ! Compute the change in the plastic potential (lambda)
         dLambda = F/(dummyVal - H_term)

         ! Compute the stress update
         Sigu = Sigu - dLambda * DE_m

         ! Compute stress invariants
         call Get_invariants(Sigu, p, q, dummyVal)

         ! Update M
         call Get_M(M_tc, dummyVal, Mu)

         ! calc dLambda * m (increment of plastic strain)
         dEpsPu = dLambda * m_vec

         ! Accumulate plastic strain
         EpsPu = EpsPu + dEpsPu

         ! Calc strain invariants
         ! call Get_strain_invariants(EpsPu, epsv_p, epsq_p)

         ! ! Calc a = dD/dEpsP
         ! call Get_dD_to_dEpsP(D_min, h, I_0, k_D, epsq_p, epsv_p, &
         !    dEpsPu, I, ApplyStrainRateUpdate, a) !a=dD/dEpsP

         ! Don't need to update the diltancy. It's a function of the total deviatoric strain. Not the plastic part 

         ! Update eta_y
         eta_yu = Mu - Du * (1.0 -No)

         ! Calc the yield function
         call YieldFunction(q, p, eta_yu, F)

         ! Update Counter to end while loop
         Counter = Counter + 1

      end do

      ! Calc the change in Sig
      dSig = Sigu - Sig ! Don't think I need this because I'm not substepping

      !Calc dEpsP
      dEpsPu = EpsPu - EpsP

      ! Update Sig, dEpsP, EpsP, eta_y, Dp
      Sig = Sigu
      ! dEpsP = dEpsPu
      EpsP = EpsPu
      eta_y = eta_yu

      dD = Du -Dp
      Dp = Du
      M = Mu

   end subroutine Ortiz_Simo_Integration

end module mod_SRMC_Ortiz_Simo
