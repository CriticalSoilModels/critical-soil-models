! Functions for calculating and updating the state parameters

module mod_state_params
   use kind_precision_module, only : real_type => dp
   use mod_yield_function, only: Get_dF_to_dSigma
   use mod_voigt_functions, only: TwoNormTensor, TensorInnerProduct
   implicit none

contains
   pure subroutine Check_Unloading(M_tc, eta_y, eta_yu, dI, Sig, dSig, LTOL, IsUnloading)
      !*********************************************************************
      ! Returns true if stress path is viscoplastic unloading              *
      !																	 *
      !*********************************************************************
      
      !input
      implicit none
      real(kind = real_type), intent(in):: M_tc, eta_y, eta_yu, dI, Sig(6), dSig(6), LTOL
      logical, intent(out):: IsUnloading

      !local variables
      real(kind = real_type):: deta, n_vec(6), n_norm, Sig_norm, dSIg_inner_n, beta, phi
      
      ! Init the unloading to false
      IsUnloading=.false.

      deta=eta_yu-eta_y!change in state parameter

      call Get_dF_to_dSigma(M_tc, eta_yu, Sig, n_vec)!Normal to surface
      call TwoNormTensor(n_vec, 6, n_norm) !norm of n_vec
      call TwoNormTensor(dSig, 6, Sig_norm) !norm of dSig
      call TensorInnerProduct(dSig, n_vec, 6,dSIg_inner_n) !inner product

      beta=acos(deta/(n_norm*Sig_norm))!conical aperture is a plane for inviscid mat.
      
      phi=acos(dSIg_inner_n/(n_norm*Sig_norm))!angle between stress rate and normal

      if (phi-beta>LTOL) IsUnloading=.true. !condition for unloading

   end subroutine Check_Unloading

   pure subroutine Get_I_coeff(D_part, G_s, p, eps_rate, I)
      !*********************************************************************
      ! Returns the inertial coefficient                                   *
      !																	 *
      !*********************************************************************
      implicit none
      !input
      real(kind = real_type), intent(in):: D_part, G_s, p, eps_rate
      !output
      real(kind = real_type), intent(out):: I

      ! TODO
      ! Currently this function is wrong. The assumption that is made is that 
      ! D_part  : [mm]
      ! eps_rate: [1/s]
      ! p'      : [kPa]
      ! G_{s}   : [-]
      ! These units don't work out when calculating I
      ! I should be unitless
      I = D_part * eps_rate * sqrt( G_s / abs(p) )

      ! The density should be used instead of the specific weight meaning
      ! D_part  : [m]
      ! eps_rate: [1/s]
      ! p'      : [kPa]
      ! \rho_{s}: [kg/m^{3}]

   end subroutine Get_I_coeff

   Subroutine check4crossing(IErate0I, IErateI, dErate_eff,RateRef, Apply)
      !******************************************************************
      ! determines if strain rate updating must occur                   *
      ! Boolean                                                         *
      !******************************************************************
      ! IErate0I: The previous (initial reference strain rate. State parameter pulled from the previous time step
      ! IErateI: The new inertial coefficient for this strain rate
      ! dErate_eff: The increment of strain rate change (This is calculated in this subroutine)
      ! RateRef: Reference strain rate other variables are compared to
      ! Apply: Boolean keeping track to determine if strain rate updates should be applied
      implicit none
      real(kind = real_type), intent(inout):: IErate0I, IErateI, dErate_eff
      real(kind = real_type), intent(in)   :: RateRef
      logical:: cond1, cond2, cond3
      logical, intent(out)::Apply
      Apply=.false.
      ! If the rate from the last time step is less than or equalt ot the reference rate, update the previous time step value to be the reference rate
      if(IErate0I<=RateRef) IErate0I=RateRef

      ! If the current rate is less than the reference rate than update the current rate to be the reference rate
      if (IErateI<=RateRef) IErateI=RateRef

      ! Cond1 - Checks if the rate has moved from slower than reference to faster than reference on this time step (Rate increased)
      cond1=(IErate0I==RateRef).and.(IErateI>RateRef)

      ! Cond2 - Checks if the rate has moved from greater than the reference to slower than the reference (Rate slowed)
      cond2=(IErate0I>RateRef).and.(IErateI==RateRef)

      ! Calculate the rate increment
      dErate_eff=IErateI-IErate0I

      ! Cond3 - Check if the current and previous value is greater than the reference rate. if they are that means that rate affects should be applied
      cond3=(IErate0I>RateRef).and.(IErateI>RateRef)

      ! Check if any of the conditions are true, if so strain rate effects need to be applied
      if (cond1.or.cond2.or.cond3) Apply=.true.
   end Subroutine check4crossing

   pure subroutine Get_M(M_tc, theta, M)
      !*********************************************************************
      ! Returns M															 *
      !																	 *
      !*********************************************************************
      implicit none
      !in
      real(kind = real_type), intent(in):: M_tc, theta
      !out
      real(kind = real_type), intent(out):: M
      !local
      real(kind = real_type):: COS_VAL
      real(kind = real_type), parameter :: pi=2*acos(0.0d0)

      COS_VAL=cos(1.5*theta+0.25*pi)

      M=M_tc*(1+0.25*COS_VAL**1.2)

   end subroutine Get_M

   pure subroutine Get_Dp(h, D_min, I, I_0, eps_q, k, ApplyRateUpdating, D)
      !*********************************************************************
      ! Returns the dilation for current inertial coefficient and dev.     *
      ! strain															 *
      !*********************************************************************
      implicit none
      logical, intent(in):: ApplyRateUpdating
      real(kind = real_type), intent(in):: h, D_min, I, I_0, eps_q, k
      real(kind = real_type), intent(out):: D

      !local variables
      real(kind = real_type):: D_mm

      if (ApplyRateUpdating) then
         ! D_mm=D_min*(I/I_0)**k !strain/rate hardening
         D_mm = update_strain_depend_param(D_min, I, I_0, k)
      else
         D_mm = D_min
      endif

      D = h * D_mm * eps_q * exp( 1.0_real_type - h * eps_q) !hardening rule
   end subroutine Get_Dp

   pure subroutine Update_GK(G_0, nu, I, I_0, k_G, k_K, G, K)
      !*********************************************************************
      ! Returns updated elastic modulus                                    *
      !																	 *
      !*********************************************************************
      implicit none
      real(kind = real_type), intent(in):: G_0, nu, I, I_0, k_G, k_K
      real(kind = real_type), intent(out):: G, K

      !local variables
      real(kind = real_type):: K_0

      ! G = G_0*(I/I_0)**k_G! updated modulus
      G = update_strain_depend_param(G_0, I, I_0, k_G)
      
      ! Update the bulk modulus
      K_0 = calc_bulk_modulus(G, nu)

      ! K = K_0*(I/I_0)**k_K! updated modulus
      K = update_strain_depend_param(K_0, I, I_0, k_K)
   end subroutine Update_GK

   pure function update_strain_depend_param(param, I, I_0, factor) result(updated_param)
      ! Function updates parameters that are dependent parameter
      ! This is for the bulk modulus and the shear modulus updates

      real(kind = real_type), intent(in) :: param, I, I_0, factor
      real(kind = real_type) :: updated_param

      updated_param = param * (I/I_0)**factor
   end function update_strain_depend_param

   pure function calc_bulk_modulus(shear_modulus, poisson_ratio) result(bulk_modulus)
      ! Calc the bulk modulus from the shear modulus and the poisson ratio
      real(kind = real_type), intent(in) :: shear_modulus, poisson_ratio
      real(kind = real_type) :: bulk_modulus

      ! Local variables
      real(kind = real_type), parameter :: ONE = 1.0_real_type, &
                                           TWO = 2.0_real_type, &
                                           THREE = 3.0_real_type 

      bulk_modulus = TWO * shear_modulus * (ONE + poisson_ratio) / ( THREE * (ONE - TWO * poisson_ratio) )
   end function calc_bulk_modulus

end module mod_state_params
