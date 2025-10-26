! Functions that calculate the derivatives of the state paramters

module mod_state_params_deriv
   use kind_precision_module, only: real_type => dp

   implicit none

contains

   pure subroutine Get_dD_to_dI(D_min, h, I_0, kD, eps_q, I, b)
      !************************************************************************
      ! Returns the derivative of the Dilation with respect to the inertial	*
      ! coefficient 															*
      ! b=dD/dI																*
      ! b is a scalar															*
      !************************************************************************
      implicit none
      real(kind = real_type), intent(in):: D_min, h, I_0, kD, eps_q, I
      !output
      real(kind = real_type), intent(out)::b
      !local variables

      b=h*D_min*eps_q*exp(1-h*eps_q)*kD*((I/I_0)**(kD-1.0))/I_0

   end subroutine Get_dD_to_dI

   subroutine Get_dD_to_dEpsP(D_min, h, I_0, k_D, epsq_p, epsv_p, &
      EpsP, I, ApplyStrainRateUpdate, a)
      !************************************************************************
      ! Returns the derivative of the Dilation with respect to the plastic    *
      ! strain																*
      ! a=dXs/dEpsp= dD/dEpsq * dEPsq/dEpsp									*
      ! a is a (1X6) vector													*
      !************************************************************************
      implicit none
      !input
      logical, intent(in):: ApplyStrainRateUpdate
      real(kind = real_type), intent(in):: D_min, h, I_0, k_D, epsq_p, epsv_p, &
         EpsP(6), I
      !output
      real(kind = real_type), intent(out):: a(6)
      !local variables
      real(kind = real_type):: D, dDdEpsq_p, dev(6),dEpsq_pdEpsp(6)

      !________________________________________________________________________
      !Get dD/dEpsq_p
      if (ApplyStrainRateUpdate) then
         D=D_min*(I/I_0)**k_D
      else
         D=D_min
      end if

      dDdEpsq_p=h*D*exp(1.0-h*epsq_p)*(1.0-h*epsq_p)
      !_______________________________________________________________________

      !_______________________________________________________________________
      !Get dEpsp_Q/dEpsp=
      dev=EpsP
      dev(1)=dev(1)-(epsv_p/3.0)
      dev(2)=dev(2)-(epsv_p/3.0)
      dev(3)=dev(3)-(epsv_p/3.0) !deviatoric stress tensor

      if (epsq_p>0.0d0) then !in case of zero plastic strain
         dEpsq_pdEpsp=(2.0/(3.0*epsq_p))*dev
      else
         dEpsq_pdEpsp=0.0d0
      endif
      !_______________________________________________________________________

      !_______________________________________________________________________
      !Get a=dXs/dEpsp

      a=dDdEpsq_p*dEpsq_pdEpsp
      !______________________________________________________________________
   end subroutine Get_dD_to_dEpsP


end module mod_state_params_deriv
