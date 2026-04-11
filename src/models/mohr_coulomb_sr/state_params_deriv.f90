! Functions that calculate the derivatives of the state paramters

module mod_state_params_deriv
   use mod_csm_kinds, only: wp

   implicit none

contains

   pure subroutine calc_ddil_by_dI(D_min, h, I_0, kD, eps_q, I, b)
      !! Returns the derivative of the dilation with respect to the inertial
      !! coefficient: b = ddil/dI (scalar)
      implicit none
      real(wp), intent(in)  :: D_min, h, I_0, kD, eps_q, I
      real(wp), intent(out) :: b

      b = h*D_min*eps_q*exp(1-h*eps_q)*kD*((I/I_0)**(kD-1.0_wp))/I_0

   end subroutine calc_ddil_by_dI

   subroutine calc_ddil_by_deps_p(D_min, h, I_0, k_D, epsq_p, epsv_p, &
      eps_p, I, apply_strain_rate_update, a)
      !! Returns the derivative of the dilation with respect to plastic strain:
      !! a = ddil/dEpsq * dEpsq/dEpsp  (6-vector)
      implicit none
      logical,  intent(in)  :: apply_strain_rate_update
      real(wp), intent(in)  :: D_min, h, I_0, k_D, epsq_p, epsv_p, eps_p(6), I
      real(wp), intent(out) :: a(6)

      real(wp) :: dil, ddil_by_depsq_p, dev(6), deps_q_by_deps_p(6)

      if (apply_strain_rate_update) then
         dil = D_min*(I/I_0)**k_D
      else
         dil = D_min
      end if

      ddil_by_depsq_p = h*dil*exp(1.0_wp - h*epsq_p)*(1.0_wp - h*epsq_p)

      dev = eps_p
      dev(1) = dev(1) - (epsv_p/3.0_wp)
      dev(2) = dev(2) - (epsv_p/3.0_wp)
      dev(3) = dev(3) - (epsv_p/3.0_wp)

      if (epsq_p > 0.0_wp) then
         deps_q_by_deps_p = (2.0_wp/(3.0_wp*epsq_p))*dev
      else
         deps_q_by_deps_p = 0.0_wp
      end if

      a = ddil_by_depsq_p * deps_q_by_deps_p

   end subroutine calc_ddil_by_deps_p


end module mod_state_params_deriv
