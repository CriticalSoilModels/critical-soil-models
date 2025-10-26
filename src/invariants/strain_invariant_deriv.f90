! Module for holding the derivatives of the strain invariants

module mod_strain_invar_deriv
   use kind_precision_module, only: real_type => dp
   use mod_strain_invariants, only: calc_dev_strain, calc_eps_vol_invariant
   implicit none

contains
   subroutine Get_dEpsq_to_dEps(Epsq, Eps, dEqdEpsq)
      !************************************************************************
      ! Returns the derivative of the deviatoric strain with respect to the   *
      ! deviatoric strain	tensor						     					*
      ! dEqdEpsq is a (1X6) vector											*
      !************************************************************************
      implicit none
      !input
      real(kind = real_type), intent(in):: Epsq, Eps(6)
      !output
      real(kind = real_type), intent(out):: dEqdEpsq(6)
      !local variables
      real(kind = real_type):: evol, dev(6)

      evol=calc_eps_vol_invariant(Eps)

      dev=calc_dev_strain(Eps, evol)

      if (Epsq>0.0d0) then !in case of zero plastic strain
         dEqdEpsq=(2.0/(3.0*Epsq))*dev
      else
         dEqdEpsq=0.0d0
      endif
   end subroutine Get_dEpsq_to_dEps

   pure function calc_inc_driver_dEpsq_to_dEps(Eps) result(dEq_dEps)
      ! Function calculates the dEpsq/dEps for incremental driver voigt ordering
      ! This is just to check the general version of the dEpsq/dEps function
      ! The values for this derivative come from python using sympy in Epsq_derivative_check.ipynb
   
      real(kind = real_type), intent(in) :: Eps(6)
      real(kind = real_type) :: dEq_dEps(6)

      ! Local variables
      real(kind = real_type), parameter :: TWO = 2.0_real_type, &
         THREE = 3.0_real_type, &
         FOUR = 4.0_real_type

      dEq_dEps(1) = FOUR*Eps(1)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)) &
         - TWO*Eps(2)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)) &
         - TWO*Eps(3)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2))


      dEq_dEps(2) = -TWO*Eps(1)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)) &
         + FOUR*Eps(2)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)) &
         - TWO*Eps(3)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) &
         + FOUR*Eps(2)**2 - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2))


      dEq_dEps(3) = -TWO*Eps(1)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)) &
         - TWO*Eps(2)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)) &
         + FOUR*Eps(3)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2))


      dEq_dEps(4) = TWO*Eps(4)/sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) &
         + FOUR*Eps(2)**2 - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 &
         + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)


      dEq_dEps(5) = TWO*Eps(5)/sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) &
         + FOUR*Eps(2)**2 - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 &
         + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)


      dEq_dEps(6) = TWO*Eps(6)/sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) &
         + FOUR*Eps(2)**2 - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 &
         + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)


   end function calc_inc_driver_dEpsq_to_dEps
end module mod_strain_invar_deriv
