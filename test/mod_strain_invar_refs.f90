! Reference implementations of strain invariant derivatives.
! These use the full Voigt-expanded sympy-generated form and exist solely
! for cross-checking the production implementations in tests.
! Do NOT use in production code.

module mod_strain_invar_refs
   use iso_fortran_env, only: dp => real64
   implicit none
   private
   public :: calc_deps_q_by_deps_full

contains

   ! deps_q/deps computed via the fully-expanded sympy formula.
   ! Used as an independent cross-check in tests.
   pure function calc_deps_q_by_deps_full(Eps) result(deps_q_by_deps)
      real(dp), intent(in) :: Eps(6)
      real(dp) :: deps_q_by_deps(6)

      ! Local variables
      real(dp), parameter :: TWO = 2.0_dp, &
         THREE = 3.0_dp, &
         FOUR = 4.0_dp

      deps_q_by_deps(1) = FOUR*Eps(1)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)) &
         - TWO*Eps(2)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)) &
         - TWO*Eps(3)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2))


      deps_q_by_deps(2) = -TWO*Eps(1)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)) &
         + FOUR*Eps(2)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)) &
         - TWO*Eps(3)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) &
         + FOUR*Eps(2)**2 - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2))


      deps_q_by_deps(3) = -TWO*Eps(1)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)) &
         - TWO*Eps(2)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)) &
         + FOUR*Eps(3)/(THREE*sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) + FOUR*Eps(2)**2 &
         - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2))


      deps_q_by_deps(4) = TWO*Eps(4)/sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) &
         + FOUR*Eps(2)**2 - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 &
         + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)


      deps_q_by_deps(5) = TWO*Eps(5)/sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) &
         + FOUR*Eps(2)**2 - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 &
         + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)


      deps_q_by_deps(6) = TWO*Eps(6)/sqrt(FOUR*Eps(1)**2 - FOUR*Eps(1)*Eps(2) - FOUR*Eps(1)*Eps(3) &
         + FOUR*Eps(2)**2 - FOUR*Eps(2)*Eps(3) + FOUR*Eps(3)**2 &
         + THREE*Eps(4)**2 + THREE*Eps(5)**2 + THREE*Eps(6)**2)


   end function calc_deps_q_by_deps_full

end module mod_strain_invar_refs
