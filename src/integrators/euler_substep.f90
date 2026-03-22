! =============================================================================
! PSEUDOCODE — does not compile
! Generic modified Euler integrator with automatic substepping and error control.
! Knows nothing about any specific model — calls only the five deferred procedures
! defined in csm_model_t.
! =============================================================================

module mod_euler_substep
   use iso_fortran_env, only: dp => real64
   use mod_csm_model
   implicit none

   real(dp), parameter :: DT_MIN_DEFAULT = 1.0e-9_dp

contains

   subroutine euler_substep(model, stress, dstrain, ftol, stol)
      class(csm_model_t), intent(inout) :: model
      real(dp),           intent(inout) :: stress(6)   ! updated in place
      real(dp),           intent(in)    :: dstrain(6)  ! strain increment (internal convention)
      real(dp),           intent(in)    :: ftol         ! yield surface tolerance
      real(dp),           intent(in)    :: stol         ! relative stress error tolerance

      real(dp) :: D(6,6)
      real(dp) :: sig_tr(6), F
      real(dp) :: alpha
      real(dp) :: T, DT, R
      real(dp) :: sig_t(6), dstrain_plastic(6)
      real(dp) :: dSig1(6), dSig2(6), deps_p1(6), deps_p2(6)
      real(dp), allocatable :: saved(:)
      logical  :: step_failed

      ! -----------------------------------------------------------------------
      ! 1. Elastic predictor
      ! -----------------------------------------------------------------------
      D      = model%elastic_stiffness()
      sig_tr = stress + matmul(D, dstrain)

      ! -----------------------------------------------------------------------
      ! 2. Check yield — purely elastic step
      ! -----------------------------------------------------------------------
      F = model%yield_fn(sig_tr)
      if (F < ftol) then
         stress = sig_tr
         return
      end if

      ! -----------------------------------------------------------------------
      ! 3. Find elastic/plastic split alpha in [0,1]
      !    stress + alpha * D * dstrain lies on the yield surface
      !    Uses only model%yield_fn — no model-specific knowledge
      ! -----------------------------------------------------------------------
      call find_yield_intersection(model, stress, dstrain, D, ftol, alpha)

      ! Stress on yield surface + remaining plastic strain increment
      stress          = stress + alpha * matmul(D, dstrain)
      dstrain_plastic = (1.0_dp - alpha) * dstrain

      ! -----------------------------------------------------------------------
      ! 4. Modified Euler substep loop with error control
      ! -----------------------------------------------------------------------
      T    = 0.0_dp
      DT   = 1.0_dp

      do while (T < 1.0_dp)

         sig_t = stress

         ! --- First Euler estimate ---
         call euler_step(model, sig_t, dstrain_plastic*DT, D, dSig1, deps_p1)

         ! --- Second Euler estimate (from end of first step) ---
         ! Snapshot state, advance to end of first step, compute second estimate,
         ! then restore. No index magic here — each model handles its own layout.
         call model%snapshot(saved)
         call model%update_hardening(deps_p1)
         call euler_step(model, sig_t + dSig1, dstrain_plastic*DT, D, dSig2, deps_p2)
         call model%restore(saved)

         ! --- Relative error estimate ---
         R = 0.5_dp * norm2(dSig1 - dSig2) &
             / max(norm2(sig_t + 0.5_dp*(dSig1 + dSig2)), 1.0e-12_dp)

         if (R > stol) then
            ! Step rejected — shrink DT and retry
            DT          = max(0.9_dp * sqrt(stol/R) * DT, DT_MIN_DEFAULT)
            step_failed = .true.

         else
            ! Step accepted — take modified Euler average
            stress = sig_t + 0.5_dp * (dSig1 + dSig2)
            call model%update_hardening(0.5_dp * (deps_p1 + deps_p2))

            ! Drift correction: push stress back onto yield surface
            F = model%yield_fn(stress)
            if (abs(F) > ftol) call correct_drift(model, stress, D, ftol)

            ! Advance pseudo-time, grow DT if step was easy
            T  = T + DT
            DT = min(0.9_dp * sqrt(stol / max(R, 1.0e-12_dp)) * DT, 1.0_dp - T)
            if (step_failed) then
               DT          = min(DT, 1.0_dp - T)   ! don't overshoot after a failure
               step_failed = .false.
            end if
         end if

      end do

   end subroutine euler_substep

   ! ---------------------------------------------------------------------------
   ! Single explicit Euler stress increment
   ! Returns dSig and deps_p for one sub-increment dstrain_sub
   ! ---------------------------------------------------------------------------
   subroutine euler_step(model, stress, dstrain_sub, D, dSig, deps_p)
      class(csm_model_t), intent(in)  :: model
      real(dp),           intent(in)  :: stress(6), dstrain_sub(6), D(6,6)
      real(dp),           intent(out) :: dSig(6), deps_p(6)

      real(dp) :: df(6), dg(6), dlambda
      real(dp) :: D_dg(6)      ! D * dg
      real(dp) :: df_D(6)      ! df^T * D (row vector)

      df  = model%flow_rule(stress)
      dg  = model%plastic_potential(stress)

      ! Plastic multiplier from consistency condition:
      !   dlambda = (df . D . dstrain) / (df . D . dg - df . dh/dlambda)
      ! (hardening denominator term omitted here — placeholder)
      D_dg   = matmul(D, dg)
      df_D   = matmul(df, D)     ! note: df is a row in the denominator
      dlambda = dot_product(df_D, dstrain_sub) / dot_product(df_D, dg)
      dlambda = max(dlambda, 0.0_dp)   ! no negative plastic multiplier

      deps_p = dlambda * dg
      dSig   = matmul(D, dstrain_sub - deps_p)
   end subroutine euler_step

   ! ---------------------------------------------------------------------------
   ! Pegasus method: find alpha such that yield_fn(stress + alpha*dSig_e) = 0
   ! Pure bisection fallback if Pegasus stalls.
   ! ---------------------------------------------------------------------------
   subroutine find_yield_intersection(model, stress, dstrain, D, ftol, alpha)
      class(csm_model_t), intent(in)  :: model
      real(dp),           intent(in)  :: stress(6), dstrain(6), D(6,6), ftol
      real(dp),           intent(out) :: alpha

      real(dp) :: dSig_e(6)
      real(dp) :: alpha0, alpha1, F0, F1, alpha_new, F_new
      integer  :: iter

      dSig_e = matmul(D, dstrain)
      alpha0 = 0.0_dp;  F0 = model%yield_fn(stress)
      alpha1 = 1.0_dp;  F1 = model%yield_fn(stress + dSig_e)

      ! If already outside at start, alpha = 0
      if (F0 >= -ftol) then
         alpha = 0.0_dp
         return
      end if

      ! Pegasus iteration
      do iter = 1, 50
         alpha_new = alpha1 - F1 * (alpha1 - alpha0) / (F1 - F0)
         F_new     = model%yield_fn(stress + alpha_new * dSig_e)
         if (abs(F_new) < ftol) exit
         if (F_new * F1 < 0.0_dp) then
            alpha0 = alpha1;  F0 = F1
         else
            F0 = F0 * F1 / (F1 + F_new)
         end if
         alpha1 = alpha_new;  F1 = F_new
      end do

      alpha = alpha_new
   end subroutine find_yield_intersection

   ! ---------------------------------------------------------------------------
   ! Drift correction: project stress back onto yield surface
   ! Simple normal correction — more sophisticated methods possible
   ! ---------------------------------------------------------------------------
   subroutine correct_drift(model, stress, D, ftol)
      class(csm_model_t), intent(inout) :: model
      real(dp),           intent(inout) :: stress(6)
      real(dp),           intent(in)    :: D(6,6), ftol

      real(dp) :: F, df(6), dg(6), dlambda
      integer  :: iter

      do iter = 1, 100
         F = model%yield_fn(stress)
         if (abs(F) < ftol) return

         df      = model%flow_rule(stress)
         dg      = model%plastic_potential(stress)
         dlambda = F / dot_product(df, matmul(D, dg))
         stress  = stress - dlambda * matmul(D, dg)
         call model%update_hardening(dlambda * dg)
      end do

   end subroutine correct_drift

end module mod_euler_substep
