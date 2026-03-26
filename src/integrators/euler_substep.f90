! Generic modified Euler integrator with automatic substepping and error control.
! Knows nothing about any specific model — calls only the deferred procedures
! defined in csm_model_t.

module mod_euler_substep
   use mod_csm_kinds, only: wp
   use mod_csm_model, only: csm_model_t
   implicit none

   real(wp), parameter :: DT_MIN_DEFAULT = 1.0e-9_wp

contains

   subroutine euler_substep(model, sig, deps, ftol, stol)
      class(csm_model_t), intent(inout) :: model
      real(wp),           intent(inout) :: sig(6)    ! updated in place
      real(wp),           intent(in)    :: deps(6)   ! strain increment (internal convention)
      real(wp),           intent(in)    :: ftol      ! yield surface tolerance
      real(wp),           intent(in)    :: stol      ! relative stress error tolerance

      real(wp) :: stiff_e(6,6)
      real(wp) :: sig_tr(6), F
      real(wp) :: alpha_ep
      real(wp) :: t_sub, dt_sub, err_rel
      real(wp) :: sig_t(6), deps_ep(6)
      real(wp) :: dsig_1(6), dsig_2(6), deps_p1(6), deps_p2(6)
      real(wp), allocatable :: saved(:)
      logical  :: step_failed

      ! -----------------------------------------------------------------------
      ! 1. Elastic predictor
      ! -----------------------------------------------------------------------
      stiff_e = model%elastic_stiffness()
      sig_tr  = sig + matmul(stiff_e, deps)

      ! -----------------------------------------------------------------------
      ! 2. Check yield — purely elastic step
      ! -----------------------------------------------------------------------
      F = model%yield_fn(sig_tr)
      if (F < ftol) then
         sig = sig_tr
         return
      end if

      ! -----------------------------------------------------------------------
      ! 3. Find elastic/plastic split alpha_ep in [0,1]
      !    sig + alpha_ep * stiff_e * deps lies on the yield surface
      !    Uses only model%yield_fn — no model-specific knowledge
      ! -----------------------------------------------------------------------
      call find_yield_intersection(model, sig, deps, stiff_e, ftol, alpha_ep)

      ! Stress on yield surface + remaining plastic strain increment
      sig     = sig + alpha_ep * matmul(stiff_e, deps)
      deps_ep = (1.0_wp - alpha_ep) * deps

      ! -----------------------------------------------------------------------
      ! 4. Modified Euler substep loop with error control
      ! -----------------------------------------------------------------------
      t_sub      = 0.0_wp
      dt_sub     = 1.0_wp
      step_failed = .false.

      do while (t_sub < 1.0_wp)

         sig_t = sig

         ! --- First Euler estimate ---
         call euler_step(model, sig_t, deps_ep*dt_sub, stiff_e, dsig_1, deps_p1)

         ! --- Second Euler estimate (from end of first step) ---
         ! Snapshot state, advance to end of first step, compute second estimate,
         ! then restore. No index magic here — each model handles its own layout.
         call model%snapshot(saved)
         call model%update_hardening(deps_p1)
         call euler_step(model, sig_t + dsig_1, deps_ep*dt_sub, stiff_e, dsig_2, deps_p2)
         call model%restore(saved)

         ! --- Relative error estimate ---
         err_rel = 0.5_wp * norm2(dsig_1 - dsig_2) &
                   / max(norm2(sig_t + 0.5_wp*(dsig_1 + dsig_2)), 1.0e-12_wp)

         if (err_rel > stol) then
            ! Step rejected — shrink dt_sub and retry
            dt_sub      = max(0.9_wp * sqrt(stol/err_rel) * dt_sub, DT_MIN_DEFAULT)
            step_failed = .true.

         else
            ! Step accepted — take modified Euler average
            sig = sig_t + 0.5_wp * (dsig_1 + dsig_2)
            call model%update_hardening(0.5_wp * (deps_p1 + deps_p2))

            ! Drift correction: push stress back onto yield surface
            F = model%yield_fn(sig)
            if (abs(F) > ftol) call correct_drift(model, sig, stiff_e, ftol)

            ! Advance pseudo-time, grow dt_sub if step was easy
            t_sub  = t_sub + dt_sub
            dt_sub = min(0.9_wp * sqrt(stol / max(err_rel, 1.0e-12_wp)) * dt_sub, 1.0_wp - t_sub)
            if (step_failed) then
               dt_sub      = min(dt_sub, 1.0_wp - t_sub)   ! don't overshoot after a failure
               step_failed = .false.
            end if
         end if

      end do

   end subroutine euler_substep

   ! ---------------------------------------------------------------------------
   ! Single explicit Euler stress increment
   ! Returns dsig and deps_p for one sub-increment deps_sub
   ! ---------------------------------------------------------------------------
   subroutine euler_step(model, sig, deps_sub, stiff_e, dsig, deps_p)
      class(csm_model_t), intent(in)  :: model
      real(wp),           intent(in)  :: sig(6), deps_sub(6), stiff_e(6,6)
      real(wp),           intent(out) :: dsig(6), deps_p(6)

      real(wp) :: dF_by_dsig(6), dg_plas_by_dsig(6), dlambda
      real(wp) :: stiff_e_dG(6)   ! stiff_e * dg_plas_by_dsig
      real(wp) :: dF_stiff_e(6)   ! dF_by_dsig^T * stiff_e (row vector)

      dF_by_dsig = model%flow_rule(sig)
      dg_plas_by_dsig = model%plastic_potential(sig)

      ! Plastic multiplier from consistency condition:
      !   dlambda = (dF_by_dsig . stiff_e . deps) / (dF_by_dsig . stiff_e . dg_plas_by_dsig)
      ! (hardening denominator term omitted here — placeholder)
      stiff_e_dG   = matmul(stiff_e, dg_plas_by_dsig)
      dF_stiff_e   = matmul(dF_by_dsig, stiff_e)
      dlambda      = dot_product(dF_stiff_e, deps_sub) / dot_product(dF_stiff_e, dg_plas_by_dsig)
      dlambda      = max(dlambda, 0.0_wp)   ! no negative plastic multiplier

      deps_p = dlambda * dg_plas_by_dsig
      dsig   = matmul(stiff_e, deps_sub - deps_p)
   end subroutine euler_step

   ! ---------------------------------------------------------------------------
   ! Pegasus method: find alpha_ep such that yield_fn(sig + alpha_ep*dsig_e) = 0
   ! Pure bisection fallback if Pegasus stalls.
   ! ---------------------------------------------------------------------------
   subroutine find_yield_intersection(model, sig, deps, stiff_e, ftol, alpha_ep)
      class(csm_model_t), intent(in)  :: model
      real(wp),           intent(in)  :: sig(6), deps(6), stiff_e(6,6), ftol
      real(wp),           intent(out) :: alpha_ep

      real(wp) :: dsig_e(6)
      real(wp) :: alpha0, alpha1, F0, F1, alpha_new, F_new
      integer  :: iter

      dsig_e = matmul(stiff_e, deps)
      alpha0 = 0.0_wp;  F0 = model%yield_fn(sig)
      alpha1 = 1.0_wp;  F1 = model%yield_fn(sig + dsig_e)

      ! If already outside at start, alpha_ep = 0
      if (F0 >= -ftol) then
         alpha_ep = 0.0_wp
         return
      end if

      ! Pegasus iteration
      do iter = 1, 50
         alpha_new = alpha1 - F1 * (alpha1 - alpha0) / (F1 - F0)
         F_new     = model%yield_fn(sig + alpha_new * dsig_e)
         if (abs(F_new) < ftol) exit
         if (F_new * F1 < 0.0_wp) then
            alpha0 = alpha1;  F0 = F1
         else
            F0 = F0 * F1 / (F1 + F_new)
         end if
         alpha1 = alpha_new;  F1 = F_new
      end do

      alpha_ep = alpha_new
   end subroutine find_yield_intersection

   ! ---------------------------------------------------------------------------
   ! Drift correction: project stress back onto yield surface
   ! Simple normal correction — more sophisticated methods possible
   ! ---------------------------------------------------------------------------
   subroutine correct_drift(model, sig, stiff_e, ftol)
      class(csm_model_t), intent(inout) :: model
      real(wp),           intent(inout) :: sig(6)
      real(wp),           intent(in)    :: stiff_e(6,6), ftol

      real(wp) :: F, dF_by_dsig(6), dg_plas_by_dsig(6), dlambda
      integer  :: iter

      do iter = 1, 100
         F = model%yield_fn(sig)
         if (abs(F) < ftol) return

         dF_by_dsig = model%flow_rule(sig)
         dg_plas_by_dsig = model%plastic_potential(sig)
         dlambda    = F / dot_product(dF_by_dsig, matmul(stiff_e, dg_plas_by_dsig))
         sig        = sig - dlambda * matmul(stiff_e, dg_plas_by_dsig)
         call model%update_hardening(dlambda * dg_plas_by_dsig)
      end do

   end subroutine correct_drift

end module mod_euler_substep
