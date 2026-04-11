!! Generic cutting-plane algorithm (CPA) stress integrator.
!!
!! Implements Algorithm 3 of Ortiz & Simo (1986): "An analysis of a new class
!! of integration algorithms for elastoplastic constitutive relations".
!!
!! Knows nothing about any specific model — calls only the deferred procedures
!! defined in `csm_model_t`. The elastic predictor is applied first; if it
!! is inside the yield surface the step is purely elastic. Otherwise Newton
!! iterations drive the stress back to the yield surface:
!!
!!   dlambda     = F / (n·De·m - H)
!!   sig         = sig - dlambda · De·m
!!   update_hardening(dlambda · m)
!!
!! where n = dF/dsig, m = dG/dsig, H = hardening modulus.
!!
!! ### Differences from euler_substep
!!
!! - No substepping: the full strain increment is applied as one elastic
!!   predictor; convergence comes from Newton iteration, not step splitting.
!! - `stol` and `dt_min` from `integrator_params_t` are not used.
!! - `max_iters` caps the Newton iteration count (default 100).
!! - For constant-modulus models `De` is computed once and held fixed (secant
!!   approximation). For pressure-dependent models (e.g. NorSand) `pre_step`
!!   is called at the start of each Newton iteration and `De` is recomputed,
!!   giving a full tangent iteration. For models with a no-op `pre_step` the
!!   behaviour is identical to the original secant algorithm.

module mod_cpa
   use mod_csm_kinds,     only: wp
   use mod_csm_model,     only: csm_model_t
   use mod_euler_substep, only: integrator_params_t, DEFAULT_INTEGRATOR_PARAMS
   implicit none
   private

   public :: cpa_step

contains

   subroutine cpa_step(model, sig, deps, iparams)
      !! Integrate one strain increment using the cutting-plane algorithm.
      class(csm_model_t),        intent(inout)        :: model
      real(wp),                  intent(inout)        :: sig(6)   !! stress, updated in place
      real(wp),                  intent(in)           :: deps(6)  !! strain increment
      type(integrator_params_t), intent(in), optional :: iparams

      type(integrator_params_t) :: params
      real(wp) :: stiff_e(6,6), De_m(6)
      real(wp) :: dF_by_dsig(6), dg_by_dsig(6)
      real(wp) :: F, H, denominator, dlambda
      integer  :: iter

      if (present(iparams)) then
         params = iparams
      else
         params = DEFAULT_INTEGRATOR_PARAMS
      end if

      ! -----------------------------------------------------------------------
      ! 1. Full elastic predictor
      ! -----------------------------------------------------------------------
      call model%pre_step(sig)
      stiff_e = model%elastic_stiffness()
      sig     = sig + matmul(stiff_e, deps)

      ! -----------------------------------------------------------------------
      ! 2. Check yield — purely elastic step
      ! -----------------------------------------------------------------------
      F = model%yield_fn(sig)
      if (F < params%ftol) return

      ! -----------------------------------------------------------------------
      ! 3. Cutting-plane Newton iterations
      !    Each iteration computes dlambda from the current stress state and
      !    moves sig toward the yield surface. update_hardening evolves the
      !    model's internal state (eps_p and any derived quantities).
      ! -----------------------------------------------------------------------
      do iter = 1, params%max_iters

         ! Update state and stiffness for the current stress — no-op for
         ! constant-modulus models; recomputes G(p) for pressure-dependent ones.
         call model%pre_step(sig)
         stiff_e = model%elastic_stiffness()

         dg_by_dsig  = model%plastic_potential(sig)
         dF_by_dsig  = model%flow_rule(sig)
         H           = model%hardening_modulus(sig, dg_by_dsig)

         De_m        = matmul(stiff_e, dg_by_dsig)
         denominator = dot_product(dF_by_dsig, De_m) - H
         dlambda     = F / denominator

         sig = sig - dlambda * De_m
         call model%update_hardening(dlambda * dg_by_dsig)

         F = model%yield_fn(sig)
         if (abs(F) <= params%ftol) return

      end do

      ! Exhausted max_iters without converging — return current stress state.
      ! Caller can detect non-convergence by checking model%yield_fn(sig) > ftol.

   end subroutine cpa_step

end module mod_cpa
