!! Generic cutting-plane return mapping (CPRM) stress integrator.
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
!! - `De` is computed once from the state at the start of the increment and
!!   held fixed throughout (secant approximation — exact for constant-modulus
!!   models such as MCSR where moduli are pre-updated before integration).

module mod_cprm
   use mod_csm_kinds,     only: wp
   use mod_csm_model,     only: csm_model_t
   use mod_euler_substep, only: integrator_params_t, DEFAULT_INTEGRATOR_PARAMS
   implicit none
   private

   public :: cprm_step

contains

   subroutine cprm_step(model, sig, deps, iparams)
      !! Integrate one strain increment using the cutting-plane return mapping.
      class(csm_model_t),        intent(inout)        :: model
      real(wp),                  intent(inout)        :: sig(6)   !! stress, updated in place
      real(wp),                  intent(in)           :: deps(6)  !! strain increment
      type(integrator_params_t), intent(in), optional :: iparams

      type(integrator_params_t) :: p
      real(wp) :: stiff_e(6,6), De_m(6)
      real(wp) :: dF_by_dsig(6), dg_by_dsig(6)
      real(wp) :: F, H, denominator, dlambda
      integer  :: iter

      if (present(iparams)) then
         p = iparams
      else
         p = DEFAULT_INTEGRATOR_PARAMS
      end if

      ! -----------------------------------------------------------------------
      ! 1. Full elastic predictor
      !    De is computed once and held fixed for the entire iteration.
      ! -----------------------------------------------------------------------
      stiff_e = model%elastic_stiffness()
      sig     = sig + matmul(stiff_e, deps)

      ! -----------------------------------------------------------------------
      ! 2. Check yield — purely elastic step
      ! -----------------------------------------------------------------------
      F = model%yield_fn(sig)
      if (F < p%ftol) return

      ! -----------------------------------------------------------------------
      ! 3. Cutting-plane Newton iterations
      !    Each iteration computes dlambda from the current stress state and
      !    moves sig toward the yield surface. update_hardening evolves the
      !    model's internal state (eps_p and any derived quantities).
      ! -----------------------------------------------------------------------
      do iter = 1, p%max_iters

         dg_by_dsig  = model%plastic_potential(sig)
         dF_by_dsig  = model%flow_rule(sig)
         H           = model%hardening_modulus(sig, dg_by_dsig)

         De_m        = matmul(stiff_e, dg_by_dsig)
         denominator = dot_product(dF_by_dsig, De_m) - H
         dlambda     = F / denominator

         sig = sig - dlambda * De_m
         call model%update_hardening(dlambda * dg_by_dsig)

         F = model%yield_fn(sig)
         if (abs(F) <= p%ftol) return

      end do

      ! Exhausted max_iters without converging — return current stress state.
      ! Caller can detect non-convergence by checking model%yield_fn(sig) > ftol.

   end subroutine cprm_step

end module mod_cprm
