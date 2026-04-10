!! UMAT wrapper for the Mohr-Coulomb Strain Rate (MCSR) model.
!!
!! Responsibilities of this wrapper — and **only** this wrapper:
!!
!! 1. Unpack PROPS/STATEV into typed model fields
!! 2. Determine problem dimensionality from NDI/NSHR
!! 3. Reorder Voigt components to internal convention
!! 4. Inflate solver vectors (NTENS) to internal 6-component form
!! 5. Compute scalar strain rate and call `mcsr_update_rate_state` (pre-step)
!! 6. For plane stress: Newton iteration to enforce σ₃₃ = 0
!! 7. Call the stress integrator
!! 8. Deflate and reorder back to solver convention
!! 9. Pack model state back into STATEV
!! 10. Return elastic tangent stiffness DDSDDE

module mod_umat_mcsr
   implicit none
   private
   public :: umat_mcsr

contains

subroutine umat_mcsr(STRESS, STATEV, DDSDDE,          &
                     SSE, SPD, SCD,                    &
                     RPL, DDSDDT, DRPLDE, DRPLDT,     &
                     STRAN, DSTRAN,                    &
                     TIME, DTIME, TEMP, DTEMP,         &
                     PREDEF, DPRED, CMNAME,            &
                     NDI, NSHR, NTENS, NSTATEV,        &
                     PROPS, NPROPS,                    &
                     COORDS, DROT, PNEWDT, CELENT,     &
                     DFGRD0, DFGRD1,                   &
                     NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

   use mod_csm_kinds,          only: wp
   use mod_mcsr_model,         only: mcsr_model_t, mcsr_from_props, &
                                     mcsr_load_state, mcsr_save_state
   use mod_mcsr_functions,     only: mcsr_update_rate_state
   use mod_euler_substep,      only: euler_substep
   use mod_integrate_stress,   only: integrate_stress, integrator_params_t
   use mod_cmname_parser,      only: integrator_name
   use mod_voigt_conventions,  only: to_internal, from_internal, ANURA3D_ORDER, &
                                     problem_type, inflate, deflate, deflate_stiffness, &
                                     PROBLEM_PLANE_STRESS
   use mod_strain_invariants,  only: calc_eps_vol_inv, calc_dev_strain, calc_eps_q_inv

   implicit none

   ! --- Standard UMAT arguments ---
   character(80) :: CMNAME
   integer :: NDI, NSHR, NTENS, NSTATEV, NPROPS
   integer :: NOEL, NPT, LAYER, KSPT, KSTEP, KINC
   real(wp) :: SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP, DTEMP, PNEWDT, CELENT
   real(wp) :: STRESS(NTENS), STATEV(NSTATEV), DDSDDE(NTENS,NTENS)
   real(wp) :: DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS)
   real(wp) :: TIME(2), PREDEF(1), DPRED(1), COORDS(3)
   real(wp) :: DROT(3,3), DFGRD0(3,3), DFGRD1(3,3)
   real(wp) :: PROPS(NPROPS)

   ! --- Local ---
   type(mcsr_model_t)        :: model
   type(integrator_params_t) :: iparams
   integer  :: ptype
   real(wp) :: sig6(6), dstran6(6), dstran6_rate(6)
   real(wp) :: D6(6,6)
   real(wp) :: e_rate_scalar, eps_vol, dev(6)

   ! -----------------------------------------------------------------------
   ! 1. Build model from PROPS + STATEV
   ! -----------------------------------------------------------------------
   model = mcsr_from_props(PROPS)
   call mcsr_load_state(model, STATEV)

   iparams%ftol   = model%params%yield_tol
   iparams%stol   = 1.0e-4_wp
   iparams%dt_min = 1.0e-9_wp

   ! -----------------------------------------------------------------------
   ! 2. Determine problem dimensionality
   ! -----------------------------------------------------------------------
   ptype = problem_type(NDI, NSHR)

   ! -----------------------------------------------------------------------
   ! 3 & 4. Reorder to internal convention then inflate to 6 components
   ! -----------------------------------------------------------------------
   sig6    = inflate(to_internal(STRESS(1:NTENS), ANURA3D_ORDER(1:NTENS)), ptype)
   dstran6 = inflate(to_internal(DSTRAN(1:NTENS), ANURA3D_ORDER(1:NTENS)), ptype)

   ! -----------------------------------------------------------------------
   ! 5. Pre-step rate state update
   !    Compute scalar strain rate magnitude from inflated dstran6 / DTIME.
   !    This must happen before the integrator so G, K, eta_y, M_lode are
   !    current for the upcoming substeps.
   ! -----------------------------------------------------------------------
   if (DTIME > 0.0_wp) then
      dstran6_rate = dstran6 / DTIME
   else
      dstran6_rate = 0.0_wp
   end if

   ! Deviatoric strain rate magnitude (eps_q of the rate tensor)
   eps_vol      = calc_eps_vol_inv(dstran6_rate)
   dev          = calc_dev_strain(dstran6_rate, eps_vol)
   e_rate_scalar = calc_eps_q_inv(dev)

   call mcsr_update_rate_state(model%params, model%state, e_rate_scalar, sig6)

   ! -----------------------------------------------------------------------
   ! 6. Plane stress: Newton iteration to find eps_33
   ! -----------------------------------------------------------------------
   if (ptype == PROBLEM_PLANE_STRESS) then
      call enforce_plane_stress(model, sig6, dstran6, model%params%yield_tol)
   end if

   ! -----------------------------------------------------------------------
   ! 7. Integrate
   ! -----------------------------------------------------------------------
   call integrate_stress(model, sig6, dstran6,             &
                         method=integrator_name(CMNAME),   &
                         iparams=iparams)

   ! -----------------------------------------------------------------------
   ! 8. Deflate and reorder back to solver convention
   ! -----------------------------------------------------------------------
   STRESS(1:NTENS) = from_internal(deflate(sig6, ptype), ANURA3D_ORDER(1:NTENS))

   ! -----------------------------------------------------------------------
   ! 9. Pack updated state
   ! -----------------------------------------------------------------------
   call mcsr_save_state(model, STATEV)

   ! -----------------------------------------------------------------------
   ! 10. Tangent stiffness: elastic tangent in 6-component space, deflated
   ! -----------------------------------------------------------------------
   D6 = model%elastic_stiffness()
   DDSDDE(1:NTENS, 1:NTENS) = deflate_stiffness(D6, ptype)

contains

   ! ==========================================================================
   ! Plane stress enforcement — wrapper-level Newton iteration.
   ! Finds eps_33 such that sig_33 = 0 after a full stress update.
   ! ==========================================================================
   subroutine enforce_plane_stress(mdl, sig6_ps, dstran6_ps, ftol)
      class(mcsr_model_t), intent(inout) :: mdl
      real(wp),            intent(inout) :: sig6_ps(6), dstran6_ps(6)
      real(wp),            intent(in)    :: ftol

      integer,  parameter :: MAX_ITER = 20
      real(wp), parameter :: TOL = 1.0e-10_wp

      real(wp) :: sig_trial(6), dstran_trial(6)
      real(wp) :: sig33, d_sig33_d_eps33
      real(wp) :: D6_ps(6,6), deps33
      real(wp), allocatable :: state_saved(:)
      integer  :: iter

      D6_ps           = mdl%elastic_stiffness()
      d_sig33_d_eps33 = D6_ps(3, 3)

      do iter = 1, MAX_ITER
         sig_trial    = sig6_ps
         dstran_trial = dstran6_ps

         call mdl%snapshot(state_saved)
         call euler_substep(mdl, sig_trial, dstran_trial, &
                            iparams=integrator_params_t(ftol=ftol, stol=1.0e-4_wp, dt_min=1.0e-9_wp))
         call mdl%restore(state_saved)

         sig33 = sig_trial(3)
         if (abs(sig33) < TOL) then
            dstran6_ps = dstran_trial
            sig6_ps    = sig_trial
            return
         end if

         deps33        = -sig33 / d_sig33_d_eps33
         dstran6_ps(3) = dstran6_ps(3) + deps33
      end do

      sig6_ps    = sig_trial
      dstran6_ps = dstran_trial
   end subroutine enforce_plane_stress

end subroutine umat_mcsr

end module mod_umat_mcsr
