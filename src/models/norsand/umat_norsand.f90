!! UMAT wrapper for the NorSand critical state model.
!!
!! Responsibilities of this wrapper — and **only** this wrapper:
!!
!! 1. Unpack PROPS/STATEV into typed model fields
!! 2. Determine problem dimensionality from NDI/NSHR
!! 3. Reorder Voigt components to internal convention
!! 4. Inflate solver vectors (NTENS) to internal 6-component form
!! 5. Call the stress integrator (pre_step called internally each substep)
!! 6. Deflate and reorder back to solver convention
!! 7. Pack model state back into STATEV
!! 8. Return elastic tangent stiffness DDSDDE
!!
!! ### PROPS layout (13 required + 2 optional entries)
!!
!! See `norsand_model.f90` module doc for the full table.
!! PROPS(14) = ftol (yield tolerance; default 1e-8 when absent or ≤ 0).
!! PROPS(15) = max_iters (substep limit; default 500 when absent or ≤ 0).
!!
!! ### STATEV layout (15 entries)
!!
!! See `norsand_model.f90` module doc for the full table.
!!
!! ### Notes
!!
!! NorSand is designed for fully 3D stress states (geotechnical).
!! Plane stress is not expected in practice but the wrapper handles it
!! via the same Newton iteration used in other model UMATs.
!!
!! `pre_step` (pressure-dependent G, K and Lode-dependent M_i) is called
!! inside `euler_substep` / `cpa_step` — no explicit call is needed here.

module mod_umat_norsand
   implicit none
   private
   public :: umat_norsand

contains

subroutine umat_norsand(STRESS, STATEV, DDSDDE,       &
                        SSE, SPD, SCD,                 &
                        RPL, DDSDDT, DRPLDE, DRPLDT,  &
                        STRAN, DSTRAN,                 &
                        TIME, DTIME, TEMP, DTEMP,      &
                        PREDEF, DPRED, CMNAME,         &
                        NDI, NSHR, NTENS, NSTATEV,     &
                        PROPS, NPROPS,                 &
                        COORDS, DROT, PNEWDT, CELENT,  &
                        DFGRD0, DFGRD1,                &
                        NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

   use mod_csm_kinds,         only: wp
   use mod_norsand_model,     only: norsand_model_t, norsand_from_props, &
                                    norsand_load_state, norsand_save_state
   use mod_euler_substep,     only: euler_substep
   use mod_integrate_stress,  only: integrate_stress, integrator_params_t
   use mod_cmname_parser,     only: integrator_name
   use mod_voigt_conventions, only: to_internal, from_internal, ANURA3D_ORDER, &
                                    problem_type, inflate, deflate, deflate_stiffness, &
                                    PROBLEM_PLANE_STRESS

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
   type(norsand_model_t)     :: model
   type(integrator_params_t) :: iparams
   integer  :: ptype
   real(wp) :: sig6(6), dstran6(6)
   real(wp) :: D6(6,6)

   real(wp), parameter :: FTOL_DEFAULT      = 1.0e-8_wp
   real(wp), parameter :: STOL_DEFAULT      = 1.0e-4_wp
   real(wp), parameter :: DT_MIN_DEFAULT    = 1.0e-9_wp
   integer,  parameter :: MAX_ITERS_DEFAULT = 500

   ! -----------------------------------------------------------------------
   ! 1. Build model from PROPS + STATEV
   ! -----------------------------------------------------------------------
   model = norsand_from_props(PROPS)
   call norsand_load_state(model, STATEV)

   ! PROPS(14) = ftol      (optional; use default when absent or <= 0)
   ! PROPS(15) = max_iters (optional; use default when absent or <= 0)
   iparams%stol   = STOL_DEFAULT
   iparams%dt_min = DT_MIN_DEFAULT

   if (NPROPS >= 14 .and. PROPS(14) > 0.0_wp) then
      iparams%ftol = PROPS(14)
   else
      iparams%ftol = FTOL_DEFAULT
   end if

   if (NPROPS >= 15 .and. PROPS(15) > 0.0_wp) then
      iparams%max_iters = int(PROPS(15))
   else
      iparams%max_iters = MAX_ITERS_DEFAULT
   end if

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
   ! 5. Integrate — pre_step is called by the integrator each substep
   ! -----------------------------------------------------------------------
   if (ptype == PROBLEM_PLANE_STRESS) then
      call enforce_plane_stress(model, sig6, dstran6, iparams)
   else
      call integrate_stress(model, sig6, dstran6,             &
                            method=integrator_name(CMNAME),   &
                            iparams=iparams)
   end if

   ! -----------------------------------------------------------------------
   ! 6. Deflate and reorder back to solver convention
   ! -----------------------------------------------------------------------
   STRESS(1:NTENS) = from_internal(deflate(sig6, ptype), ANURA3D_ORDER(1:NTENS))

   ! -----------------------------------------------------------------------
   ! 7. Pack updated state
   ! -----------------------------------------------------------------------
   call norsand_save_state(model, STATEV)

   ! -----------------------------------------------------------------------
   ! 8. Tangent stiffness
   ! -----------------------------------------------------------------------
   D6 = model%consistent_tangent()
   DDSDDE(1:NTENS, 1:NTENS) = deflate_stiffness(D6, ptype)

contains

   ! ==========================================================================
   ! Plane stress enforcement — Newton iteration to find eps_33 such that
   ! sig_33 = 0 after a full stress update.
   ! ==========================================================================
   subroutine enforce_plane_stress(mdl, sig6_ps, dstran6_ps, ps_iparams)
      class(norsand_model_t),    intent(inout) :: mdl
      real(wp),                  intent(inout) :: sig6_ps(6), dstran6_ps(6)
      type(integrator_params_t), intent(in)    :: ps_iparams

      integer,  parameter :: MAX_ITER = 20
      real(wp), parameter :: TOL = 1.0e-10_wp

      real(wp) :: sig_trial(6), dstran_trial(6)
      real(wp) :: sig33, d_sig33_d_eps33, deps33
      real(wp) :: D6_ps(6,6)
      real(wp), allocatable :: state_saved(:)
      integer  :: iter

      D6_ps           = mdl%elastic_stiffness()
      d_sig33_d_eps33 = D6_ps(3, 3)

      do iter = 1, MAX_ITER
         sig_trial    = sig6_ps
         dstran_trial = dstran6_ps

         call mdl%snapshot(state_saved)
         call euler_substep(mdl, sig_trial, dstran_trial, iparams=ps_iparams)
         call mdl%restore(state_saved)

         sig33 = sig_trial(3)
         if (abs(sig33) < TOL) then
            sig6_ps    = sig_trial
            dstran6_ps = dstran_trial
            return
         end if

         deps33        = -sig33 / d_sig33_d_eps33
         dstran6_ps(3) = dstran6_ps(3) + deps33
      end do

      sig6_ps    = sig_trial
      dstran6_ps = dstran_trial

   end subroutine enforce_plane_stress

end subroutine umat_norsand

end module mod_umat_norsand
