module mod_umat_mcss
   implicit none
   private
   public :: umat_mcss

contains

! UMAT wrapper for the MCSS model under the new architecture.
!
! Responsibilities of this wrapper — and ONLY this wrapper:
!   1. Unpack PROPS/STATEV into typed model fields
!   2. Determine problem dimensionality from NDI/NSHR
!   3. Reorder Voigt components to internal convention
!   4. Inflate solver vectors (NTENS) to internal 6-component form
!   5. For plane stress: Newton iteration to find eps_33 (sig_33 = 0 constraint)
!   6. Call the integrator
!   7. Deflate and reorder back to solver convention
!   8. Pack model state back into STATEV
!   9. Compute DDSDDE
!
! Everything below this layer uses named fields and 6-component internal vectors.
!
! enforce_plane_stress is an internal subroutine (after contains) so it:
!   - inherits all use-associated names including dp
!   - has an automatic explicit interface (needed for the polymorphic argument)
! =============================================================================

subroutine umat_mcss(STRESS, STATEV, DDSDDE,       &
                                     SSE, SPD, SCD,                 &
                                     RPL, DDSDDT, DRPLDE, DRPLDT,  &
                                     STRAN, DSTRAN,                 &
                                     TIME, DTIME, TEMP, DTEMP,     &
                                     PREDEF, DPRED, CMNAME,        &
                                     NDI, NSHR, NTENS, NSTATEV,    &
                                     PROPS, NPROPS,                 &
                                     COORDS, DROT, PNEWDT, CELENT, &
                                     DFGRD0, DFGRD1,               &
                                     NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

   use mod_csm_kinds,          only: wp
   use mod_mcss_model,        only: mcss_model_t, mcss_from_props, &
                                    mcss_load_state, mcss_save_state
   use mod_euler_substep,     only: euler_substep
   use mod_integrate_stress,  only: integrate_stress
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
   type(mcss_model_t) :: model
   integer  :: ptype
   real(wp) :: sig6(6), dstran6(6)
   real(wp) :: D6(6,6)

   ! -----------------------------------------------------------------------
   ! 1. Build model from PROPS + STATEV
   ! -----------------------------------------------------------------------
   model = mcss_from_props(PROPS)
   call mcss_load_state(model, STATEV)

   ! -----------------------------------------------------------------------
   ! 2. Determine problem dimensionality
   ! -----------------------------------------------------------------------
   ptype = problem_type(NDI, NSHR)

   ! -----------------------------------------------------------------------
   ! 3 & 4. Reorder to internal convention then inflate to 6 components
   !        to_internal handles the Voigt component ordering (e.g. Anura3D swap)
   !        inflate handles the 4->6 or 3->6 padding
   ! -----------------------------------------------------------------------
   sig6    = inflate(to_internal(STRESS(1:NTENS), ANURA3D_ORDER(1:NTENS)), ptype)
   dstran6 = inflate(to_internal(DSTRAN(1:NTENS), ANURA3D_ORDER(1:NTENS)), ptype)

   ! -----------------------------------------------------------------------
   ! 5. Plane stress: Newton iteration to find eps_33
   !    For plane strain and axisymmetric eps_33 is zero — inflate already
   !    set component 3 of dstran6 correctly (zero or solver-provided).
   !    For plane stress eps_33 is unknown; we solve for it so that sig_33 = 0
   !    after the stress update. The model is called repeatedly — the integrator
   !    and model are unaware this iteration is happening.
   ! -----------------------------------------------------------------------
   if (ptype == PROBLEM_PLANE_STRESS) then
      call enforce_plane_stress(model, sig6, dstran6, PROPS(11))
   end if

   ! -----------------------------------------------------------------------
   ! 6. Integrate — integrator and model work entirely in 6-component space
   ! -----------------------------------------------------------------------
   call integrate_stress(model, sig6, dstran6,              &
                         ftol=model%yield_tol,              &
                         stol=1.0e-4_wp,                    &
                         method=integrator_name(CMNAME))

   ! -----------------------------------------------------------------------
   ! 7. Deflate and reorder back to solver convention
   ! -----------------------------------------------------------------------
   STRESS(1:NTENS) = from_internal(deflate(sig6, ptype), ANURA3D_ORDER(1:NTENS))

   ! -----------------------------------------------------------------------
   ! 8. Pack updated state
   ! -----------------------------------------------------------------------
   call mcss_save_state(model, STATEV)

   ! -----------------------------------------------------------------------
   ! 9. Tangent stiffness: compute in 6-component space, deflate to NTENS x NTENS
   ! -----------------------------------------------------------------------
   D6 = model%elastic_stiffness()
   DDSDDE(1:NTENS, 1:NTENS) = deflate_stiffness(D6, ptype)

contains

   ! ==========================================================================
   ! Plane stress enforcement — wrapper-level Newton iteration.
   ! Finds eps_33 such that sig_33 = 0 after a full stress update.
   ! Internal subroutine: inherits dp and use-associations from the host.
   !
   ! Algorithm:
   !   Given sig_33 = f(eps_33), solve f(eps_33) = 0 by Newton's method.
   !   Jacobian df/d(eps_33) approximated from the elastic stiffness column.
   ! ==========================================================================
   subroutine enforce_plane_stress(mdl, sig6_ps, dstran6_ps, ftol)
      class(mcss_model_t), intent(inout) :: mdl
      real(wp),            intent(inout) :: sig6_ps(6), dstran6_ps(6)
      real(wp),            intent(in)    :: ftol

      integer,  parameter :: MAX_ITER = 20
      real(wp), parameter :: TOL = 1.0e-10_wp

      real(wp) :: sig_trial(6), dstran_trial(6)
      real(wp) :: sig33, d_sig33_d_eps33
      real(wp) :: D6_ps(6,6)
      real(wp) :: deps33
      real(wp), allocatable :: state_saved(:)
      integer  :: iter

      ! Elastic stiffness column 3 gives dsig_33/deps_33 (approximate Jacobian)
      D6_ps           = mdl%elastic_stiffness()
      d_sig33_d_eps33 = D6_ps(3, 3)

      ! Initial guess: eps_33 = 0 (already set by inflate)
      do iter = 1, MAX_ITER

         ! Trial update with current eps_33 guess
         sig_trial    = sig6_ps
         dstran_trial = dstran6_ps

         call mdl%snapshot(state_saved)
         call euler_substep(mdl, sig_trial, dstran_trial, ftol=ftol, stol=1.0e-4_wp)
         call mdl%restore(state_saved)

         ! Check sig_33
         sig33 = sig_trial(3)
         if (abs(sig33) < TOL) then
            dstran6_ps = dstran_trial
            sig6_ps    = sig_trial
            return
         end if

         ! Newton update: eps_33 -= sig_33 / (dsig_33/deps_33)
         deps33           = -sig33 / d_sig33_d_eps33
         dstran6_ps(3)    = dstran6_ps(3) + deps33

      end do

      ! Converged solution already stored in sig_trial/dstran_trial on last accepted iter
      sig6_ps    = sig_trial
      dstran6_ps = dstran_trial

   end subroutine enforce_plane_stress

end subroutine umat_mcss

end module mod_umat_mcss
