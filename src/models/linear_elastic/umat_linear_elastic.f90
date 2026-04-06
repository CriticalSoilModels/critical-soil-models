! UMAT wrapper for the linear elastic model under the new architecture.
!
! PROPS layout:
!   1: G   shear modulus [kPa]
!   2: nu  Poisson's ratio [-]
!
! No STATEV — linear elastic has no internal state.

module mod_umat_linear_elastic
   implicit none
   private
   public :: umat_linear

contains

subroutine umat_linear(STRESS, STATEV, DDSDDE,          &
                        SSE, SPD, SCD,                   &
                        RPL, DDSDDT, DRPLDE, DRPLDT,    &
                        STRAN, DSTRAN,                   &
                        TIME, DTIME, TEMP, DTEMP,        &
                        PREDEF, DPRED, CMNAME,           &
                        NDI, NSHR, NTENS, NSTATEV,       &
                        PROPS, NPROPS,                   &
                        COORDS, DROT, PNEWDT, CELENT,    &
                        DFGRD0, DFGRD1,                  &
                        NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

   use mod_csm_kinds,             only: wp
   use mod_linear_elastic_model,  only: linear_elastic_model_t, le_from_props
   use mod_integrate_stress,      only: integrate_stress, integrator_params_t, DEFAULT_INTEGRATOR_PARAMS
   use mod_cmname_parser,         only: integrator_name
   use mod_voigt_conventions,     only: to_internal, from_internal, ANURA3D_ORDER, &
                                        problem_type, inflate, deflate, deflate_stiffness

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
   type(linear_elastic_model_t) :: model
   integer  :: ptype
   real(wp) :: sig6(6), dstran6(6), D6(6,6)

   ! 1. Build model from PROPS (no state to load)
   model = le_from_props(PROPS)

   ! 2. Problem dimensionality
   ptype = problem_type(NDI, NSHR)

   ! 3 & 4. Reorder + inflate to 6-component internal form
   sig6    = inflate(to_internal(STRESS(1:NTENS),  ANURA3D_ORDER(1:NTENS)), ptype)
   dstran6 = inflate(to_internal(DSTRAN(1:NTENS),  ANURA3D_ORDER(1:NTENS)), ptype)

   ! 5. Integrate (F = -1 always → exits after elastic predictor)
   call integrate_stress(model, sig6, dstran6, method=integrator_name(CMNAME), &
                         iparams=DEFAULT_INTEGRATOR_PARAMS)

   ! 6. Deflate + reorder back to solver convention
   STRESS(1:NTENS) = from_internal(deflate(sig6, ptype), ANURA3D_ORDER(1:NTENS))

   ! 7. Tangent stiffness
   D6 = model%elastic_stiffness()
   DDSDDE(1:NTENS, 1:NTENS) = deflate_stiffness(D6, ptype)

end subroutine umat_linear

end module mod_umat_linear_elastic
