! Central UMAT entry point — Abaqus standard interface.
!
! This is the single subroutine Abaqus links against.  Its only job is
! to route the call to the correct model implementation based on the
! material name portion of CMNAME.
!
! CMNAME convention:  MATERIALNAME[_INTEG_INTEGRATORNAME]
!
!   LINEAR_ELASTIC                → linear elastic (default integrator)
!   LINEAR_ELASTIC_INTEG_EULER    → linear elastic, Euler substepping
!   MCSS_INTEG_EULER              → Mohr-Coulomb SS, Euler substepping
!   MCSS_INTEG_ORTIZ_SIMO         → Mohr-Coulomb SS, Ortiz-Simo
!
! Abaqus uppercases CMNAME before passing it; model names are matched
! in uppercase.  The integrator suffix is lowercased inside each model UMAT.

subroutine UMAT(STRESS, STATEV, DDSDDE,          &
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

   use mod_cmname_parser,       only: material_name
   use mod_umat_linear_elastic, only: umat_linear
   use mod_umat_mcss,           only: umat_mcss
   use mod_csm_kinds,           only: wp

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

   select case(trim(material_name(CMNAME)))

   case("LINEAR_ELASTIC")
      call umat_linear(STRESS, STATEV, DDSDDE,          &
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

   case("MCSS")
      call umat_mcss(STRESS, STATEV, DDSDDE,          &
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

   case default
      error stop "UMAT: unknown material name in CMNAME"

   end select

end subroutine UMAT
