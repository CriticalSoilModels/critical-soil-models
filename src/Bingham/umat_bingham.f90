module mod_bingham
   use stdlib_kinds, only: dp
   ! use mod_strain_invariants, only: calc_eps_vol_invariant

contains

   subroutine esm_bingham(npt,noel,idset,stress,eunloading,plasticmultiplier, &
      dstran,nstatev,statev,naddvar,additionalvar,cmname,nprops,props,numberofphases,ntens)


      ! implicit real(dp) (a-h, o-z)
      integer :: ntens, nstatev, naddvar, nprops, npt, noel, idset, numberofphases
      real(dp) :: eunloading, plasticmultiplier
      character(len=80), intent(in):: cmname
      real(dp), intent(inout) :: stress(ntens), dstran(ntens),statev(nstatev),additionalvar(naddvar),props(nprops)


!---local variables required in standard umat
      integer :: istep, timestep
      real(dp), dimension(:), allocatable :: ddsddt ! only for fully coupled thermal analysis: variation of stress increment due to temperature
      real(dp), dimension(:), allocatable :: drplde ! only for fully coupled thermal analysis: variation of volumetric heat generation due to strain increment
      real(dp), dimension(:), allocatable :: stran
      real(dp), dimension(:), allocatable :: time
      real(dp), dimension(:), allocatable :: predef
      real(dp), dimension(:), allocatable :: dpred
      real(dp), dimension(:), allocatable :: coords
      real(dp), dimension(:,:), allocatable :: ddsdde ! jacobian matrix of the constitutive model (tangent stiffness matrix in case of mc)
      real(dp), dimension(:,:), allocatable :: drot
      real(dp), dimension(:,:), allocatable :: dfgrd0
      real(dp), dimension(:,:), allocatable :: dfgrd1
      real(dp) :: sse, spd, scd ! specific elastic strain energy, plastic dissipation, creep dissipation
      real(dp) :: rpl ! only for fully coupled thermal analysis: volumetric heat generation
      real(dp) :: drpldt ! only for fully coupled thermal analysis: variation of volumetric heat generation due to temperature
      real(dp) :: pnewdt, dtime, temp, dtemp, celent
      real(dp) :: value ! auxiliary variable holding any real valued number
      real(dp) :: porosity, waterpressure, waterpressure0, gaspressure, gaspressure0, degreesaturation


      integer :: ndi, nshr, layer, kspt, kstep, kinc

!---local variables defned by the user
! e.g. integer :: var_local
!---user can define here additional variables

      allocate( ddsddt(ntens), drplde(ntens), stran(ntens), time(2), predef(1), dpred(1),  &
         coords(3), ddsdde(ntens,ntens), drot(3,3), dfgrd0(3,3), dfgrd1(3,3) )

!initialization
      eunloading = 0.0
      plasticmultiplier = 0.0

!rename additional variables
      porosity = additionalvar(1)
      waterpressure = additionalvar(2)
      waterpressure0 = additionalvar(3)
      gaspressure = additionalvar(4)
      gaspressure0 = additionalvar(5)
      degreesaturation = additionalvar(6)
      time(1) = additionalvar(7)   !totalrealtime
      time(2) = additionalvar(8)   !overalltotaltime
      dtime = additionalvar(9)     !timeincrement
      istep = additionalvar(10)
      timestep = additionalvar(11)   !note: very first time and load step: istep=1 and timestep=1
!call the umat
      call umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, &
         dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatev, props, nprops, coords, drot, pnewdt, celent, dfgrd0, &
         dfgrd1, noel, npt, layer, kspt, kstep, kinc)


!---definition of eunloading -> required to define the max time step
      eunloading = props(5)
!---always define this value to run the simulation

! plasticmultiplier can be given as an output because plastic points can be plotted as a result
   end subroutine esm_bingham

   subroutine umat(stress,statev,ddsdde,sse,spd,scd, &
      rpl,ddsddt,drplde,drpldt, &
      stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
      ndi,nshr,ntens,nstatev,props,nprops,coords,drot,pnewdt, &
      celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

      implicit none

      character(len=80) cmname
      integer :: ntens,nstatev,nprops,ndi,nshr,noel, &
         npt,layer,kspt,kstep,kinc
      real(dp) :: sse,spd,scd,rpl,drpldt,dtime,temp,dtemp, &
         pnewdt,celent
      real(dp) :: stress(ntens),statev(nstatev), &
         ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens), &
         stran(ntens),dstran(ntens),time(2),predef(1),dpred(1), &
         props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)


      ! arguments:
      !          i/o  type
      !  props    i   r()  : list with model parameters
      !  dstran   i   r()  : strain increment
      !  ddsdde   o   r(,) : material stiffness matrix
      !  stress  i/o  r()  : stresses
      !  statev  i/o  r()  : state variables

      !---  local variables
      real(dp) :: dsig(ntens), sig(ntens), strainrate(ntens), devstrainrate(ntens), sigma(ntens)
      real(dp) :: dstranvol
      real(dp):: eel, enu, yieldstress, viscosity, bulkmodulusliquid, liquidpressurecavitation, one, two, g
      real(dp):: shearstress, fac, d1,d2
      real(dp):: volstrainrate, volstrainratecomponent, liquidpressure, liquidpressureincrement
      logical ::  computepressure
      integer :: i


      eel = props(1)
      enu = props(2)
      yieldstress=props(3)
      viscosity = props(4)
      bulkmodulusliquid=props(5)
      liquidpressurecavitation=props(6)
      one = 1.0d0
      two = 2.0d0
      g = eel/(two*(one + enu))

      if (ntens == 6) then ! 3d
         shearstress =sqrt (( (stress(1) + &
            -stress(2))**2 + &
            (stress(1) + &
            -stress(3))**2 + &
            (stress(2)+ &
            -stress(3))**2 ) / 3.0 + &
            ( stress(4)**2 + &
            stress(5)**2 + &
            stress(6)**2 ) * 2.0)
      else if (ntens == 4) then ! 2d
         shearstress =sqrt (( (stress(1) + &
            -stress(2))**2 + &
            (stress(1)+ &
            -stress(3))**2 + &
            (stress(2)+ &
            -stress(3))**2 ) / 3.0 + &
            ( stress(4)**2) * 2.0)
      end if

      if   (shearstress.le.sqrt(2.0)*yieldstress) then
         ! calculate elastic stress increment (dsige = elastic stiffness d * strain increment deps)
         fac = two * g / ( one - two * enu )
         d1 = fac * ( one - enu )
         d2 = fac * enu
         dstranvol = dstran(1) + dstran(2) + dstran(3)
         dsig(1) = (d1 - d2) * dstran(1) + d2 * dstranvol
         dsig(2) = (d1 - d2) * dstran(2) + d2 * dstranvol
         dsig(3) = (d1 - d2) * dstran(3) + d2 * dstranvol
         dsig(4) = g * dstran(4)

         if (ntens == 6) then
            dsig(5) = g * dstran(5)
            dsig(6) = g * dstran(6)
         end if

         ! elastic stress
         sig = stress + dsig

         ! stress state parameters update
         do i = 1, ntens
            stress(i) = sig(i)
         end do

         ddsdde = 0.0
         ddsdde(1:3,1:3) = d2
         ddsdde(1,1) = d1
         ddsdde(2,2) = d1
         ddsdde(3,3) = d1
         ddsdde(4,4) = g
         if (ntens == 6) then
            ddsdde(5,5) = g
            ddsdde(6,6) = g
         end if

         return
      else

         sigma = stress

         dstranvol = dstran(1) + dstran(2) + dstran(3)
         strainrate = dstran / dtime

         volstrainrate = strainrate(1) + strainrate(2) + strainrate(3)
         volstrainratecomponent = volstrainrate / 3.0
         devstrainrate(1:3) = strainrate(1:3) - volstrainratecomponent
         devstrainrate(4:ntens) = strainrate(4:ntens)


         ! consider stress state of previous time step
         liquidpressure = (stress(1)+stress(2)+stress(3))/3

         computepressure = (liquidpressure.lt.liquidpressurecavitation)

         if(computepressure) then

            liquidpressureincrement = bulkmodulusliquid  * dstranvol

            liquidpressure = liquidpressure + liquidpressureincrement


            if(liquidpressure .gt. liquidpressurecavitation) then
               liquidpressure=liquidpressurecavitation
               stress(1:3) = liquidpressure
               stress(4:ntens) = 0.0
               return
            end if

            sigma(1:3) = liquidpressure
            sigma(4:ntens) = 0.0

            sigma(1:3) = liquidpressure +  &
               two * viscosity * devstrainrate(1:3)
            sigma(4:ntens) = viscosity * devstrainrate(4:ntens)

         else
            liquidpressure=liquidpressurecavitation
            stress(1:3) = liquidpressure
            stress(4:ntens) = 0.0
            return

         end if

         stress = sigma
      end if
   end subroutine umat

end module mod_bingham
