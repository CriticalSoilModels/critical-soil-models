! =============================================================================
! PSEUDOCODE — does not compile
! Voigt convention translation utilities.
! Internal convention: [11, 22, 33, 12, 13, 23]  (standard mechanics order)
! Translation and inflation/deflation happen only at the UMAT boundary —
! never inside models or integrators.
!
! Problem types and their solver-side component counts:
!
!   3D           NDI=3, NSHR=3, NTENS=6   [11, 22, 33, 12, 13, 23]
!   Plane strain NDI=3, NSHR=1, NTENS=4   [11, 22, 33, 12]
!   Axisymmetric NDI=3, NSHR=1, NTENS=4   [11, 22, 33, 12]  (33 = hoop)
!   Plane stress NDI=2, NSHR=1, NTENS=3   [11, 22, 12]      (33 = 0)
!
! Internally every model and integrator always works with 6 components.
! The wrapper inflates solver vectors to 6 on entry and deflates back on exit.
! Plane stress requires a Newton iteration in the wrapper to find eps_33
! such that sig_33 = 0 — handled in the UMAT wrapper, not here.
! =============================================================================

module mod_voigt_conventions
   use iso_fortran_env, only: dp => real64
   implicit none

   ! ---------------------------------------------------------------------------
   ! Problem type identifiers
   ! ---------------------------------------------------------------------------
   integer, parameter :: PROBLEM_3D           = 1
   integer, parameter :: PROBLEM_PLANE_STRAIN = 2
   integer, parameter :: PROBLEM_AXISYMMETRIC = 3
   integer, parameter :: PROBLEM_PLANE_STRESS = 4

   ! ---------------------------------------------------------------------------
   ! Named reorder arrays for each supported solver convention
   ! Internal (this library): [11, 22, 33, 12, 13, 23]
   ! ---------------------------------------------------------------------------

   integer, parameter :: INTERNAL_ORDER(6) = [1, 2, 3, 4, 5, 6]

   ! Anura3D: [11, 22, 33, 12, 23, 31]  — positions 5 and 6 swapped
   integer, parameter :: ANURA3D_ORDER(6) = [1, 2, 3, 4, 6, 5]

   ! Abaqus:  [11, 22, 33, 12, 13, 23]  — same as internal
   integer, parameter :: ABAQUS_ORDER(6) = [1, 2, 3, 4, 5, 6]

contains

   ! ---------------------------------------------------------------------------
   ! Determine problem type from NDI and NSHR (as passed by the solver UMAT)
   ! ---------------------------------------------------------------------------
   function problem_type(ndi, nshr) result(ptype)
      integer, intent(in) :: ndi, nshr
      integer :: ptype

      if (ndi == 3 .and. nshr == 3) then
         ptype = PROBLEM_3D
      else if (ndi == 3 .and. nshr == 1) then
         ptype = PROBLEM_PLANE_STRAIN   ! or axisymmetric — identical inflation rule
      else if (ndi == 2 .and. nshr == 1) then
         ptype = PROBLEM_PLANE_STRESS
      else
         error stop "voigt_conventions: unrecognised NDI/NSHR combination"
      end if
   end function problem_type

   ! ---------------------------------------------------------------------------
   ! Inflate a solver-side Voigt vector (NTENS components) to internal 6-component
   !
   ! Plane strain / axisymmetric (NTENS=4):
   !   solver [11, 22, 33, 12]  ->  internal [11, 22, 33, 12, 0, 0]
   !   Components 5 (13) and 6 (23) are zero — no out-of-plane shear.
   !
   ! Plane stress (NTENS=3):
   !   solver [11, 22, 12]  ->  internal [11, 22, 0, 12, 0, 0]
   !   Component 3 (33) is set to zero here. For strain inflation the caller
   !   is responsible for inserting the eps_33 found by the plane-stress Newton
   !   iteration before passing to the integrator.
   !
   ! 3D (NTENS=6): identity (after convention reorder)
   ! ---------------------------------------------------------------------------
   function inflate(v_solver, ptype) result(v6)
      real(dp), intent(in) :: v_solver(:)   ! NTENS components, already in internal order
      integer,  intent(in) :: ptype
      real(dp) :: v6(6)

      v6 = 0.0_dp

      select case(ptype)
       case(PROBLEM_3D)
         v6 = v_solver(1:6)

       case(PROBLEM_PLANE_STRAIN, PROBLEM_AXISYMMETRIC)
         ! solver: [11, 22, 33, 12]  ->  internal positions 1,2,3,4
         v6(1) = v_solver(1)   ! 11
         v6(2) = v_solver(2)   ! 22
         v6(3) = v_solver(3)   ! 33
         v6(4) = v_solver(4)   ! 12
         ! v6(5) = 0  (13 — no out-of-plane shear)
         ! v6(6) = 0  (23 — no out-of-plane shear)

       case(PROBLEM_PLANE_STRESS)
         ! solver: [11, 22, 12]  ->  internal positions 1,2,4
         ! v6(3) = 0 for stress; for strain the wrapper must insert eps_33 afterwards
         v6(1) = v_solver(1)   ! 11
         v6(2) = v_solver(2)   ! 22
         v6(4) = v_solver(3)   ! 12
         ! v6(3) = 0  (sig_33 = 0 by definition; eps_33 set by caller)
         ! v6(5) = 0, v6(6) = 0

       case default
         error stop "inflate: unknown problem type"
      end select
   end function inflate

   ! ---------------------------------------------------------------------------
   ! Deflate an internal 6-component Voigt vector back to solver-side NTENS
   ! ---------------------------------------------------------------------------
   function deflate(v6, ptype) result(v_solver)
      real(dp), intent(in) :: v6(6)
      integer,  intent(in) :: ptype
      real(dp), allocatable :: v_solver(:)

      select case(ptype)
       case(PROBLEM_3D)
         allocate(v_solver(6))
         v_solver = v6

       case(PROBLEM_PLANE_STRAIN, PROBLEM_AXISYMMETRIC)
         allocate(v_solver(4))
         v_solver(1) = v6(1)   ! 11
         v_solver(2) = v6(2)   ! 22
         v_solver(3) = v6(3)   ! 33
         v_solver(4) = v6(4)   ! 12

       case(PROBLEM_PLANE_STRESS)
         allocate(v_solver(3))
         v_solver(1) = v6(1)   ! 11
         v_solver(2) = v6(2)   ! 22
         v_solver(3) = v6(4)   ! 12 (skips sig_33 — solver knows it is zero)

       case default
         error stop "deflate: unknown problem type"
      end select
   end function deflate

   ! ---------------------------------------------------------------------------
   ! Deflate a 6x6 internal stiffness matrix to NTENS x NTENS solver-side form
   ! ---------------------------------------------------------------------------
   function deflate_stiffness(C6, ptype) result(C_out)
      real(dp), intent(in) :: C6(6,6)
      integer,  intent(in) :: ptype
      real(dp), allocatable :: C_out(:,:)

      ! Index map: which internal indices appear in the solver output
      integer :: idx(6), n, i, j
      n = 0

      select case(ptype)
       case(PROBLEM_3D)
         n = 6; idx(1:6) = [1, 2, 3, 4, 5, 6]
       case(PROBLEM_PLANE_STRAIN, PROBLEM_AXISYMMETRIC)
         n = 4; idx(1:4) = [1, 2, 3, 4]
       case(PROBLEM_PLANE_STRESS)
         n = 3; idx(1:3) = [1, 2, 4]
       case default
         error stop "deflate_stiffness: unknown problem type"
      end select

      allocate(C_out(n, n))
      do i = 1, n
         do j = 1, n
            C_out(i,j) = C6(idx(i), idx(j))
         end do
      end do
   end function deflate_stiffness

   ! ---------------------------------------------------------------------------
   ! Convention reorder (unchanged from before)
   ! ---------------------------------------------------------------------------

   function to_internal(v, from_order) result(v_out)
      real(dp), intent(in) :: v(6)
      integer,  intent(in) :: from_order(6)
      real(dp) :: v_out(6)
      integer  :: i

      do i = 1, 6
         v_out(from_order(i)) = v(i)
      end do
   end function to_internal

   function from_internal(v, to_order) result(v_out)
      real(dp), intent(in) :: v(6)
      integer,  intent(in) :: to_order(6)
      real(dp) :: v_out(6)
      integer  :: i

      do i = 1, 6
         v_out(i) = v(to_order(i))
      end do
   end function from_internal

   function reorder_stiffness(C, from_order, to_order) result(C_out)
      real(dp), intent(in) :: C(6,6)
      integer,  intent(in) :: from_order(6), to_order(6)
      real(dp) :: C_out(6,6)
      real(dp) :: C_internal(6,6)
      integer  :: i, j

      do i = 1, 6
         do j = 1, 6
            C_internal(from_order(i), from_order(j)) = C(i,j)
         end do
      end do
      do i = 1, 6
         do j = 1, 6
            C_out(i,j) = C_internal(to_order(i), to_order(j))
         end do
      end do
   end function reorder_stiffness

end module mod_voigt_conventions
