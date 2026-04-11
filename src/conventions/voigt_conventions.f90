!! Voigt convention translation and inflate/deflate utilities.
!!
!! The library's **internal** Voigt ordering is `[11, 22, 33, 12, 13, 23]`.
!! Convention translation and inflate/deflate happen **only** at the UMAT
!! boundary — never inside models or integrators.
!!
!! ### Supported problem types
!!
!! | Type            | NDI | NSHR | NTENS | Solver components          |
!! |-----------------|-----|------|-------|----------------------------|
!! | 3D              | 3   | 3    | 6     | [11, 22, 33, 12, 13, 23]   |
!! | Plane strain    | 3   | 1    | 4     | [11, 22, 33, 12]           |
!! | Axisymmetric    | 3   | 1    | 4     | [11, 22, 33, 12] (33=hoop) |
!! | Plane stress    | 2   | 1    | 3     | [11, 22, 12] (33=0)        |
!!
!! Internally every model always works with 6 components. `inflate` pads solver
!! vectors to 6 on entry; `deflate` strips them back on exit.
!! Plane stress requires a Newton iteration in the UMAT wrapper to satisfy
!! σ₃₃ = 0 — this is handled in the wrapper, not here.

module mod_voigt_conventions
   use mod_csm_kinds, only: wp
   implicit none

   ! ---------------------------------------------------------------------------
   ! Problem type identifiers
   ! ---------------------------------------------------------------------------
   integer, parameter :: PROBLEM_3D           = 1  !! Full 3D (NTENS=6)
   integer, parameter :: PROBLEM_PLANE_STRAIN = 2  !! Plane strain (NTENS=4)
   integer, parameter :: PROBLEM_AXISYMMETRIC = 3  !! Axisymmetric (NTENS=4)
   integer, parameter :: PROBLEM_PLANE_STRESS = 4  !! Plane stress (NTENS=3)

   ! ---------------------------------------------------------------------------
   ! Named reorder arrays for each supported solver convention.
   ! Internal (this library): [11, 22, 33, 12, 13, 23]
   ! ---------------------------------------------------------------------------

   integer, parameter :: INTERNAL_ORDER(6) = [1, 2, 3, 4, 5, 6]
   !! Identity permutation — internal library order.

   integer, parameter :: ANURA3D_ORDER(6) = [1, 2, 3, 4, 6, 5]
   !! Anura3D convention: `[11, 22, 33, 12, 23, 31]` — positions 5 and 6 swapped.

   integer, parameter :: ABAQUS_ORDER(6) = [1, 2, 3, 4, 5, 6]
   !! Abaqus convention: `[11, 22, 33, 12, 13, 23]` — same as internal.

contains

   function problem_type(ndi, nshr) result(ptype)
      !! Determine the problem type from the NDI and NSHR arguments passed by the solver.
      integer, intent(in) :: ndi   !! Number of direct stress components
      integer, intent(in) :: nshr  !! Number of shear stress components
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

   function inflate(v_solver, ptype) result(v6)
      !! Pad a solver-side Voigt vector (NTENS components) to the internal 6-component form.
      !!
      !! Missing components are set to zero:
      !! - Plane strain/axisymmetric: positions 5 (13) and 6 (23) → 0
      !! - Plane stress: position 3 (33) → 0 for stress;
      !!   for strain the caller inserts ε₃₃ after calling this function.
      real(wp), intent(in) :: v_solver(:)  !! NTENS components, already in internal order
      integer,  intent(in) :: ptype        !! Problem type (use `problem_type`)
      real(wp) :: v6(6)

      v6 = 0.0_wp

      select case(ptype)
       case(PROBLEM_3D)
         v6 = v_solver(1:6)

       case(PROBLEM_PLANE_STRAIN, PROBLEM_AXISYMMETRIC)
         v6(1) = v_solver(1)   ! 11
         v6(2) = v_solver(2)   ! 22
         v6(3) = v_solver(3)   ! 33
         v6(4) = v_solver(4)   ! 12

       case(PROBLEM_PLANE_STRESS)
         v6(1) = v_solver(1)   ! 11
         v6(2) = v_solver(2)   ! 22
         v6(4) = v_solver(3)   ! 12

       case default
         error stop "inflate: unknown problem type"
      end select
   end function inflate

   function deflate(v6, ptype) result(v_solver)
      !! Strip an internal 6-component Voigt vector back to the solver-side NTENS form.
      real(wp), intent(in) :: v6(6)          !! Internal 6-component vector
      integer,  intent(in) :: ptype          !! Problem type (use `problem_type`)
      real(wp), allocatable :: v_solver(:)

      select case(ptype)
       case(PROBLEM_3D)
         allocate(v_solver(6))
         v_solver = v6

       case(PROBLEM_PLANE_STRAIN, PROBLEM_AXISYMMETRIC)
         allocate(v_solver(4))
         v_solver(1) = v6(1)
         v_solver(2) = v6(2)
         v_solver(3) = v6(3)
         v_solver(4) = v6(4)

       case(PROBLEM_PLANE_STRESS)
         allocate(v_solver(3))
         v_solver(1) = v6(1)
         v_solver(2) = v6(2)
         v_solver(3) = v6(4)   ! 12 (skips σ₃₃ — solver knows it is zero)

       case default
         error stop "deflate: unknown problem type"
      end select
   end function deflate

   function deflate_stiffness(C6, ptype) result(C_out)
      !! Extract the NTENS×NTENS solver-side stiffness block from the internal 6×6 matrix.
      real(wp), intent(in) :: C6(6,6)          !! Internal 6×6 stiffness matrix
      integer,  intent(in) :: ptype            !! Problem type (use `problem_type`)
      real(wp), allocatable :: C_out(:,:)

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

   function to_internal(v, from_order) result(v_out)
      !! Reorder a Voigt vector from a solver convention to internal order.
      real(wp), intent(in) :: v(6)           !! Input vector in solver convention
      integer,  intent(in) :: from_order(6)  !! Solver order (e.g. `ANURA3D_ORDER`)
      real(wp) :: v_out(6)
      integer  :: i

      do i = 1, 6
         v_out(from_order(i)) = v(i)
      end do
   end function to_internal

   function from_internal(v, to_order) result(v_out)
      !! Reorder a Voigt vector from internal order to a solver convention.
      real(wp), intent(in) :: v(6)         !! Input vector in internal order
      integer,  intent(in) :: to_order(6)  !! Target order (e.g. `ANURA3D_ORDER`)
      real(wp) :: v_out(6)
      integer  :: i

      do i = 1, 6
         v_out(i) = v(to_order(i))
      end do
   end function from_internal

   function reorder_stiffness(C, from_order, to_order) result(C_out)
      !! Reorder a 6×6 stiffness matrix between two Voigt conventions.
      real(wp), intent(in) :: C(6,6)          !! Input stiffness in source convention
      integer,  intent(in) :: from_order(6)   !! Source order
      integer,  intent(in) :: to_order(6)     !! Target order
      real(wp) :: C_out(6,6)
      real(wp) :: C_internal(6,6)
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
