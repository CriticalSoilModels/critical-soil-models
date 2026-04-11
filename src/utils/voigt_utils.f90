!! Low-level Voigt vector operations (norms, tensor products, matrix↔Voigt conversion).
!!
!! All routines assume the internal Voigt ordering `[11, 22, 33, 12, 13, 23]`.
!! Shear components occupy positions 4–6; the symmetric doubling factor (×2 for
!! stress, ×0.5 for engineering strain) is applied explicitly in each function.

module mod_voigt_utils
   use mod_csm_kinds, only: wp
   implicit none
   private
   public :: calc_dev_stress, calc_two_norm_tensor, calc_two_norm_tensor_strain, &
             calc_voigt_to_matrix, calc_matrix_to_voigt, &
             calc_voigt_square, calc_tensor_inner_product
contains

   pure function calc_dev_stress(stress, mean_stress) result(dev_stress)
      !! Compute the deviatoric stress tensor: s = σ - p I.
      real(kind=wp), intent(in) :: stress(6)       !! Stress vector (Voigt)
      real(kind=wp), intent(in) :: mean_stress      !! Mean stress p [same units]
      real(kind=wp) :: dev_stress(6)

      integer :: i

      dev_stress = stress
      do i = 1, 3
         dev_stress(i) = dev_stress(i) - mean_stress
      end do

   end function calc_dev_stress

   pure function calc_two_norm_tensor(v) result(norm)
      !! Frobenius norm of a symmetric stress tensor in Voigt-6 form: √(A:A).
      !!
      !! Normal components (1:3) enter once; shear components (4:6) enter twice
      !! to account for the off-diagonal symmetry.
      real(kind=wp), intent(in) :: v(6)  !! Voigt stress vector
      real(kind=wp) :: norm
      integer :: i

      norm = 0.0_wp
      do i = 1, 3
         norm = norm + v(i)*v(i)
      end do
      do i = 4, 6
         norm = norm + 2.0_wp*(v(i)*v(i))
      end do
      norm = sqrt(norm)

   end function calc_two_norm_tensor

   pure function calc_two_norm_tensor_strain(v) result(norm)
      !! Frobenius norm of a symmetric strain tensor stored with engineering shear strains.
      !!
      !! Shear components (4:6) store 2εᵢⱼ (engineering convention), so a factor
      !! of 0.5 is applied before squaring to recover the tensorial values.
      real(kind=wp), intent(in) :: v(6)  !! Voigt strain vector (engineering shear)
      real(kind=wp) :: norm
      integer :: i

      norm = 0.0_wp
      do i = 1, 3
         norm = norm + v(i)*v(i)
      end do
      do i = 4, 6
         norm = norm + 0.5_wp*(v(i)*v(i))
      end do
      norm = sqrt(norm)

   end function calc_two_norm_tensor_strain

   pure function calc_voigt_to_matrix(voigt_vector) result(matrix)
      !! Expand a Voigt-6 stress vector to a full 3×3 symmetric matrix.
      !!
      !! Voigt order: `[11, 22, 33, 12, 13, 23]`
      real(kind=wp), intent(in) :: voigt_vector(6)  !! Voigt stress vector
      real(kind=wp) :: matrix(3, 3)

      integer :: i

      do i = 1, 3
         matrix(i,i) = voigt_vector(i)
      end do

      matrix(1, 2) = voigt_vector(4)
      matrix(2, 1) = voigt_vector(4)
      matrix(1, 3) = voigt_vector(5)
      matrix(3, 1) = voigt_vector(5)
      matrix(2, 3) = voigt_vector(6)
      matrix(3, 2) = voigt_vector(6)
   end function calc_voigt_to_matrix

   pure function calc_matrix_to_voigt(matrix) result(voigt_vector)
      !! Compress a symmetric 3×3 matrix to a Voigt-6 stress vector.
      !!
      !! Voigt order: `[11, 22, 33, 12, 13, 23]`
      !!
      !! **Note:** for strain tensors, halve the shear components before calling
      !! this function since it stores raw matrix entries (not engineering shear).
      real(kind=wp), intent(in) :: matrix(3,3)       !! Symmetric 3×3 tensor
      real(kind=wp) :: voigt_vector(6)

      integer :: i

      do i = 1, 3
         voigt_vector(i) = matrix(i,i)
      end do

      voigt_vector(4) = matrix(1, 2)
      voigt_vector(5) = matrix(1, 3)
      voigt_vector(6) = matrix(2, 3)
   end function calc_matrix_to_voigt

   pure function calc_voigt_square(voigt) result(voigt2)
      !! Compute the matrix square A² and return it in Voigt-6 form.
      !!
      !! Equivalent to contracting A with itself: (A²)ᵢⱼ = Aᵢₖ Aₖⱼ.
      !! Assumes internal Voigt ordering `[11, 22, 33, 12, 13, 23]`.
      real(kind=wp), intent(in) :: voigt(6)   !! Input tensor (Voigt)
      real(kind=wp) :: voigt2(6)

      voigt2(1) = voigt(1)**2 + voigt(4)**2 + voigt(5)**2
      voigt2(2) = voigt(2)**2 + voigt(4)**2 + voigt(6)**2
      voigt2(3) = voigt(3)**2 + voigt(5)**2 + voigt(6)**2
      voigt2(4) = voigt(4) * (voigt(1) + voigt(2)) + voigt(5)*voigt(6)
      voigt2(5) = voigt(5) * (voigt(1) + voigt(3)) + voigt(4)*voigt(6)
      voigt2(6) = voigt(6) * (voigt(2) + voigt(3)) + voigt(4)*voigt(5)
   end function calc_voigt_square

   pure function calc_tensor_inner_product(a, b) result(inner)
      !! Symmetric tensor double contraction A:B in Voigt-6 form.
      !!
      !! Normal components (1:3) contribute once; shear components (4:6) contribute
      !! twice to account for both off-diagonal pairs.
      real(kind=wp), intent(in) :: a(6)   !! First tensor (Voigt)
      real(kind=wp), intent(in) :: b(6)   !! Second tensor (Voigt)
      real(kind=wp) :: inner
      integer :: i

      inner = 0.0_wp
      do i = 1, 3
         inner = inner + a(i)*b(i)
      end do
      do i = 4, 6
         inner = inner + 2.0_wp*(a(i)*b(i))
      end do
   end function calc_tensor_inner_product

end module mod_voigt_utils
