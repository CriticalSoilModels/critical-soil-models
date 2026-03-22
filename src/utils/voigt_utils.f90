module mod_voigt_utils
   use stdlib_kinds, only: dp, i32 => int32
   implicit none
   private
   public :: calc_dev_stress, calc_two_norm_tensor, calc_two_norm_tensor_strain, calc_voigt_to_matrix, calc_matrix_to_voigt, &
   calc_voigt_square, calc_tensor_inner_product
contains

   pure function calc_dev_stress(stress, mean_stress) result(dev_stress)
      real(kind = dp), intent(in) :: stress(6), mean_stress
      real(kind = dp) :: dev_stress(6)

      ! Local variables
      integer :: i

      ! Initialize the deviatoric stress to the input stress
      dev_stress = stress

      ! Subtract mean stress from normal components
      do i = 1, 3
         dev_stress(i) = dev_stress(i) - mean_stress
      end do

   end function calc_dev_stress

   pure function calc_two_norm_tensor(v) result(norm)
      ! Frobenius norm of a symmetric tensor in Voigt-6 form: sqrt(A:A)
      ! Normal components (1:3) enter once; shear components (4:6) enter twice.
      real(kind = dp), intent(in) :: v(6)
      real(kind = dp) :: norm
      integer :: i

      norm = 0.0_dp
      do i = 1, 3
         norm = norm + v(i)*v(i)
      end do
      do i = 4, 6
         norm = norm + 2.0_dp*(v(i)*v(i))
      end do
      norm = sqrt(norm)

   end function calc_two_norm_tensor

   pure function calc_two_norm_tensor_strain(v) result(norm)
      ! Frobenius norm of a symmetric strain tensor in Voigt-6 form.
      ! Uses engineering shear strains (factor 0.5 on shear components).
      real(kind = dp), intent(in) :: v(6)
      real(kind = dp) :: norm
      integer :: i

      norm = 0.0_dp
      do i = 1, 3
         norm = norm + v(i)*v(i)
      end do
      do i = 4, 6
         norm = norm + 0.5_dp*(v(i)*v(i))
      end do
      norm = sqrt(norm)

   end function calc_two_norm_tensor_strain

   pure function calc_voigt_to_matrix(voigt_vector) result(matrix)

      real(kind = dp), intent(in) :: voigt_vector(6)
      real(kind = dp) :: matrix(3, 3)

      ! Local variables
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

      real(kind = dp), intent(in) :: matrix(3,3)
      real(kind = dp) :: voigt_vector(6)

      ! Local variables
      integer :: i

      ! This works for stress voigt
      ! Warning: For strain need to half the shear terms before passing in
      ! Incremental driver Voigt order
      !={
      !    11 (xx),
      !    22 (yy),
      !    33 (zz),
      !    12 (xy),
      !    13 (xz),
      !    23 (yz)
      !}

      ! Fill the diagonal
      do i = 1, 3
         voigt_vector(i) = matrix(i,i)
      end do

      ! Store the shear terms
      voigt_vector(4) = matrix(1, 2)
      voigt_vector(5) = matrix(1, 3)
      voigt_vector(6) = matrix(2, 3)
   end function calc_matrix_to_voigt

   pure function calc_voigt_square(voigt) result(voigt2)
      ! Multies a voigt vector using matrix multiplication against itself
      ! Assumes incremental driver ordering for the voigt vector
      real(kind = dp), intent(in) :: voigt(6)
      real(kind = dp) :: voigt2(6)

      voigt2(1)=voigt(1)**2 + voigt(4)**2 + voigt(5)**2
      voigt2(2)=voigt(2)**2 + voigt(4)**2 + voigt(6)**2
      voigt2(3)=voigt(3)**2 + voigt(5)**2 + voigt(6)**2
      voigt2(4)=voigt(4) * ( voigt(1) + voigt(2) ) + voigt(5)*voigt(6)
      voigt2(5)=voigt(5) * ( voigt(1) + voigt(3) ) + voigt(4)*voigt(6)
      voigt2(6)=voigt(6) * ( voigt(2) + voigt(3) ) + voigt(4)*voigt(5)
   end function calc_voigt_square

   pure function calc_tensor_inner_product(a, b) result(inner)
      ! Symmetric tensor inner product A:B in Voigt-6 form.
      ! Normal components (1:3) enter once; shear components (4:6) enter twice.
      real(kind = dp), intent(in) :: a(6), b(6)
      real(kind = dp) :: inner
      integer :: i

      inner = 0.0_dp
      do i = 1, 3
         inner = inner + a(i)*b(i)
      end do
      do i = 4, 6
         inner = inner + 2.0_dp*(a(i)*b(i))
      end do
   end function calc_tensor_inner_product

end module mod_voigt_utils
