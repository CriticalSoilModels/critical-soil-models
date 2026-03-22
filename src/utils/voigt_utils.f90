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
      integer(i32) :: i

      ! Initialize the deviatoric stress to the input stress
      dev_stress = stress

      ! Subtract mean stress from normal components
      do i = 1, 3
         dev_stress(i) = dev_stress(i) - mean_stress
      end do

   end function calc_dev_stress

   pure Subroutine calc_two_norm_tensor(Tensor, N, TwoNorm)
      !***********************************************************************
      !     ! WaveHello: This is actually the frobenius norm
      !     Calculate 2NormTensor = sqrt(A:A)
      !
      ! I   Tensor  : (Square or vector of dimension N)
      ! I   N     :   Number of elements
      ! O   2Norm : Resulting norm
      !
      !***********************************************************************
      implicit none
      real(kind = dp), intent(in) :: Tensor(N)
      integer, intent(in):: N

      real(kind = dp), intent(out) :: TwoNorm
      !***********************************************************************
      ! Local variables
      integer :: X, I

      X=N/2
      TwoNorm=0.0d0
      Do I=1,X
         TwoNorm=TwoNorm+Tensor(I)*Tensor(I)
      end Do
      Do I=X+1,N
         TwoNorm=TwoNorm+2*(Tensor(I)*Tensor(I))
      end do
      TwoNorm=sqrt(TwoNorm)

   end subroutine calc_two_norm_tensor

   pure Subroutine calc_two_norm_tensor_strain(Tensor, N, TwoNorm)
      !***********************************************************************
      !     ! WaveHello: This is actually the frobenius norm
      !     Calculate 2NormTensor = sqrt(A:A)
      !
      ! I   Tensor  : (Square or vector of dimension N)
      ! I   N     :   Number of elements
      ! O   2Norm : Resulting norm
      !
      !***********************************************************************
      implicit none
      real(kind = dp), intent(in)  :: Tensor(N)
      integer, intent(in) :: N

      real(kind = dp), intent(out) :: TwoNorm
      !***********************************************************************
      ! Local variables
      integer :: X, I

      X=N/2
      TwoNorm=0.0d0
      Do I=1,X
         TwoNorm=TwoNorm+Tensor(I)*Tensor(I)
      end Do
      Do I=X+1,N
         !The convention in UMAT is to use engineering shear strains
         TwoNorm=TwoNorm+0.5*(Tensor(I)*Tensor(I))
      end do
      TwoNorm=sqrt(TwoNorm)

   end subroutine calc_two_norm_tensor_strain

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

   pure Subroutine calc_tensor_inner_product(TensorA, TensorB, N, Re)
      !***********************************************************************
      !
      !     Calculate 2NormTensor = sqrt(A:A)
      !
      ! I   Tensor  : (Square or vector of dimension N)
      ! I   N     :   Number of elements
      ! O   2Norm : Resulting norm
      !
      !***********************************************************************
      implicit none

      real(kind = dp),intent(in) :: TensorA(N), TensorB(N)
      real(kind = dp),intent(out) :: Re
      integer, intent(in) :: N
      !***********************************************************************

      ! Local variables
      integer :: X, I

      X=N/2
      Re=0.0d0
      Do I=1,X
         Re=Re+TensorA(I)*TensorB(I)
      end Do
      Do I=X+1,N
         Re=Re+2*(TensorA(I)*TensorB(I))
      end do
   end subroutine calc_tensor_inner_product

end module mod_voigt_utils
