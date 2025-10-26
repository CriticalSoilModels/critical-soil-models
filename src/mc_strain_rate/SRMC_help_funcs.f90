module mod_SRMC_funcs
   use kind_precision_module, only: real_type => dp

   implicit none

contains
   ! Instead of this fucntion you can use matmul. It's a fortran intrinsic function for multiplying matrices of any size
   Subroutine MatVec(xMat,IM,Vec,N,VecR)
      !***********************************************************************
      !
      !     Calculate VecR = xMat*Vec
      !
      ! I   xMat  : (Square) Matrix (IM,*)
      ! I   Vec   : Vector
      ! I   N     : Number of rows/colums
      ! O   VecR  : Resulting vector
      !
      !***********************************************************************
      implicit none
      real(real_type), intent(in)  :: xMat(N, N), Vec(N)
      integer, intent(in)          :: IM, N
      real(real_type), intent(out) :: VecR(N)

      !***********************************************************************
      ! Local variables
      integer :: I, J
      real(real_type) :: X

      Do I=1,N
         X=0
         Do J=1,N
            X=X+xMat(I,J)*Vec(J)
         End Do
         VecR(I)=X
      End Do
      Return
   End Subroutine MatVec

   ! Instead of this function you should use dot_product. It's a fortran instrinsic
   Subroutine DotProduct_2(VecA, VecB,N, Dp)
      !***********************************************************************
      !
      !     Calculate the dot product of A(Nx1) and B(1xN)
      !
      ! I   VecA VecB  : Vectors
      ! I   N     :   Dimension
      ! O   Dp : Dot product
      !
      !***********************************************************************
      implicit none
      real(real_type), intent(in)  :: VecA(N), VecB(N)
      integer, intent(in)          :: N
      real(real_type), intent(out) :: Dp

      !***********************************************************************
      ! Local variables
      integer :: I
      Dp=0.0d0
      Do I=1,N
         Dp=Dp+VecA(I)*VecB(I)
      end do

   end subroutine DotProduct_2


end module mod_SRMC_funcs
