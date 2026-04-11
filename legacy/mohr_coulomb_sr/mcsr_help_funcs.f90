!! **Legacy** — Low-level matrix/vector helpers for the MCSR model.
!!
!! Provides `MatVec` (matrix-vector product) and `DotProduct_2` (dot product).
!! These duplicate Fortran intrinsics (`matmul`, `dot_product`) and exist
!! only to support the legacy MCSR code. New code should use intrinsics directly.
module mod_SRMC_funcs
   use mod_csm_kinds, only: wp

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
      real(wp), intent(in)  :: xMat(N, N), Vec(N)
      integer, intent(in)          :: IM, N
      real(wp), intent(out) :: VecR(N)

      !***********************************************************************
      ! Local variables
      integer :: I, J
      real(wp) :: X

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
   Subroutine DotProduct_2(VecA, VecB,N, res)
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
      real(wp), intent(in)  :: VecA(N), VecB(N)
      integer, intent(in)          :: N
      real(wp), intent(out) :: res

      !***********************************************************************
      ! Local variables
      integer :: I
      res=0.0_wp
      Do I=1,N
         res=res+VecA(I)*VecB(I)
      end do

   end subroutine DotProduct_2


end module mod_SRMC_funcs
