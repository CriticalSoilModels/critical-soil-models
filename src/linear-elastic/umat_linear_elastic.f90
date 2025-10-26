module mod_lin_elastic
   use stdlib_kinds, only: dp
   implicit none
   private
   public :: UMAT

contains
   SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD, &
      RPL,DDSDDT,DRPLDE,DRPLDT, &
      STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME, &
      NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT, &
      CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      implicit none

      CHARACTER*80 CMNAME
      integer :: NTENS,NSTATEV,NPROPS,NDI,NSHR,NOEL, &
         NPT,LAYER,KSPT,KSTEP,KINC
      real(8) :: SSE,SPD,SCD,RPL,DRPLDT,DTIME,TEMP,DTEMP, &
         PNEWDT,CELENT
      real(8) :: STRESS(NTENS),STATEV(NSTATEV), &
         DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS), &
         STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1), &
         PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

      ! local variables
      real(dp), PARAMETER :: ONE=1.0D0, TWO=2.0D0
      integer :: I,J
      real(dp) :: E,ANU,ALAMDA,AMU

      E=PROPS(1)
      ANU=PROPS(2)
      ALAMDA=ANU*E/ (ONE+ANU)/(ONE-TWO*ANU)
      AMU=E/TWO/(ONE+ANU)
      DO I=1,NTENS
         DO J=1,NTENS
            DDSDDE(I,J)=0.0D0
         ENDDO
      ENDDO
      DDSDDE(1,1)=ALAMDA+TWO*AMU
      DDSDDE(2,2)=DDSDDE(1,1)
      DDSDDE(3,3)=DDSDDE(1,1)
      DDSDDE(4,4)=AMU
      DDSDDE(5,5)=AMU
      DDSDDE(6,6)=AMU
      DDSDDE(1,2)=ALAMDA
      DDSDDE(1,3)=ALAMDA
      DDSDDE(2,3)=ALAMDA
      DDSDDE(2,1)=DDSDDE(1,2)
      DDSDDE(3,1)=DDSDDE(1,3)
      DDSDDE(3,2)=DDSDDE(2,3)

      DO I=1,NTENS
         DO J=1,NTENS
            STRESS(I)=STRESS(I)+DDSDDE(I,J)*DSTRAN(J)
         ENDDO
      ENDDO
      RETURN
   END
end module mod_lin_elastic

