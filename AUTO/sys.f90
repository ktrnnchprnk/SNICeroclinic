!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!          SNICeroclinic system Oct2024
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION x,y,a,b,c,funcx,funcy,eps
     
      
      x = U(1)
      y = U(2)
      
      a = PAR(1)
      b = PAR(2)
      c = PAR(3)
      eps = PAR(4)
      
      funcy =y**3-3*y
      funcx = a*x**2+b*x+c
      F(1) = eps*(funcx-y)
      F(2) = (x-funcy)
      
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
      !     ----------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
       ! Setup that works for central SNICeroclinic
       PAR(1)=0.042 !a
       PAR(2)=0.2 !b
       PAR(3)=-1.8 !c
       PAR(4)=1 !epsilon
       U(1)=-2.35127
       U(2)=-2.038058
       


       
       END SUBROUTINE STPNT

      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
	!periodic
       FB(1)=U0(1)-U1(1) 
       FB(2)=U0(2)-U1(2)
       

      END SUBROUTINE BCND

      SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
!      ---------- ----

       IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
       INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
       DOUBLE PRECISION PAR(*)
       DOUBLE PRECISION U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
       DOUBLE PRECISION FI(NINT)
       DOUBLE PRECISION DINT(NINT,*)

       	fi(1) = UPOLD(1)*(U(1)-UOLD(1))+UPOLD(2)*(U(2)-UOLD(2))		! phase condition for u

       !el2=u(1)*U(1)+U(2)*U(2)+U(3)*U(3)+U(4)*U(4)
       !FI(1)=(1/PAR(11)*(el2)-PAR(9))!a^2=parametro

       RETURN
       END 

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT
      SUBROUTINE PVLS(NDIM,U,PAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

      DOUBLE PRECISION, EXTERNAL :: GETP,GETU2
      INTEGER NDX,NCOL,NTST
      PAR(5)=GETP('STA',0,U)
      PAR(6)=GETP('MIN',2,U)

     
      END SUBROUTINE PVLS
