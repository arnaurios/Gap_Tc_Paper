      MODULE LEGENDRE_POLS
       USE precision_definition

       CONTAINS
!       ========================================================
!       Purpose: This program computes the Legendre polynomials
!                Pn(x) and their derivatives Pn'(x) using
!                subroutine LPN
!       Input :  x --- Argument of Pn(x)
!                n --- Degree of Pn(x) ( n = 0,1,...)
!       Output:  PN(n) --- Pn(x)
!                PD(n) --- Pn'(x)
!       Example:    x = 0.5
!                  n          Pn(x)            Pn'(x)
!                ---------------------------------------
!                  0       1.00000000        .00000000
!                  1        .50000000       1.00000000
!                  2       -.12500000       1.50000000
!                  3       -.43750000        .37500000
!                  4       -.28906250      -1.56250000
!                  5        .08984375      -2.22656250
!       ========================================================

       SUBROUTINE LPN(N,X,PN)
!
!       ===============================================
!       Purpose: Compute Legendre polynomials Pn(x)
!                and their derivatives Pn'(x)
!       Input :  x --- Argument of Pn(x)
!                n --- Degree of Pn(x) ( n = 0,1,...)
!       Output:  PN(n) --- Pn(x)
!                PD(n) --- Pn'(x)
!       ===============================================
          USE precision_definition
          IMPLICIT NONE
          REAL(long), INTENT(IN), DIMENSION(:) :: X
          INTEGER(ilong), INTENT(IN) :: N
          REAL(long), INTENT(OUT), DIMENSION(0:N,size(X)) :: PN

          REAL(long), DIMENSION(size(X)) :: P0,P1,PF
          INTEGER(ilong) :: K
          REAL(long) :: XK

          PN(0,:)=1.0D0
          
          if(N > 0) then
             PN(1,:)=X
             P0(:)=1.0D0
             P1(:)=X
             DO 10 K=2,N
                XK=real(K,long)
                PF=(2.0D0*XK-1.0D0)/XK*X*P1-(XK-1.0D0)/XK*P0
                PN(K,:)=PF
                P0=P1
                P1=PF
10           ENDDO
          endif
      END SUBROUTINE LPN

    END MODULE LEGENDRE_POLS
