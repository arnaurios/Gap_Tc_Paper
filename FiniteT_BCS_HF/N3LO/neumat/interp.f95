! ... MODULE THAT INTERFACES ALL NEEDED DATA FOR INTERPOLATION
      MODULE interpolation
        CONTAINS
! ... ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... !D LINEAR INTERPOLATION
      PURE SUBROUTINE LIN_INT(xa,ya,x,y)
        USE precision_definition
!        USE interpolation
        IMPLICIT NONE
        REAL (long), INTENT(IN), DIMENSION(:) :: xa,ya
        REAL (long), INTENT(IN), DIMENSION(:) :: x
        REAL (long), INTENT(OUT), DIMENSION(size(x)) :: y
        INTEGER (ilong) :: ii,nf
        INTEGER (ilong), DIMENSION(size(x)) :: klo,khi

        nf=size(x)
        FORALL( ii=1:nf )
           klo(ii)=locate(xa,x(ii))
           khi(ii)=klo(ii)+1
           y(ii) = ya(klo(ii)) + ( ya(khi(ii))-ya(klo(ii)) )/( xa(khi(ii))-xa(klo(ii)) )*( x(ii)-xa(klo(ii)) )
        END FORALL        
      END SUBROUTINE LIN_INT

! ... ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... 2D LINEAR INTERPOLATION
      PURE SUBROUTINE LIN_INT2D(xa,ya,za,x,y,z)

        USE precision_definition
!        USE interpolation

        IMPLICIT NONE
        REAL (long), INTENT(IN), DIMENSION(:) :: xa,ya
        REAL (long), INTENT(IN), DIMENSION(:,:)  :: za
        REAL (long), INTENT(IN), DIMENSION(:) :: x,y
        REAL (long), INTENT(OUT), DIMENSION(size(x),size(y)) :: z

        INTEGER (ilong) :: ix,nsx,iy,nsy
        INTEGER (ilong), DIMENSION(size(x),size(y)) :: kxlo,kxhi,kylo,kyhi
        
        REAL (long), DIMENSION(size(x),size(y)) :: z1,z2,z3,z4,t,u

        nsx=size(x)
        nsy=size(y)

        FORALL( iy=1:nsy, ix=1:nsx) 
           kylo(ix,iy)=locate(ya,y(iy))
           kyhi(ix,iy)=kylo(ix,iy)+1
           u(ix,iy)=(y(iy)-ya(kylo(ix,iy)))/(ya(kyhi(ix,iy))-ya(kylo(ix,iy)))
           
           kxlo(ix,iy)=locate(xa,x(ix))
           kxhi(ix,iy)=kxlo(ix,iy)+1
           t(ix,iy)=(x(ix)-xa(kxlo(ix,iy)))/(xa(kxhi(ix,iy))-xa(kxlo(ix,iy)))

           z1(ix,iy)=za(kxlo(ix,iy),kylo(ix,iy))
           z2(ix,iy)=za(kxhi(ix,iy),kylo(ix,iy))
           z3(ix,iy)=za(kxhi(ix,iy),kyhi(ix,iy))
           z4(ix,iy)=za(kxlo(ix,iy),kyhi(ix,iy))

           z(ix,iy)=(1.d0-t(ix,iy))*(1.d0-u(ix,iy))*z1(ix,iy) &
                + t(ix,iy)*(1.d0-u(ix,iy))*z2(ix,iy) + t(ix,iy)*u(ix,iy)*z3(ix,iy) &
                + (1.d0-t(ix,iy))*u(ix,iy)*z4(ix,iy)
        END FORALL

      END SUBROUTINE LIN_INT2D

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... Given an array xx(1:N), and given a value x, returns a value j such that x is between
! ... xx(j) and xx(j + 1). xx must be monotonic, either increasing or decreasing. j = 0 or
! ... j = N-1 is returned to indicate that x is out of range.
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      INTEGER(ilong) PURE FUNCTION LOCATE(xx,x)
        use precision_definition
        IMPLICIT NONE
        REAL (long), INTENT(IN), DIMENSION(:) :: xx
        REAL (long), INTENT(IN) :: x
        INTEGER (ilong) :: n,jl,jm,ju
        LOGICAL :: ascnd

        n=size(xx)

        if( x < xx(1) ) then
           locate=1
        elseif( x > xx(n) ) then
           locate=n
        else
           ascnd = (xx(n) >= xx(1) ) ! True if ascending, false if descending
           jl=0
           ju=n+1
           do while (ju - jl > 1) 
              jm=(ju+jl)/2
              if (ascnd .eqv. (x >= xx(jm))) then
                 jl=jm
              else
                 ju=jm
              end if
              locate=jl
           end do
        end if

        if(locate==n) locate=n-1

      END FUNCTION LOCATE

!     ******************************************************************
!     CUBIC SPLINE, STARTING WITH ZERO SECOND DERIVATIVES AT THE        
!     BOUNDARIES OF THE APPROXIMATION INTERVAL.                         
!     IGO = F      BUILD UP SPLINE AND INTERPOLATE FOR M POINTS.        
!     IGO = D      BUILD UP SPLINE AND COMPUTE DERIVATIVES AT M POINTS. 
!     REAL*8 VERSION.        J.GALONSKA, 15.12.1971           
!     ******************************************************************
      PURE SUBROUTINE SPL_INT(X,Y,XI,FI,IGO)                            
        USE precision_definition
        IMPLICIT NONE
        REAL (long), INTENT(IN), DIMENSION(:) :: X,Y
        REAL (long), INTENT(IN), DIMENSION(:) :: XI
        REAL (long), INTENT(OUT), DIMENSION(size(XI)) :: FI
        CHARACTER*1, INTENT(IN) :: IGO
!                                   
        REAL (long), DIMENSION(size(X)) :: AU,Q
        REAL(long) :: FACT,HK,YK,YSAVE,AUX,HX,DIVQ
        REAL(long) :: DIJ,DIM1J,HI,HI2,DIJ3,DIM1J3
        INTEGER (ilong) :: N,NN,NN2,K,M,M1,M2,M3,J,KK
!
        N=size(x)
        M=size(xi)
!
        FACT=1d0/6d0
!  
        AU(1) = 0_long                                                      
        AU(N) = 0_long
        Q(1) = 0_long
        HK = X(2) - X(1)                                                  
        YSAVE = (Y(2)-Y(1)) / HK                                          
        AUX = 0_long                                     
        NN = N - 1                                                        
        DO 10  K = 2,NN                                                   
           HX = X(K+1) - X(K-1)                                              
           DIVQ = (HK*Q(K-1)+HX+HX)                                          
           HK = X(K+1) - X(K)                                                
           YK = (Y(K+1)-Y(K)) / HK                                           
           Q(K) = - HK / DIVQ                                                
           AU(K) = (6_long*(YK-YSAVE)-AUX) / DIVQ                               
           YSAVE = YK                                                        
           AUX = AU(K) * HK  
10      ENDDO
!                                                           
        NN2 = NN + 2                                                      
        DO 20  KK = 2,NN                                                  
           K = NN2 - KK                                                      
           AU(K) = Q(K) * AU(K+1) + AU(K)                                    
20      ENDDO
!     ****************************************************************** 
!     INTERPOLATION OR COMPUTATION OF DERIVATIVES. 
!     IGO = 1      INTERPOLATE FOR M POINTS.                            
!     IGO = 2      COMPUTE DERIVATIVES AT M POINTS.                    
!     ******************************************************************                              
        DO 100 J = 1,M                                                   
           IF (X(1) > XI(J))  THEN                                         
              M1=1
              M2=2
              GO TO 90
           ENDIF
           IF (XI(J) > X(N))  THEN                                       
              M1=N-1
              M2=N
              GO TO 90
           ENDIF
           M1 = 1   
           M2 = N         
50         M3 = (M2+M1)/2  
           IF (XI(J) >= X(M3))  GO TO 70 
           M2 = M3     
           GO TO 80  
70         M1 = M3   
80         IF (M1+1-M2 /= 0)  GO TO 50    
90         DIJ = X(M2) - XI(J)
           DIM1J = X(M1) - XI(J) 
           HI = X(M2) - X(M1) 
           HI2 = HI * HI
           IF (IGO == 'D' )  GO TO 95
           DIJ3 = DIJ * DIJ * DIJ 
           DIM1J3 = DIM1J * DIM1J * DIM1J
           FI(J) = FACT * (AU(M1)*DIJ3-AU(M2)*DIM1J3+(6_long*Y(M1)-HI2*AU(M1)) &
                *DIJ-(6_long*Y(M2)-HI2*AU(M2))*DIM1J) / HI  
           GO TO 100                                                         
95         FI(J) = FACT * (3_long*(AU(M2)*DIM1J*DIM1J-AU(M1)*DIJ*DIJ) &        
                -6_long*(Y(M1)-Y(M2))+HI2*(AU(M1)-AU(M2))) / HI               
100     END DO
      END SUBROUTINE SPL_INT


      END MODULE interpolation


