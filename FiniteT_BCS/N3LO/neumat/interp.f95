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
! ... 1D SPLINE INTERPOLATION
! ... THESE INTERPOLATES BY FINDING CUBIC POLYNOMIALS THAT
! ... REPRODUCE THE VALUES OF Y AT EACH X VALUE
! ... AND HAVE CONTINUOUS DERIVATIVES AT THESE POINTS
! ... cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE SPL_INT(xa,ya,x,y)
        USE precision_definition

        IMPLICIT NONE
        REAL (long), INTENT(IN), DIMENSION(:) :: xa,ya
        REAL (long), INTENT(IN), DIMENSION(:) :: x
        REAL (long), INTENT(OUT), DIMENSION(size(x)) :: y
        REAL (long), DIMENSION(size(xa)) :: y2,u
        REAL (long) :: sig,p,qn,un,a,b,h
        INTEGER (ilong) :: ii,nf,nx
        INTEGER (ilong), DIMENSION(size(x)) :: klo,khi

        nf=size(xa)
        nx=size(x)

        y2(1)=0_long
        u(1)=0_long

        DO 11 ii=2,nf-1
           sig=(xa(ii)-xa(ii-1))/(xa(ii+1)-xa(ii-1))
           p=sig*y2(ii-1)+2_long
           y2(ii)=(sig-1_long)/p
           u(ii)=(6_long*((ya(ii+1)-ya(ii))/(xa(ii+1)-xa(ii))-(ya(ii)-ya(ii-1)) &
                /(xa(ii)-xa(ii-1)))/(xa(ii+1)-xa(ii-1))-sig*u(ii-1))/p
11      ENDDO

        qn=0_long
        un=0_long

        y2(nf)=(un-qn*u(nf-1))/(qn*y2(nf-1)+1_long)
! ... INTRODUCTION OF FORALL BELOW LEADS TO ERRORS - CHECK LOOP STRUCTURE
        DO ii=nf-1,1,-1
          y2(ii)=y2(ii)*y2(ii+1)+u(ii)
        ENDDO

        DO 13 ii=1,nx
           klo(ii)=locate(xa,x(ii))
           khi(ii)=klo(ii)+1
           h=xa(khi(ii))-xa(klo(ii))
           if(h==0.d0) then
              write(*,*) 'bad xa input'
              stop
           endif
           a=( xa(khi(ii))-x(ii) )/h
           b=( x(ii)-xa(klo(ii)) )/h
           y(ii)=a*ya( klo(ii) )+b*ya( khi(ii) )+ &
                ( (a**3-a)*y2( klo(ii) )+(b**3-b)*y2( khi(ii) ) )*(h**2)/6_long
13      ENDDO

      END SUBROUTINE SPL_INT


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

           z(ix,iy)=(1_long-t(ix,iy))*(1_long-u(ix,iy))*z1(ix,iy) &
                + t(ix,iy)*(1_long-u(ix,iy))*z2(ix,iy) + t(ix,iy)*u(ix,iy)*z3(ix,iy) &
                + (1_long-t(ix,iy))*u(ix,iy)*z4(ix,iy)
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

        locate=1
        if( x < xx(1) ) then
           locate=1
        elseif( x > xx(n) ) then
           locate=n-1
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
      END FUNCTION LOCATE


      END MODULE interpolation
