      MODULE mesh_generator

        CONTAINS

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... Creates a mesh of Gauss-Legendre points
! ... INPUT: x1 (x2) - Initial (final) points - n - # of points
! ... OUTPUT: (x,w) vectors of length n - points and weights, respectively
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE GAUSS(x1,x2,n,x,w)
        USE precision_definition
        IMPLICIT NONE
        REAL (long), INTENT(IN) :: x1,x2
        INTEGER (ilong), INTENT (IN) :: n
        REAL (long), INTENT(OUT), DIMENSION(:) :: x,w
        
        REAL (long) :: eps,p1,p2,p3,pp,xm,xl,z1,z,pi
        INTEGER (ilong) :: i,m,j

        pi=4.d0*datan(1.d0)

        eps=3e-14
        m=(n+1)/2
        xm=0.5d0*(x2+x1)
        xl=0.5d0*(x2-x1)
        z1=100000.d0

        do i=1,m
! ... Aproximation of the ith root
         z=cos(pi*(real(i,long)-0.25d0)/(real(n,long)+0.5d0))
! ... Main loop for finding the root with Newtons approximation         
 1       continue
         p1=1.d0
         p2=0.d0
         do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*real(j,long)-1.d0)*z*p2-(real(j,long)-1.d0)*p3)/real(j,long)
 11      enddo
! ... p1 is the Legendre polynomial 
! ... we compute pp, its derivative, by a relation involving
! ... p2, the polynomial of one lower order
         pp=n*(z*p1-p2)/(z*z-1.d0)
         z1=z
! ... Newtons method
         z=z1-p1/pp

         if(abs(z-z1) > eps) go to 1
         
         x(i)=xm-xl*z
         x(n+1-i)=xm+xl*z
         w(i)=2.d0*xl/((1-z*z)*pp*pp)
         w(n+1-i)=w(i)
      enddo

      END SUBROUTINE GAUSS

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... Creates a mesh of linearly spaced points
! ... INPUT: x1 (x2) - Initial (final) points - n - # of points
! ... OUTPUT: (x,w) vectors of length n - points and weights, respectively
! ... WEIGHTS AS FOR TRAPEZOIDAL RULE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE LINEAR(x1,x2,n,x,w)
        USE precision_definition
        IMPLICIT NONE
        REAL (long), INTENT(IN) :: x1,x2
        INTEGER (ilong), INTENT (IN) :: n
        REAL (long), INTENT(OUT), DIMENSION(:) :: x,w
        
        REAL (long) :: dx
        INTEGER (ilong) :: i

        dx=(x2-x1)/real(n-1,long)
        do i=1,N
           x(i) = x1+real(i-1,long)*dx
           w(i) = dx
           if(i == 1 .or. i == n) w(i) = dx/2d0
        enddo       
      END SUBROUTINE LINEAR

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... Creates a mesh of logarithmic points
! ... INPUT: x1 (x2) - Initial (final) points - n - # of points
! ... updw - CHARACTER TO KNOW WHETHER POINTS GO UPWARDS OR DOWNWARDS
! ... OUTPUT: (x,w) vectors of length n - points and weights, respectively
! ... WEIGHTS AS FOR TRAPEZOIDAL RULE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE LOGARITHMIC(x1,x2,n,x,w,updw)
        USE precision_definition
        IMPLICIT NONE
        REAL (long), INTENT(IN) :: x1,x2
        INTEGER (ilong), INTENT (IN) :: n
        REAL (long), INTENT(OUT), DIMENSION(:) :: x,w
        CHARACTER(len=1) :: updw

        REAL (long) :: fac,x0
        REAL (long), DIMENSION(n) :: xaux,waux
        INTEGER (ilong) :: i,iset

        call GAUSS(0d0,1d0,n,xaux,waux)

        fac= (x2-x1)/log( (1.d0+xaux(N))/(1.d0-xaux(N)) )
        if(updw == 'U') then
           x0=x1
        elseif(updw == 'D') then
           x0=x2
           fac=-fac
        else
           write(*,*) 'Uknown updw in LOGARITHMIC'
           stop
        endif

        do i=1,N
           if(updw == 'U') iset=i
           if(updw == 'D') iset=n-i+1
           x(iset) = x0 + fac*log( (1d0+xaux(i))/(1d0-xaux(i)) )
           w(iset) = waux(i)*abs(fac)*2d0/( 1d0 - xaux(i)**2 )
        enddo
   
      END SUBROUTINE LOGARITHMIC


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... Creates a mesh of logarithmic points
! ... INPUT: x1 (x2) - Initial (final) points - n - # of points
! ... updw - CHARACTER TO KNOW WHETHER POINTS GO UPWARDS OR DOWNWARDS
! ... OUTPUT: (x,w) vectors of length n - points and weights, respectively
! ... WEIGHTS AS FOR TRAPEZOIDAL RULE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE TANGENTIAL(x1,x2,n,x,w,updw)
        USE precision_definition
        IMPLICIT NONE
        REAL (long), INTENT(IN) :: x1,x2
        INTEGER (ilong), INTENT (IN) :: n
        REAL (long), INTENT(OUT), DIMENSION(:) :: x,w
        CHARACTER(len=1) :: updw

        REAL (long) :: fac,x0,pi
        REAL (long), DIMENSION(n) :: xaux,waux
        INTEGER (ilong) :: i,iset

        call GAUSS(0d0,1d0,n,xaux,waux)

        pi=4._long*atan(1._long)
        fac=(x2-x1)/tan( pi/2d0*xaux(N) )
        if(updw == 'U') then
           x0=x1
        elseif(updw == 'D') then
           x0=x2
           fac=-fac
        else
           write(*,*) 'Uknown updw in LOGARITHMIC'
           stop
        endif

        do i=1,N
           if(updw == 'U') iset=i
           if(updw == 'D') iset=n-i+1
           x(iset) = x0 + fac*tan( pi/2d0*xaux(i)  )
           w(iset) = waux(i)*abs(fac)*pi/2d0/(cos(pi/2d0*xaux(i))**2)
        enddo
   
      END SUBROUTINE TANGENTIAL

    END MODULE mesh_generator
