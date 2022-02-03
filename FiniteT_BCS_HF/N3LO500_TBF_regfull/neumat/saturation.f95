! ... COMPUTES SATURATION PROPERTIES
    SUBROUTINE SATURATION(rhosat)
      USE precision_definition
      USE physical_constants
      USE gogny_force
      USE gognyfunctions

      IMPLICIT NONE
      REAL(long) :: r0,xkf,xkf2,tkin,xk0,xq0
      REAL(long) :: efm,ep0,eps,et0,fdf,ff,r1,vcomp,vdpot,vpot,hmm,dukn
      REAL(long) :: vc2,vq
      REAL(long), INTENT(OUT) :: rhosat
      REAL (long), DIMENSION(2) :: xf,xfm
        
      INTEGER(ilong) :: ii

      hmm=(sum(htm)/2d0)/2d0

! ... INITIAL GUESS      
      r0=0.14d0
      do ii=1,20
         xkf=(6.d0*pi**2/4.d0*r0)**(1.d0/3.d0)
         xkf2=xkf**2
         tkin=hmm*xkf2

! ... Pressure/rho function and derivative
         call v_pot(r0,ep0,vpot,vdpot,vcomp,vc2,vq)

! ... Pressure/rho
         ff=2.d0/5.d0*tkin+vpot

! ... D(Pressure/rho)/D rho
         fdf=4.d0/15.d0*tkin/r0+vdpot

         r1=r0
         r0=r0-ff/fdf

         eps=dabs(r1-r0)/r0
      enddo

      xf=xkf*mu
      call f_effective_mass(xf,xf,xfm)
      dukn=-1d0/2d0/sqrt(pi)*sum( ( (W+2d0*B-2d0*H-4d0*M)*xfm )*mu**2 )
      efm=1d0/(1d0 + dukn/hmm/2d0)                          

! ... SATURATION ENERGY
      et0=3.d0/5.d0*tkin + ep0
! ... COMPRESSIBILITY
      xk0=9.d0*2.d0/3.d0*hmm*xkf**2+vcomp
! ... ANOTHER COMPRESSIBILITY
      xk0=-6d0/5d0*hmm*xkf**2 + vc2

! ... SKEWNESS
      xq0=24d0/5d0*hmm*xkf**2 + vq

      if(eps.gt.1e-6) then
         write(*,*) 'No convergence in saturation?'
         stop
      endif 

 909  format('Saturation properties' &
           ,/,2x,'r0=',f8.4,' [fm-3]  e0=',f9.4,' [MeV]  K0=' &
           ,f9.4,' [MeV]  m*/m=',f6.4,'  Q=',f9.4,' [MeV]')
      write(*,909) r0,et0,xk0,efm,xq0
      rhosat=r0

    END SUBROUTINE SATURATION

     
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     POTENTIAL ENERGY FOR SYMMETRIC NUCLEAR MATTER WITH GOGNY
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE V_POT(rr,vv,vd,v2d,vc,vc2,vq)
        USE precision_definition
        USE physical_constants
        USE gogny_force
        USE gognyfunctions

        IMPLICIT NONE

        REAL(long), INTENT(IN) :: rr
        REAL(long), INTENT(OUT) :: vv,vd,v2d,vc,vc2,vq
        
        REAL(long) :: xkf,ezr,edf,eef,dzr,ddf,def
        REAL(long) :: d2zr,d2df,d2ef,c2zr,c2df,c2ef,q2zr,q2df,q2ef


        REAL (long), DIMENSION(2) :: T0NN,AAN,BBN,xxx,ffg,d1g,d2g,d3g

        xkf=(6.d0*pi**2*rr/4d0)**(1.d0/3.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... PRESSURE OVER DENSITY - DIFFERRENT COMPONENTS
! ... ZERO RANGE TERM
        T0NN= 3d0*t0/8d0
        ezr=sum(T0NN*rr**(gamma+1d0))

! ... DIRECT FINITE RANGE
        AAN= pi**1.5d0*mu**3*( 4d0*W + 2d0*B - 2d0*H - M )/4d0
        edf=sum( AAN )*rr/2d0

! ... EXCHANGE FINITE RANGE
        xxx=mu*xkf
        call potential(xxx,ffg)
        call derpotential(xxx,d1g,d2g,d3g)

        BBN= -(W + B*2d0 - H*2d0 - M*4d0 )/sqrt(pi)
        eef=sum( BBN*ffg )/2d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... DERIVATIVE OF POTENTIAL ENERGY TIMES THE DENSITY
        dzr=sum((gamma+1d0)*T0NN*rr**(gamma+1d0))
        ddf=edf
        def=sum( BBN*d1g*xxx )/2d0/3d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... DERIVATIVE OF DERIVATIVE OF POTENTIAL ENERGY TIMES THE DENSITY
        d2zr=sum((gamma+1d0)**2*T0NN*rr**(gamma+1d0))/rr
        d2df=edf/rr
        d2ef=sum( BBN*xxx*(d1g + xxx*d2g) )/2d0/9d0/rr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... COMPRESSIBILITY
        c2zr=9d0*sum((gamma+1d0)*gamma*T0NN*rr**(gamma+1d0))
        c2df=0d0
        c2ef=sum( BBN*xxx*(-2d0*d1g + xxx*d2g) )/2d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... SKEWNESS
        q2zr=27d0*sum((gamma+1d0)*gamma*(gamma-1d0)*T0NN*rr**(gamma+1d0))
        q2df=0d0
        q2ef=sum( BBN*xxx*(10d0*d1g - 6d0*xxx*d2g + xxx**2*d3g) )/2d0

        vv=ezr+edf+eef
        vd=dzr+ddf+def
        v2d=d2zr+d2df+d2ef

        vc=9.d0*(vd+v2d*rr)
        vc2=c2zr+c2df+c2ef
        vq =q2zr+q2df+q2ef

    END SUBROUTINE V_POT

