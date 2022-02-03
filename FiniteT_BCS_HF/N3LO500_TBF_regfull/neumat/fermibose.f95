MODULE fermibose
  USE precision_definition

  CONTAINS

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... FERMI-DIRAC FUNCTION WITH PROPER LIMITS 
      FUNCTION fermi(t,xmu,x)   
        USE precision_definition
        REAL(long) :: t,xmu

        REAL(long), INTENT(IN), ALLOCATABLE, DIMENSION(:) :: x
        REAL(long), DIMENSION(size(x)) :: z
        REAL(long), DIMENSION(size(x)) :: fermi

        z=(x-xmu)/t

        where (z > 60._long) 
           fermi=0._long
        elsewhere (z < -60.d0) 
           fermi=1._long
        elsewhere
           fermi=1._long/(1._long+dexp(z))
        endwhere
        
      END FUNCTION fermi

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... FERMI-DIRAC FUNCTION WITH PROPER LIMITS 
      FUNCTION ttanh(t,x)   
        USE precision_definition
        REAL(long) :: t

        REAL(long), INTENT(IN), ALLOCATABLE, DIMENSION(:) :: x
        REAL(long), DIMENSION(size(x)) :: z
        REAL(long), DIMENSION(size(x)) :: ttanh

        z=x/t/2._long

        where (z > 30._long) 
           ttanh= 1._long
        elsewhere (z < -30.d0) 
           ttanh=-1._long
        elsewhere
           ttanh=(dexp(2._long*z) - 1._long)/( dexp(2._long*z) + 1._long)
        endwhere
        
      END FUNCTION ttanh

      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... BOSE-EINSTEN FUNCTION WITH PROPER LIMITS 
      FUNCTION bose(t,xmu,x)
        USE precision_definition
        REAL(long) :: t,xmu
        REAL(long), INTENT(IN), ALLOCATABLE, DIMENSION(:) :: x
        REAL(long), DIMENSION(size(x)) :: z
        REAL(long), DIMENSION(size(x)) :: bose

        z=(x-2._long*xmu)/t

        where (z > 60._long) 
           bose=0._long
        elsewhere (z < -60._long) 
           bose=-1._long
        elsewhere
           bose=1._long/(dexp(z)-1._long)
        endwhere

      END FUNCTION bose

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... DERIVATIVE OF FERMI-DIRAC FUNCTION WITH PROPER LIMITS
      FUNCTION derfermiT(t,xmu,x)   
        USE precision_definition
        REAL(long) :: t,xmu

        REAL(long), INTENT(IN), ALLOCATABLE, DIMENSION(:) :: x
        REAL(long), DIMENSION(size(x)) :: z
        REAL(long), DIMENSION(size(x)) :: derfermiT

        z=(x-xmu)/t

        where(dabs(z) > 60._long) 
           derfermiT=0_long
        elsewhere
           derfermiT=z*dexp(z)/(1._long+dexp(z))**2/t
        endwhere

      END FUNCTION derfermiT

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... SIGMA FUNCTION WITH PROPER LIMITS 
! ... THIS IS USED IN ENTROPY CALCULATIONS OF FERMIONIC SYSTEMS
      FUNCTION sigma(t,xmu,x)
        USE precision_definition
        REAL(long) :: t,xmu

        REAL(long), INTENT(IN), ALLOCATABLE, DIMENSION(:) :: x
        REAL(long), DIMENSION(size(x)) :: z,ff,s1,s2
        REAL(long), DIMENSION(size(x)) :: sigma

        z=(x-xmu)/t
        
        where (dabs(z) > 36._long) 
           sigma=0._long
        elsewhere
           ff=1._long/(1._long + dexp(z))
           s1=ff*dlog(ff)
           s2=(1._long-ff)*dlog(1._long-ff)
           sigma=-(s1+s2)
        endwhere

      END FUNCTION sigma

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... LOG FERMI FUNCTION WITH PROPER LIMITS 
! ... THIS IS USED IN RESPONSE
      FUNCTION fermi_log(t,xmu,x)
        USE precision_definition
        REAL(long) :: t,xmu

        REAL(long), INTENT(IN), ALLOCATABLE, DIMENSION(:) :: x
        REAL(long), DIMENSION(size(x)) :: z
        REAL(long), DIMENSION(size(x)) :: fermi_log

        z=(x-xmu)/t
        
        where (z > 60._long) 
           fermi_log = 0._long
        elsewhere (z < -60._long)
           fermi_log = z
        elsewhere
           fermi_log=dlog( 1._long + dexp(-z) )
        endwhere

      END FUNCTION fermi_log
      
    END MODULE fermibose
