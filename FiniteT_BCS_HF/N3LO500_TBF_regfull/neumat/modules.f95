! ... DEFINES PRECISIONS
      MODULE precision_definition
        IMPLICIT NONE
        INTEGER, parameter :: ilong=selected_int_kind(10)
        INTEGER, parameter :: long=selected_real_kind(15,307)
        !INTEGER, parameter :: long=selected_real_kind(15,310)
        !	     INTEGER (ilong), parameter :: Nx=100
        LOGICAL :: rkspace ! -> TRUE IF R SPACE, FALSE IF K SPACE
      END MODULE precision_definition 

! ... PHYSICAL CONSTANTS TO BE USED
      MODULE physical_constants
        USE precision_definition
        REAL (long) :: pi,pi2,hbc,hbc2,hbc3,deg
        REAL (long) :: xmassn,xmassp,xmass,htmm
        REAL (long), DIMENSION(2) :: htm
        COMPLEX (long) :: zi
      END MODULE physical_constants
