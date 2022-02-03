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
        REAL (long) :: pi,pi2,hbc,hbc2,hbc3,deg,fact
        REAL (long) :: xmassn,xmassp,xmass,htmm
        REAL (long), DIMENSION(2) :: htm
        COMPLEX (long) :: zi
      END MODULE physical_constants

! ... MESH IN X AND K SPACE
      MODULE mesh
        USE precision_definition
        REAL (long) :: delta,cutoff
        INTEGER (ilong) :: Nkmesh
        REAL (long), ALLOCATABLE, DIMENSION(:) :: xk,wk,xkn,xkw
        REAL (long), ALLOCATABLE, DIMENSION(:) :: xkmesh,wkmesh
      END MODULE mesh

! ... TD PROPERTIES
      MODULE thermodynamics
        USE precision_definition
        INTEGER (ilong) :: nrho,nxkf,nasy,ntemp
        REAL (long) :: xxkf,xkfin,dxkf
        REAL (long) :: rho,rhoin,drho
        REAL (long) :: asy,asyin,asyfi,dasy
        REAL (long) :: xkfermi,ef
        REAL (long) :: rhon,rhop,xkfn,xkfp,efn,efp,xmun,xmup
        REAL (long) :: temp
!        REAL (long),ALLOCATABLE,DIMENSION(:) :: xmu,xkf,ef

      END MODULE thermodynamics
