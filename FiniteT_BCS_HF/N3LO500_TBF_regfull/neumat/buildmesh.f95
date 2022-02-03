! ... MESH IN X AND K SPACE
      MODULE mesh
        USE precision_definition
        REAL (long) :: delta,cutoff
        INTEGER (ilong) :: Nkmesh,Nkmeshhf,Nkout,Nth
        REAL (long), ALLOCATABLE, DIMENSION(:) :: xk,wk,xkn,xkw
        REAL (long), ALLOCATABLE, DIMENSION(:) :: xkmesh,wkmesh
        REAL (long), ALLOCATABLE, DIMENSION(:) :: xkhf,wkhf
        REAL (long), ALLOCATABLE, DIMENSION(:) :: xkout,wkout
        REAL (long), ALLOCATABLE, DIMENSION(:) :: theta,wtheta
        INTEGER(ilong), ALLOCATABLE, DIMENSION(:) :: N
        
      CONTAINS

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... DEFINES MESHES 
! ... IN THE FIRST ITERATION, MESHES ARE DEFINED AND DATA
! ... IS INTERPOLATED FROM WHAT HAS BEEN READ IN READDAT.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE BUILDMESHES(xkfermi)

        USE precision_definition
        USE physical_constants
        USE mesh_generator
        USE interpolation
        USE fermibose

        IMPLICIT NONE
        REAL(long), INTENT(IN) :: xkfermi

        REAL(long), ALLOCATABLE, DIMENSION(:) :: xmsh
        REAL(long) :: xi,xf,cmaxk
        INTEGER(ilong) :: Ni,Nf
        INTEGER(ilong) :: ireg,Nreg
        
        !     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! ... ALLOCATE ARRAYS
        Nreg=6 ! +1 on tangential
        if( ALLOCATED(N) ) DEALLOCATE( N )
        ALLOCATE( N(Nreg+1),xmsh(Nreg+1) )
        
! ... MAGIC SET OF NUMBERS THAT MAKE IT CONVERGE!        
        N(1)=48
        xmsh(1)=0.6_long
        N(2)=100
        xmsh(2)=0.99_long
        N(3)=100
        xmsh(3)=1.01_long
        N(4)=48
        xmsh(4)=1.4_long
        N(5)=24
        xmsh(5)=2.5_long
        N(6)=24
        xmsh(6)=5_long
        N(7)=24
        cmaxk=1e6

        Nkmesh=sum(N)
        if(ALLOCATED(xkmesh)) DEALLOCATE( xkmesh,wkmesh )
        ALLOCATE( xkmesh(Nkmesh),wkmesh(Nkmesh) )

        if(ALLOCATED(xk)) DEALLOCATE( xk,wk )
        ALLOCATE( xk(Nkmesh),wk(Nkmesh) )
        
        ! ... cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! ... LOOP OVER PARTICLE TYPE             
        ! ... MESH OF EXTERNAL MOMENTA
        Nf=0
        do ireg=1,Nreg
           xi=real(0,long)
           if(ireg>1) xi=xmsh(ireg-1)*xkfermi
           xf=xmsh(ireg)*xkfermi
           Ni=Nf+1
           Nf=Ni-1+N(ireg)
           call GAUSS(xi,xf,N(ireg),xkmesh(Ni:Nf),wkmesh(Ni:Nf))
        enddo
        
        xi=xf
        xf=cmaxk/hbc
        Ni=Nf+1
        Nf=Ni-1+N(Nreg+1)
        call TANGENTIAL( xi,xf,N(Nreg+1),xkmesh(Ni:Nf),wkmesh(Ni:Nf),'U')
        !call LOGARITHMIC( xi,xf,N(5),xkmesh(Ni:Nf),wkmesh(Ni:Nf),'U')
        
        ! ... CALL NNFORCE IN THE SAME MESH
        xk=xkmesh*hbc

        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! ... MESH TO STORE VNN FOR HF CALCULATION
        ! ... MESH OF MOMENTA TO STORE VHF
        Nreg=2
        N(1)=42
        xmsh(1)=1._long
        N(2)=42
        xmsh(2)=2._long
        N(3)=24
        Nkmeshhf=sum(N(1:3))

        if( ALLOCATED(xkhf)) DEALLOCATE( xkhf,wkhf )
        ALLOCATE( xkhf(Nkmeshhf),wkhf(Nkmeshhf) )
        
        Nf=0
        do ireg=1,Nreg
           xi=real(0,long)
           if(ireg>1) xi=xmsh(ireg-1)*xkfermi*hbc
           xf=xmsh(ireg)*xkfermi*hbc
           Ni=Nf+1
           Nf=Ni-1+N(ireg)
           call GAUSS(xi,xf,N(ireg),xkhf(Ni:Nf),wkhf(Ni:Nf))
        enddo
        xi=xf
        cmaxk=100_long*hbc
        xf=cmaxk
        Ni=Nf+1
        Nf=Ni-1+N(Nreg+1)
        call TANGENTIAL( xi,xf,N(Nreg+1),xkhf(Ni:Nf),wkhf(Ni:Nf),'U')

        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! ... FOR 2D DATA, INTERPOLATE TO kf x kout mesh
        Nkout=100
        if( ALLOCATED( xkout ) ) DEALLOCATE( xkout,wkout )
        ALLOCATE( xkout(Nkout),wkout(Nkout) )
        xi=real(0d0,long)
        xf=real(10d0,long)
        CALL LINEAR(xi,xf,Nkout,xkout,wkout)
        
        
        !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! ... GAUSSIAN MESH-POINTS FOR ANGULAR INTEGRATION
        Nth=32
        if(.not. ALLOCATED(theta) ) ALLOCATE( theta(Nth),wtheta(Nth) )
        call GAUSS(-1._long,1._long,Nth,theta,wtheta)

        DEALLOCATE( xmsh )
        
      END SUBROUTINE BUILDMESHES
      
    END MODULE MESH
