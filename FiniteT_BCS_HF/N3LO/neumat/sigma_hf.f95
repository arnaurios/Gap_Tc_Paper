      MODULE selfenergy
        USE precision_definition
        REAL (long), ALLOCATABLE, DIMENSION(:) :: sigma_hf,dgap,xme

      CONTAINS

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computes the generalized HF self-energy               c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE HARTREE_FOCK(v,shf)

        USE precision_definition
        USE physical_constants
        USE NN_interaction
        USE mesh
        USE mesh_generator
        USE interpolation
        USE thermodynamical
        USE fermibose

        IMPLICIT NONE
        REAL(long), ALLOCATABLE, DIMENSION(:,:,:,:), INTENT(IN) :: v
        REAL(long), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: shf
        REAL(long) :: vchan,vkint
        REAL(long), ALLOCATABLE, DIMENSION(:) :: xq,vpp
        REAL(long), ALLOCATABLE, DIMENSION(:) :: vaux
        REAL(long), ALLOCATABLE, DIMENSION(:) :: xintegrand,speint,xmint,gapint
        INTEGER(ilong) :: ikmesh,ikint,ij,ich,ik,isos

!        open(unit=99,file='hfself.dat',status='unknown')

        if(.NOT. ALLOCATED( shf ) ) ALLOCATE( shf(Nkmesh) )
        if(.not. ALLOCATED( xme ) ) ALLOCATE( xme(Nkmesh) )

! ... ANGULAR INTEGRAL
        ALLOCATE( xq(Nth),vpp(Nth) )


! ... INTERNAL MOMENTUM INTEGRAL
        ALLOCATE( vaux(Nkmeshhf) )
        ALLOCATE( xintegrand(Nkmeshhf) )
        ALLOCATE( speint(Nkmeshhf),gapint(Nkmeshhf),xmint(Nkmeshhf) )

        xme = MOMDIS(t,xmu,spe,dgap)

        call SPL_INT( xk,spe,xkhf,speint,'F')
        call SPL_INT( xk,dgap,xkhf,gapint,'F')
        xmint = MOMDIS(t,xmu,speint,gapint)

!199     format(/,/,'#',4x,'k [MeV]',4x,'S_HF [MeV]',4x,'n(k)',/,'# xkf=',f10.6,' MeV')
!        write(99,199) xkfermi*hbc

! ... EXTERNAL LOOP
        do 701 ikmesh=1,Nkmesh
           do 702 ikint=1,Nkmeshhf

              xq = sqrt( xk(ikmesh)**2 + xkhf(ikint)**2 &
                       + 2d0*xk(ikmesh)*xkhf(ikint)*theta )/2._long

              vkint=0._long
              do 703 ij=jminhf,jmaxhf
                 do 704 ich=1,4
                    isos=isospin(ij,ich)

                    if( isos == 1 ) then

                       FORALL( ik=1:Nkmeshhf )
                          vaux(ik)= v(ik,ik,ij,ich)
                       END FORALL

                       call SPL_INT( xkhf,vaux,xq,vpp,'F' )

                       WHERE( xq > xkhf(Nkmeshhf) ) vpp=0._long

                       vchan = real(2*ij+1,long)*(2_long)*sum(vpp*wtheta)/(8d0*pi)

                       vkint = vkint + vchan
                    endif
704              enddo ! channel loop

703           enddo ! Partial wave loop
              xintegrand(ikint) = xkhf(ikint)**2*vkint*xmint(ikint)
!              write(*,*) xkhf(ikint),xintegrand(ikint)
702        enddo ! Momentum Sum loop
           shf(ikmesh) = 2._long*pi*sum(xintegrand*wkhf)/hbc2

!           write(99,'(3es16.6)') xk(ikmesh),shf(ikmesh),xme(ikmesh)

701     enddo ! External momentum loop

        DEALLOCATE( xq,vpp,xintegrand,vaux,speint,gapint,xmint )

      END SUBROUTINE HARTREE_FOCK


    END MODULE selfenergy
