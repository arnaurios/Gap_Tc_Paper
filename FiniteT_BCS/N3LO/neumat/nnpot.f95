! ... NN interaction passed around
      MODULE NN_INTERACTION
        USE precision_definition
        INTEGER (ilong) :: Nch,iz1,iz2,Nmsh
        INTEGER (ilong) :: Jminhf,Jmaxhf,Jmin,Jmax
        REAL (long) :: cmaxk
        REAL (long), ALLOCATABLE, DIMENSION(:,:,:,:) :: vNN,vHF
        CHARACTER(len=100) :: p_waves

        CONTAINS

! ... cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... Reads the NN interaction as needed in the program     c
! ... cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE NNPOT(jmn,jmx,xk,vv)

        USE precision_definition
        USE physical_constants
!        USE mesh

        IMPLICIT NONE
        INTEGER(ilong), INTENT(IN) :: jmn,jmx
        REAL(long), ALLOCATABLE, DIMENSION(:),INTENT(IN) :: xk
        REAL(long), ALLOCATABLE, DIMENSION(:,:,:,:),INTENT(OUT) :: vv
        INTEGER :: jc
        INTEGER(ilong) :: ik,jk
        INTEGER :: inn

! ... POTENTIAL DEPENDENT SECTION
! ... N3LO POTENTIAL COMMONS
        DOUBLE PRECISION  :: v,xmev,ymev
        LOGICAL :: heform,sing,trip,coup,endep
        CHARACTER*4 :: label
        INTEGER :: kda,kwrite,kread,kpunch
        common /cpot/v(6),xmev,ymev
        common /cstate/jc,heform,sing,trip,coup,endep,label
        common /crdwrtNN/ kread,kwrite,kpunch,kda(9)
        common /cnn/ inn
! ... END OF POTENTIAL DEPENDENT SECTION

        Nmsh=size(xk)
        ALLOCATE( vv(Nmsh,Nmsh,jmn:jmx,Nch) )

! ... DATA OF THE CDBONN POTENTIAL
        heform=.false.
        sing=.true.
        trip=.true.
        coup=.true.
        kread=5
        kwrite=6

! ... CCCCCCCCCCCCCC  HERE WE CALL THE NN POTENTIAL CCCCCCCCCCCCCCCCCC

!     inn=1  means pp potential,
!     inn=2  means np potential, and
        !     inn=3  means nn potential.
        if( iz1*iz2 == 1 .or. iz1*iz2 ==-1) then
           if (iz1.eq.-1 .and. iz2.eq.-1) then
              inn=3
           elseif(iz1.eq.1 .and. iz2.eq.1) then
              inn=1
           else
              inn=2
           endif
        else
           write(*,*) 'Wrong isospins!',iz1,iz2
           stop
        endif

! ... NN potential
      do jc=jmn,jmx

         do ik=1,Nmsh

            xmev=xk(ik)

            do jk=1,Nmsh

               ymev=xk(jk)

               call N3LO

! ... UNCOUPLED STATES
! ... Singlet
               vv(ik,jk,jc,1)=v(1)

! ... Uncoupled triplet
               vv(ik,jk,jc,2)=v(2)

! ... COUPLED STATES
! ...  Coupled triplet V--
               vv(ik,jk,jc,3)=v(4)
! ...  Coupled triplet V++
               vv(ik,jk,jc,4)=v(3)
! ... 3 four states are the diagonal waves

! ...   Coupled triplet V-+
               vv(ik,jk,jc,5)=v(6)
! ...   Coupled triplet V+-
               vv(ik,jk,jc,6)=v(5)

! ... Empty sub boxes
               vv(ik,jk,jc,7)=0.d0
               vv(ik,jk,jc,8)=0.d0

               if(jc == 0) then
                  vv(ik,jk,jc,2)=v(3)
                  vv(ik,jk,jc,4)=0d0
               endif

            enddo               ! End of loop over jk
         enddo                  ! End of loop over ik
      enddo                     ! End of loop over jc
!      write(*,*) 'Finished loading the NN potential!'

      END SUBROUTINE NNPOT

      INTEGER(ilong) FUNCTION ISOSPIN(J,ICH)
        USE precision_definition
        IMPLICIT NONE
        INTEGER(ilong) :: J,ICH

        if(J == 0) then
           isospin=1
        else
           ! ... If the channel is the singlet or the coupled triplet, the isospin is given by J modulus 2
           if(ich==1 .or. ich==3 .or. ich==4) then
              isospin=mod(j+1,2)
              ! ... Otherwise it differs by a unit
           elseif(ich==2) then
              isospin=mod(j,2)
           endif
        endif

      END FUNCTION ISOSPIN


    END MODULE NN_INTERACTION
