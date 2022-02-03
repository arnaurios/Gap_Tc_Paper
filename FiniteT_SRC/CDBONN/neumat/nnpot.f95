! ... NN interaction passed around
      MODULE NN_INTERACTION
        USE precision_definition
        INTEGER (ilong) :: Jmin,Jmax,Nch,iz1,iz2,Nmsh
        REAL (long) :: cmaxk
        REAL (long), ALLOCATABLE, DIMENSION(:,:,:,:) :: vNN
        CHARACTER(len=100) :: p_waves
        
        CONTAINS
        
! ... cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... Reads the NN interaction as needed in the program     c
! ... cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE NNPOT

        USE precision_definition
        USE physical_constants
        USE mesh

        IMPLICIT NONE

        INTEGER :: jc
        INTEGER(ilong) :: ik,jk,inn

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

        write(*,*) 
        write(*,*) 'Loading the NN potential...'

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
      if (iz1.eq.-1 .and. iz2.eq.-1) then 
         inn=3
      elseif(iz1.eq.1 .and. iz2.eq.1) then 
         inn=1
      else
         inn=2
      endif

! ... NN potential 
      do jc=Jmin,Jmax

         do ik=1,Nmsh

            xmev=xk(ik)
            
            do jk=1,Nmsh

               ymev=xk(jk)
               
                  call cdbonnNN

! ... UNCOUPLED STATES
! ... Singlet
               vNN(ik,jk,jc,1)=v(1)

! ... Uncoupled triplet
               vNN(ik,jk,jc,2)=v(2)

! ... COUPLED STATES
! ...  Coupled triplet V--
               vNN(ik,jk,jc,3)=v(4)
! ...  Coupled triplet V++
               vNN(ik,jk,jc,4)=v(3)
! ... 3 four states are the diagonal waves

! ...   Coupled triplet V-+
               vNN(ik,jk,jc,5)=v(6)
! ...   Coupled triplet V+-
               vNN(ik,jk,jc,6)=v(5)

! ... Empty sub boxes
               vNN(ik,jk,jc,7)=0.d0
               vNN(ik,jk,jc,8)=0.d0

!               if(jc == 0) then
!                  vNN(ik,jk,jc,2)=v(3)
!                  vNN(ik,jk,jc,4)=0d0
!               endif

            enddo               ! End of loop over jk      
         enddo                  ! End of loop over ik
      enddo                     ! End of loop over jc

      write(*,*) vNN(1,1,0,1)
      
      write(*,*) 'Finished loading the NN potential!'

      END SUBROUTINE NNPOT


    END MODULE NN_INTERACTION
            
