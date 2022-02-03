   MODULE read_write

     USE precision_definition

     INTEGER (ilong), ALLOCATABLE, DIMENSION(:,:) :: iunit1,iunit2,iunit3

     CONTAINS

! ... cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... WRITES OUT AN NN INTERACTION IN AND OUT OF THE DIAGONAL
! ... cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE WRITEVNN
        USE physical_constants, ONLY : hbc,hbc3,pi,xmass
        USE NN_interaction
        USE mesh
        IMPLICIT NONE
        REAL(long) :: fact
        INTEGER(ilong) :: jj,iju,ija,ik,jk,ll
        CHARACTER(1) :: Ja

129     format('#',3x,'k[fm]',10x,'V_NN(ich=1,8) [fm]')

        fact=1d0
        fact=hbc3 * 2d0/pi/xmass/hbc

! ... WRITE OUT NN INTERACTION
        iju=129
        ija=139
        do jj=Jmin,Jmax
           iju=iju+1
           ija=ija+1

           write(Ja,'(SS,i1)') jj
           open(unit=iju,file='diag_vNN_J'//trim(Ja)//'.dat',status='unknown')
           open(unit=ija,file='2d_vNN_J'//trim(Ja)//'.dat',status='unknown')
           write(*,*) '2d_vNN_J'//trim(Ja)//'.dat'

           write(iju,129)
           do ik=1,Nmsh
              write(iju,'(10es15.5)') xk(ik)/hbc,(vNN(ik,ik,jj,ll)/fact,ll=1,8)
              do jk=1,Nmsh
!                 write(ija,'(10e15.5)') (xkn(ik)/hbc)**2,(xkn(jk)/hbc)**2 &
                 write(ija,'(12es15.5)') xk(ik)/hbc,xk(jk)/hbc &
                      ,(xk(ik)/hbc)**2,(xk(jk)/hbc)**2 &
                      ,(vNN(ik,jk,jj,ll)/fact,ll=1,6)
              enddo
              write(ija,*)
           enddo
           write(ija,'(/,/)')
           write(iju,'(/,/)')
        enddo

      END SUBROUTINE WRITEVNN

   END MODULE read_write
