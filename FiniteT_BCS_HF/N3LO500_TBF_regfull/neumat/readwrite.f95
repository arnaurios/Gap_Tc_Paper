   MODULE read_write

     USE precision_definition

     INTEGER (ilong), ALLOCATABLE, DIMENSION(:,:) :: iunit1,iunit2,iunit3

     CONTAINS
       
! ... cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... WRITES OUT AN NN INTERACTION IN AND OUT OF THE DIAGONAL
! ... cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE WRITEVNN
        USE physical_constants, ONLY : hbc,xmass,pi
        USE NN_interaction
        USE mesh
        USE interpolation
        USE thermodynamical, ONLY : xkfermi
        IMPLICIT NONE

        INTEGER(ilong) :: jj,iju,ija,ik,jk,ll
        CHARACTER(1) :: Ja
        REAL (long), DIMENSION(1) :: ax !,bx
        REAL (long), DIMENSION(8) :: bx
        REAL (long) :: fact
        
129     format('#',3x,'k[fm]',10x,'V_NN(ich=1,8) [fm]')
139     format('#',3x,'k1[Mev]',3x,'k2[Mev]',10x,'V_NN(ich=1,8) [fm]')
        
!        fact=1d0
        fact=4_long*pi/xmass/(2_long*pi**2)/hbc !2/m/pi       

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
           write(ija,139)
           
           do ik=1,Nmsh
              write(iju,'(10e15.5)') xk(ik),(vNN(ik,ik,jj,ll)/fact,ll=1,8)
              do jk=1,Nmsh
!                 write(ija,'(10e15.5)') (xkn(ik)/hbc)**2,(xkn(jk)/hbc)**2 &
                 write(ija,'(12e15.5)') xk(ik)/hbc,xk(jk)/hbc &
!                      ,(xk(ik)/hbc)**2,(xk(jk)/hbc)**2 &
                      ,(vNN(ik,jk,jj,ll)/fact,ll=1,6)
              enddo
              write(ija,*)
           enddo
           write(ija,'(/,/)')
           write(iju,'(/,/)')
        enddo
        

        ! ... WRITE OUT PAIRING INTERACTION
229     format('#',3x,'k[fm]',10x,'V_NN(k=kf,k,ich=1,8) [fm] kf=',f10.4)
        
        iju=229
        do jj=Jmin,Jmax
           iju=iju+1

           write(Ja,'(SS,i1)') jj
           open(unit=iju,file='1dkf_vNN_J'//trim(Ja)//'.dat',status='unknown')
           write(*,*) '1dkf_vNN_J'//trim(Ja)//'.dat',xkfermi

           write(iju,229) xkfermi

           ax=xkfermi*hbc                 
           do ik=1,Nmsh
              do ll=1,8
                 call LIN_INT(xk,vNN(ik,:,jj,ll),ax,bx(ll) )
              enddo
              write(iju,'(10e15.5)') xk(ik)/hbc,(bx(ll)/fact,ll=1,8)
           enddo
           write(iju,'(/,/)')
        enddo

      END SUBROUTINE WRITEVNN

   END MODULE read_write
