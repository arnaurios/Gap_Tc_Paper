      MODULE PAIRING
        USE precision_definition

        REAL(long), ALLOCATABLE, DIMENSION(:) :: gap,gap0,auxgap,newgap,xkaux,spaux,dgap2,xid
        REAL(long), ALLOCATABLE, DIMENSION(:,:) :: vnngap
        LOGICAL :: coupled

        CONTAINS

          SUBROUTINE P_GAP
            USE precision_definition
            USE physical_constants
            USE NN_interaction
            USE mesh
            USE interpolation
            USE thermodynamical
            USE KSPACE_MODULE
            USE fermibose

            IMPLICIT NONE

            REAL (long), DIMENSION(1) :: ax,bx
            REAL(long) :: vpp,vpf,ffac,acc,test,sgn
            INTEGER(ilong) :: ij,ikmesh,jkmesh,iter,maxiter,nitersig
            INTEGER(ilong) :: ic,Nkaux,i1,i2,i3,i4
            CHARACTER(5) :: pwname
            LOGICAL :: writegap

            write(*,*) 'Pairing...'

            if(.not. ALLOCATED( xme ) ) ALLOCATE( xme(Nkmesh) )

! ... NUMERICAL PARAMETERS
            ffac=0.12345_long
            acc=1e-7
            maxiter=10000
            nitersig=2000

! ... All the partial waves (up to 9) are included
            do 999 ij=Jmin,Jmax

               do 998 ic=1,3
                  call AVFACNAME(iz1,iz2,ij,ic,pwname,coupled)

! ... IF UNPHYSICAL CHANNEL OR INPUT PWs IS NOT
                  if( ss == 2 .or. index(p_waves,trim(pwname)) .eq. 0 ) CYCLE

                  write(*,*) ij,pwname,ss,coupled

                  ! ... UNCOUPLED CHANNELS
                  if( .not.coupled) then
                     Nkaux=Nkmesh
                     i1=ic
                  elseif( coupled ) then
                     Nkaux=2*Nkmesh
                     i1=3
                     i2=5
                     i3=6
                     i4=4
                  endif
! ... ALLOCATE ALL ARRAYS
                  if(ALLOCATED(gap0)) DEALLOCATE( gap0 )
                  if(.not.ALLOCATED(gap0)) ALLOCATE( gap0(Nkmesh) )

                  if(ALLOCATED(gap)) DEALLOCATE(gap)
                  if(.not.ALLOCATED(gap)) ALLOCATE(gap(Nkaux))

                  if(ALLOCATED(xid)) DEALLOCATE(xid)
                  if(.not.ALLOCATED(xid)) ALLOCATE(xid(Nkaux))

                  if(ALLOCATED(auxgap)) DEALLOCATE( auxgap )
                  if(.not.ALLOCATED(auxgap)) ALLOCATE( auxgap(Nkaux) )

                  if(ALLOCATED(newgap)) DEALLOCATE( newgap )
                  if(.not.ALLOCATED(newgap)) ALLOCATE( newgap(Nkaux) )

                  if(ALLOCATED(xkaux)) DEALLOCATE( xkaux )
                  if(.not.ALLOCATED(xkaux)) ALLOCATE( xkaux(Nkaux) )

                  if(ALLOCATED(spaux)) DEALLOCATE( spaux )
                  if(.not.ALLOCATED(spaux)) ALLOCATE( spaux(Nkaux) )

                  if(ALLOCATED(dgap2)) DEALLOCATE(dgap2)
                  if(.not.ALLOCATED(dgap2)) ALLOCATE( dgap2(Nkaux) )

                  if(ALLOCATED(vnngap)) DEALLOCATE(vnngap)
                  if(.not.ALLOCATED(vnngap)) ALLOCATE( vnngap(Nkaux,Nkaux) )
! ... DONE ALLOCATING ALL ARRAYS

                  ! ... Diagonal channels only
                  gap=0_long

                  ! ... BUILD NkxNk matrix
                  FORALL(ikmesh=1:Nkmesh,jkmesh=1:Nkmesh)
                     vnngap(ikmesh,jkmesh) = vNN(ikmesh,jkmesh,ij,i1)*wkmesh(jkmesh)*xkmesh(jkmesh)**2
                  END FORALL
                  xkaux(1:Nkmesh)=xkmesh
                  spaux(1:Nkmesh)=spe

! ... FIRST GUESS OF THE GAP
! ... BASED ON FACT THAT GAP(k) ~ v(k,k_F), as EXPLAINED IN KUCHAREK ZPA 1989
                  ax=xkfermi
                  call LIN_INT2D(xkmesh,xkmesh,vNN(:,:,ij,i1),ax,ax,bx)
                  vpf=bx(1)
                  test=xkfermi/2_long*sqrt(abs(vpf))
                  do ikmesh=1,Nkmesh
                     call LIN_INT(xkmesh,vNN(ikmesh,:,ij,i1),ax,bx)
                     vpp=bx(1)
                     sgn=1_long
                     if(vpp<0) sgn=-1_long
!                     gap0(ikmesh)=vpp !-sgn*abs(vpp)
                     gap0(ikmesh)=-test*sgn*sqrt(abs(vpp/vpf))
                  enddo
! ... DANGEROUS TO TAKE gap=gap0, DYNAMICAL RESIZING OCCURS
                  gap(1:Nkmesh)=gap0

                  IF( coupled ) THEN
                     ! ... FILL IN THE MATRIX FOR THE COUPLED CHANNEL CASE
                     FORALL(ikmesh=1:Nkmesh,jkmesh=1:Nkmesh)
                        ! ... MINUS SIGN ACCORDING TO DEAN AND HJORTH-JENSEN
                        vnngap(ikmesh,Nkmesh+jkmesh) = -vNN(ikmesh,jkmesh,ij,i2)*wkmesh(jkmesh)*xkmesh(jkmesh)**2
                        vnngap(Nkmesh+ikmesh,jkmesh) = -vNN(ikmesh,jkmesh,ij,i3)*wkmesh(jkmesh)*xkmesh(jkmesh)**2
                        vnngap(Nkmesh+ikmesh,Nkmesh+jkmesh) = vNN(ikmesh,jkmesh,ij,i4)*wkmesh(jkmesh)*xkmesh(jkmesh)**2
                     END FORALL
                     xkaux(Nkmesh+1:Nkaux)=xkmesh
                     spaux(Nkmesh+1:Nkaux)=spe

                     call LIN_INT2D(xkmesh,xkmesh,vNN(:,:,ij,i4),ax,ax,bx)
                     vpf=bx(1)
                     test=xkfermi/2_long*sqrt(abs(vpf))
                     do ikmesh=1,Nkmesh
                        call LIN_INT(xkmesh,vNN(ikmesh,:,ij,i4),ax,bx)
                        vpp=bx(1)
                        sgn=1_long
                        if(vpp<0) sgn=-1_long
                        gap0(ikmesh)=-test*sgn*sqrt(abs(vpp/vpf))
                     enddo

                     gap(Nkmesh+1:Nkaux)=gap0
                  ENDIF

                  if( .not. coupled ) then
                     dgap=gap
                     dgap2=dgap**2
                  else
                     dgap=sqrt( gap(1:Nkmesh)**2 + gap(Nkmesh+1:Nkaux)**2 )
                     dgap2(1:Nkmesh)=dgap**2
                     dgap2(Nkmesh+1:Nkaux)=dgap**2
                  endif

!                  do ikmesh=1,Nkmesh
!                     write(188,'(5es18.6,2x,i1)') xkmesh(ikmesh),dgap(ikmesh),0._long,xid(ikmesh),0._long,ic
!                  enddo
!                  write(188,'(/,/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... ITERATION LOOP
                  ! spe=xkmesh**2*htm(1)/2_long
                  do iter=0,maxiter

                     !if( mod(iter,1) == 0 ) then
                        call ADJUST_MU( dgap )
                        !if( mod(iter,10) == 0 ) write(*,*) xmu
                        xme = MOMDIS(t,xmu,spe,dgap)
                     !endif

                     spaux(1:Nkmesh)=spe
                     if(coupled) spaux(Nkmesh+1:Nkaux)=spe

                     xid=sqrt( ( spaux-xmu )**2 + dgap2 )
                     auxgap=gap/2._long/xid*ttanh( t,xid )

                     ! ... CUTOFF AT VERY LARGE K
                     WHERE(xkaux > cutoff) auxgap=auxgap*exp(-xkaux/cutoff)

                     newgap=-MATMUL(vnngap,auxgap)

                     if( maxval( abs(newgap-gap) ) < acc) EXIT

                     gap=(1_long-ffac)*gap+ffac*newgap

                     if( coupled ) then
                        dgap=sqrt( gap(1:Nkmesh)**2 + gap(Nkmesh+1:Nkaux)**2 )
                        dgap2(1:Nkmesh)=dgap**2
                        dgap2(Nkmesh+1:Nkaux)=dgap**2
                     else
                        dgap=gap
                        dgap2=dgap**2
                     endif

                     if( any(isnan(dgap)) ) then
                        write(*,*) 'NAN IN DGAP'
                        stop
                     endif

                     if( any(dgap==dgap+1) ) then
                        write(*,*) 'INF IN DGAP'
                        stop
                     endif

                     if( maxval( abs(dgap) ) < acc) EXIT

!                     if(mod(iter,100) == 0) then
!                        write(188,'(a,2x,i6)') '#',iter
!                        write(188,*) '#',iter
!                        do ikmesh=1,Nkmesh
!                           write(188,'(5es18.6,2x,i1)') xkmesh(ikmesh),dgap(ikmesh),auxgap(ikmesh),xid(ikmesh),newgap(ikmesh),ic
!                        enddo
!                        write(188,'(/,/)')
!                     endif

                  enddo ! ITERATION LOOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  writegap=.true.

                  if( maxval( abs(newgap-gap) ) > acc .or. iter==maxiter) then
                     write(*,'(a,1x,2i2,1x,a)') 'bad gap',iz1,iz2,pwname
                     writegap=.false.
                     dgap=0_long
!                     RETURN
                  endif

                  if( maxval( abs(gap) ) < acc ) then
                     write(*,*) 'zero gap ',iz1,iz2,pwname
                     dgap=0_long
                  endif

                  if(writegap) call WRITE_GAP_OUT(ij,ic,iter)

                  if( ALLOCATED(gap) ) DEALLOCATE(gap)
                  if( ALLOCATED(xid) ) DEALLOCATE(xid)
                  if( ALLOCATED(gap0) ) DEALLOCATE(gap0)
                  if( ALLOCATED(auxgap) ) DEALLOCATE(auxgap)
                  if( ALLOCATED(newgap) ) DEALLOCATE(newgap)
                  if( ALLOCATED(xkaux) ) DEALLOCATE(xkaux)
                  if( ALLOCATED(dgap2) ) DEALLOCATE(dgap2)
                  if( ALLOCATED(vnngap) ) DEALLOCATE(vnngap)

998            enddo  ! BOXES LOOP

999         enddo ! PW LOOP

          END SUBROUTINE P_GAP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... SUBROUTINE TO WRITE OUT GAP INFO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SUBROUTINE WRITE_GAP_OUT(ij,ic,iter)
            USE precision_definition
            USE thermodynamical, ONLY : xkfermi,rho,t,dgap,xme
            USE mesh, ONLY : xkmesh,Nkmesh,xk,xkout,Nkout
            USE read_write
            USE interpolation

            IMPLICIT NONE

            INTEGER(ilong), INTENT(IN) :: ij,ic,iter
            REAL (long), ALLOCATABLE, DIMENSION(:) :: gkout

            INTEGER(ilong) :: ik
            REAL (long), DIMENSION(1) :: ax,bx
            REAL (long) :: gkf
            INTEGER :: un1,un2,un3

            un1=int(iunit1(ij,ic),4)
            un2=int(iunit2(ij,ic),4)
            un3=int(iunit3(ij,ic),4)

            ax=xkfermi
            call LIN_INT(xkmesh,dgap,ax,bx)
            gkf=bx(1)

            write(*,'(3f10.4,x,1es16.4,1x,i7)') &
                 rho,xkfermi,t,gkf,iter

96          format(10es16.6)

            write(un1,96) rho,t,xkfermi,gkf

            do ik=1,Nkmesh
               if( .not. coupled ) then
                  write(un2,96) &
                       rho,t,xkfermi,xk(ik),dgap(ik),xid(ik),0d0,xme(ik),gap(ik)
               else
                  write(un2,96) &
                       rho,t,xkfermi,xk(ik),dgap(ik),xid(ik),0d0,xme(ik),gap(ik),gap(ik+Nkmesh)
               endif
            enddo

101         format(/,/)
            write(un2,101)

            ALLOCATE( gkout(Nkout) )

            call LIN_INT(xkmesh,dgap,xkout,gkout)
            do ik=1,Nkout
               write(un3,96) &
                          rho,t,xkfermi,xkout(ik),gkout(ik)
            enddo
            write(un3,*)

          END SUBROUTINE WRITE_GAP_OUT

     END MODULE PAIRING
