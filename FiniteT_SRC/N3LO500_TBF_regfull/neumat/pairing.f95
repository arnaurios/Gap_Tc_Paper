      MODULE PAIRING
        USE precision_definition
        REAL(long), ALLOCATABLE, DIMENSION(:) :: gap,gap0,auxgap,newgap,xkaux,dgap,dgap2,xid
        REAL(long), ALLOCATABLE, DIMENSION(:,:) :: vnngap

        REAL (long), ALLOCATABLE, DIMENSION(:) :: xkout,gkout

        INTEGER (ilong) :: Nkout
        LOGICAL :: xi_inp

        CONTAINS

          SUBROUTINE P_GAP
            USE precision_definition
            USE physical_constants
            USE NN_interaction
            USE mesh
            USE interpolation
            USE thermodynamics
            USE KSPACE_MODULE

            IMPLICIT NONE

            REAL (long), DIMENSION(1) :: ax,bx
            REAL(long), ALLOCATABLE, DIMENSION(:) :: xidr,xidd

            REAL(long) :: vpp,vpf,ffac,acc,test,sgn
            INTEGER(ilong) :: ij,ikmesh,jkmesh,iter,maxiter,nitersig
            INTEGER(ilong) :: ic,Nkaux,i1,i2,i3,i4
            CHARACTER(5) :: pwname
            LOGICAL :: coupled,writegap

            write(*,*) 'Pairing...'

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
                  if(ALLOCATED(dgap)) DEALLOCATE(dgap)
                  if(.not.ALLOCATED(dgap)) ALLOCATE( dgap(Nkmesh) )

                  if(ALLOCATED(gap0)) DEALLOCATE( gap0 )
                  if(.not.ALLOCATED(gap0)) ALLOCATE( gap0(Nkmesh) )

                  if(ALLOCATED(gap)) DEALLOCATE(gap)
                  if(.not.ALLOCATED(gap)) ALLOCATE(gap(Nkaux))

                  if(ALLOCATED(xidr)) DEALLOCATE(xidr)
                  if(.not.ALLOCATED(xidr)) ALLOCATE(xidr(Nkaux))

                  if(ALLOCATED(xidd)) DEALLOCATE(xidd)
                  if(.not.ALLOCATED(xidd)) ALLOCATE(xidd(Nkaux))

                  if(ALLOCATED(auxgap)) DEALLOCATE( auxgap )
                  if(.not.ALLOCATED(auxgap)) ALLOCATE( auxgap(Nkaux) )

                  if(ALLOCATED(newgap)) DEALLOCATE( newgap )
                  if(.not.ALLOCATED(newgap)) ALLOCATE( newgap(Nkaux) )

                  if(ALLOCATED(xkaux)) DEALLOCATE( xkaux )
                  if(.not.ALLOCATED(xkaux)) ALLOCATE( xkaux(Nkaux) )

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
                     xidr=xid
                  else
                     dgap=sqrt( gap(1:Nkmesh)**2 + gap(Nkmesh+1:Nkaux)**2 )
                     dgap2(1:Nkmesh)=dgap**2
                     dgap2(Nkmesh+1:Nkaux)=dgap**2
                     xidr(1:Nkmesh)=xid
                     xidr(Nkmesh+1:Nkaux)=xid
                  endif

!                  do ikmesh=1,Nkmesh
!                     write(188,'(5es18.6,2x,i1)') xkmesh(ikmesh),dgap(ikmesh),0_long,xid(ikmesh),0_long,ic
!                  enddo
!                  write(188,'(/,/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... ITERATION LOOP
                  do iter=0,maxiter

                     if(.not. xi_inp) xidd=2_long*sqrt( ( (xkaux(:)**2-xkfermi**2)*htm(1)/2_long )**2 + dgap2 )
                     if( xi_inp ) xidd=2_long*sqrt( ( xidr )**2 + dgap2 )
                     auxgap=gap/xidd

                     ! ... CUTOFF AT VERY LARGE K
                     WHERE(xkaux > cutoff) auxgap=auxgap*exp(-xkaux/cutoff)

                     newgap=-MATMUL(vnngap,auxgap)

                     if( maxval( abs(newgap-gap) ) < acc) EXIT

!                     ffac=1_long-1_long/(1_long+exp((real(iter-nitersig/2))/(nitersig/10)))
!                     write(*,*) iter,ffac
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

 !                    if(mod(iter,100) == 0) then
                        !                        write(188,'(a,2x,i6)') '#',iter
 !                       write(188,*) '#',iter
 !                       do ikmesh=1,Nkmesh
 !                          write(188,'(5es18.6,2x,i1)') xkmesh(ikmesh),dgap(ikmesh),auxgap(ikmesh),xid(ikmesh),newgap(ikmesh),ic
 !                       enddo
 !                       write(188,'(/,/)')
 !                    endif

                  enddo ! ITERATION LOOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  writegap=.true.

                  if( maxval( abs(newgap-gap) ) > acc .or. iter==maxiter) then
                     write(*,'(a,1x,2i2,1x,a)') 'bad gap',iz1,iz2,pwname
                     write(*,*)  maxval( abs(newgap-gap) ), iter
                     writegap=.false.
                     dgap=0_long
                  endif

                  if( maxval( abs(gap) ) < acc ) then
                     write(*,*) 'zero gap ',iz1,iz2,pwname
                     dgap=0_long
                  endif

                  if(writegap) call WRITE_GAP_OUT(ij,ic,iter)

                  if( ALLOCATED(dgap) ) DEALLOCATE(dgap)
                  if( ALLOCATED(gap) ) DEALLOCATE(gap)
                  if( ALLOCATED(xidr) ) DEALLOCATE(xidr)
                  if( ALLOCATED(xidd) ) DEALLOCATE(xidd)
                  if( ALLOCATED(gap0) ) DEALLOCATE(gap0)
                  if( ALLOCATED(auxgap) ) DEALLOCATE(auxgap)
                  if( ALLOCATED(newgap) ) DEALLOCATE(newgap)
                  if( ALLOCATED(xkaux) ) DEALLOCATE(xkaux)
                  if( ALLOCATED(dgap2) ) DEALLOCATE(dgap2)
                  if( ALLOCATED(vnngap) ) DEALLOCATE(vnngap)

998            enddo  ! BOXES LOOP
999      enddo ! PW LOOP

          END SUBROUTINE P_GAP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... SUBROUTINE TO WRITE OUT GAP INFO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SUBROUTINE WRITE_GAP_OUT(ij,ic,iter)
            USE precision_definition
            USE thermodynamics, ONLY : xkfermi,rho,temp
            USE mesh, ONLY : xkmesh,Nkmesh
            USE read_write
            USE interpolation

            IMPLICIT NONE

            INTEGER(ilong), INTENT(IN) :: ij,ic,iter

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
                 rho,xkfermi,temp,gkf,iter

96          format(10es16.6)

            !                  write(un1,'(2es16.6,1x,i5)') &
            !                       rho,ax,bx,iter

            write(un1,96) xkfermi,temp,gkf

            do ik=1,Nkmesh
               write(un2,96) &
                    rho,temp,xkfermi,xkmesh(ik),dgap(ik),xid(ik),gap0(ik)
            enddo

101         format(/,/)
            write(un2,101)
!            write(188,101)

            call LIN_INT(xkmesh,dgap,xkout,gkout)
            do ik=1,Nkout
               write(un3,96) &
                          rho,temp,xkfermi,xkout(ik),gkout(ik)
            enddo
            write(un3,*)

          END SUBROUTINE WRITE_GAP_OUT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... READS FROM g2a_rho FILES THE ENERGY DENOMINATORS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SUBROUTINE READXI
            USE physical_constants, ONLY : xmass,hbc
            USE NN_interaction
            USE mesh
            USE thermodynamics
            USE interpolation

            IMPLICIT NONE
            CHARACTER(5) :: rhoc
            CHARACTER(4) :: tchar
            CHARACTER(100) :: filename,filename2
            INTEGER(ilong) :: ik,Nkread
            REAL(long), DIMENSION(:), ALLOCATABLE :: xkr,xir
            REAL(long), DIMENSION(:), ALLOCATABLE :: xkrr,xirr
            REAL :: x0,y0,x1,y1,ps0,ps1,ps,effm
            LOGICAL :: file_exists
            REAL(long) :: sr,edkr,xiqpr,xkfmv,rel_width_gap
            INTEGER(ilong) :: ikk,Nkk

            write(rhoc,'(SS,f4.2)') rho

            if(rho.lt.0.02d0) then
               write(rhoc,'(SS,F4.3)') rho
            endif

            if(rho .ge. 0.02d0 .and. mod(rho*1000d0,10d0) .ne. 0d0) then
               write(rhoc,'(SS,F4.3)') rho
            endif

            if( temp > 0.d0 ) then
              write(tchar,'(SS,f3.1)') temp
              if( mod(temp*100d0,10d0) .ne. 0d0) then
                write(tchar,'(SS,f4.2)') temp
              endif

              filename='g2a_rho'//trim(rhoc)//'_t'//trim(tchar)//'/wim_t0.dat'
              filename2='g2a_rho'//trim(rhoc)//'_t'//trim(tchar)//'/wim_mesh.dat'

            else
              filename='g2a_rho'//trim(rhoc)//'/wim_t0.dat'
              filename2='g2a_rho'//trim(rhoc)//'/wim_mesh.dat'
            endif

            write(*,*) filename
            inquire(FILE=trim(filename), EXIST=file_exists)
            if( .not.file_exists ) then
               write(*,*) 'INPUT FILE NOT AVAILABLE AT THIS RHO AND T'
               write(*,*) filename
               stop
            else
               open(unit=200,file=trim(filename),status='unknown')
               open(unit=100,file=trim(filename2),status='unknown')

               read(200,'(/,/,/)')

96             format(6e16.6)

               Nkread=100
               ALLOCATE( xkrr(Nkread),xirr(Nkread) )
               write(*,*) "first read",Nkread
               ik=0
               xkfmv=xkfermi*hbc
               rel_width_gap=0.0d0 !2d0 !2d0
               do ikk=1,Nkread
                   read(200,96) xkrr(ikk),xirr(ikk),sr,edkr,xiqpr
                   if( abs(xkrr(ikk) - xkfmv)/xkfmv > rel_width_gap ) ik=ik+1
               enddo
               Nkk=ik

               if( temp == 0.d0 ) then
                 ALLOCATE( xkr(Nkread),xir(Nkread) )
                 xkr=xkrr
                 xir=xirr
              else
                 Nkk=ik
                 write(*,*) "second read",Nkk
                 ALLOCATE( xkr(Nkk),xir(Nkk) )

                 rewind(200)
                 read(200,'(/,/,/)')
                 ik=0
                 do ikk=1,Nkread
                   if( abs(xkrr(ikk) - xkfmv)  /xkfmv > rel_width_gap ) then
                     ik=ik+1
                     read(200,96) xkr(ik),xir(ik),sr,edkr,xiqpr
                   else
                     read(200,96) xkrr(ikk),xirr(ikk),sr,edkr,xiqpr
                        !write(*,*) "ala",xkrr(ikk),xkrr(ikk)/xkfmv
                   endif
                 enddo
               endif

               DEALLOCATE( xkrr,xirr )
               close(200)

               !            call LIN_INT( xkr,xir,xk,xid )
               call SPL_INT( xkr,xir,xk,xid )

               ! READ DATA HAS A REAL GAP BETWEEK KF-delta AND KF+DELTA
               write(*,*) minval( abs(xkr - xkfermi*hbc) )/xkfermi/hbc !,(xkr(41)-xkfermi*hbc)/xkfermi*hbc

               ! ... EXTRAPOLATE TO HIGH MOMENTA NEEDED FOR BCS EQUATION
               if( maxval(xk) > maxval(xkr) ) then
                  write(*,*) 'Xi mesh small, extrapolate',maxval(xk),maxval(xkr)
                  x0=xkr(Nkk-1)
                  y0=xir(Nkk-1)
                  x1=xkr(Nkk)
                  y1=xir(Nkk)
                  !           write(*,*) x0,y0,x1,y1
               !           WHERE( xk > maxval(xkr) ) xid=y0 + (y1-y0)/(x1-x0)*(xk - x0)
                  ps0=y0/x0**2
                  ps1=y1/x1**2
                  ps=(ps0+ps1)/2d0
                  !           write(*,*) ps0,ps1,ps
                  effm=1d0/2d0/ps/xmass
                  write(*,*) 'ef',effm
                  WHERE( xk > maxval(xkr) ) xid=ps*xk**2
               endif

               do ik=1,Nkmesh
                  write(100,96) xk(ik),xid(ik)
               enddo

               close(100)
               close(200)
            endif

          END SUBROUTINE READXI

        END MODULE PAIRING
