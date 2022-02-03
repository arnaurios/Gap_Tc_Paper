! ... ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! ... ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      PROGRAM BCS_PAIRING

        USE precision_definition
        USE physical_constants
        USE NN_INTERACTION
        USE NNN_INTERACTION
        USE mesh
        USE mesh_generator
        USE pairing
        USE interpolation
        USE KSPACE_MODULE
        USE read_write
        USE thermodynamical
        USE selfenergy

        IMPLICIT NONE
        REAL(long) :: fact
        INTEGER(ilong) ::ikf,itemp
        INTEGER(ilong) :: jc,ich
        CHARACTER(5) :: pwname
        LOGICAL :: cc,is_hf
        INTEGER :: pw_size

        ! ... MY CONSTANTS 2012
        hbc=197.3269718_long
        hbc2=hbc*hbc
        hbc3=hbc*hbc*hbc
        pi=4._long*datan(1._long)
        pi2=pi*pi

        xmassn=939.56536_long
        !xmassp=938.27203_long
        xmass=xmassn !(xmassn+xmassp)/2_long

        htm(1)=hbc2/xmassn
        htm(2)=hbc2/xmassn
        htmm=(htm(1)+htm(2))/2_long
        zi=(0._long,1._long)

        Nch=8

        deg=2._long

        iz1=-1
        iz2=-1

        jminhf=0
        jmaxhf=10

!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... GET THE DENSITIES AND ASYMMETRIES
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... GET THE DENSITIES AND ASYMMETRIES
        write(6,*) 'Give me rhoin, incr de rho and nrho'
!        read(5,*) rhoin,drho,nrho
!        write(6,*) rhoin,drho,nrho

        read(5,*) xkfin,dxkf,nxkf
        write(6,*) xkfin,dxkf,nxkf

        write(6,*) 'Give me tin, incr temp and ntemp'
        read(5,*) t_i,dt,nt
        write(6,*) t_i,dt,nt

        write(*,*) 'PWAVES? 1S0, 3PF2, etc'
        read (5, "(A)", SIZE=pw_size,ADVANCE='NO',EOR=10) p_waves
10      continue

        write(*,*) p_waves

! ... FIND THE JMIN AND JMAX FROM REQUESTS
        Jmin=7
        Jmax=0
        do jc=0,7
           do ich=1,3
              call AVFACNAME(iz1,iz2,jc,ich,pwname,cc)
              if( ss /=2 ) then
                 if( index(p_waves,trim(pwname)) .ne. 0 ) then
                    write(*,*) pwname
                    if( jc > Jmax ) Jmax=jc
                    if( jc < Jmin ) Jmin=jc
                 endif
              endif
           enddo
        enddo

        if(Jmin == 7 .and. Jmax == 0) then
           write(*,*) 'I did not understand partial wave input'
           stop
        endif
        write(*,'(a,i2)') 'Jmin=',Jmin
        write(*,'(a,i2)') 'Jmax=',Jmax

!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... FILE FORMATTING
        ALLOCATE( iunit1(Jmin:Jmax,4),iunit2(Jmin:Jmax,4),iunit3(Jmin:Jmax,4) )
        iunit1=200
        iunit2=300
        iunit3=400
        do jc=Jmin,Jmax
           do ich=1,3
              call AVFACNAME(iz1,iz2,jc,ich,pwname,cc)
              if( ss /= 2 .and. index(p_waves,trim(pwname)) .ne. 0 ) then
                 iunit1(jc,ich)=maxval(iunit1)+1
                 iunit2(jc,ich)=maxval(iunit2)+1
                 iunit3(jc,ich)=maxval(iunit3)+1

                 open(unit=iunit1(jc,ich),file='gap_kf_'//trim(pwname)//'.dat',status='unknown')
                 open(unit=iunit2(jc,ich),file='gap_'//trim(pwname)//'.dat',status='unknown')
                 open(unit=iunit3(jc,ich),file='2dgap_'//trim(pwname)//'.dat',status='unknown')

                 write(*,*) jc,ich,iunit1(jc,ich)
                 write(*,*) 'gap_'//trim(pwname)//'.dat'
              endif
           enddo
        enddo

        !     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! ... FILE FORMATTING
130     format(/,'#',3x,'rho [fm-3]',6x,'T [MeV]',6x,'kf [fm-1]',6x,'Gap(kf) [MeV]')
131     format(/,'#',3x,'rho [fm-3]',6x,'T [MeV]',6x,'kf [fm-1]',7x,'k [MeV]',7x,'av Gap(k) [MeV]' &
        ,7x,'Xi (k)',7x,'SigmaHF(k)',7x,'Momdis n(k)',7x,'Gap_L1',7x,'Gap_L2')

        do jc=Jmin,Jmax
           do ich=1,3
              call AVFACNAME(iz1,iz2,jc,ich,pwname,cc)
              if( ss /= 2 .and. index(p_waves,trim(pwname)) .ne. 0 ) then
                 if(iunit1(jc,ich) /= 200 ) write(iunit1(jc,ich),130)
                 if(iunit2(jc,ich) /= 300 ) write(iunit2(jc,ich),131)
                 if(iunit3(jc,ich) /= 400 ) write(iunit3(jc,ich),131)
              endif
           enddo
        enddo

        fact=hbc3
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... LOOP OVER DENSITY
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do 100 ikf=1,nxkf

           xxkf=xkfin + real(ikf-1,long)*dxkf
           rho=deg/6._long/pi2*xxkf**3

           ! ... FERMI MOMENTUM - NEUTRONS
           xkfermi=(6_long*pi2*rho/deg)**(1d0/3d0)
           ef=htm(1)*xkfermi**2/2_long
           cutoff=30_long*xkfermi

           call BUILDMESHES( xxkf )
           if(ikf==1) write(*,'(100i3)') N

           if( .not.ALLOCATED(spe) ) ALLOCATE( spe(Nkmesh) )
           if( .not.ALLOCATED(dgap) ) ALLOCATE( dgap(Nkmesh) )
           if( .not.ALLOCATED(vHF) ) ALLOCATE( vHF(Nkmeshhf,Nkmeshhf,jminhf:jmaxhf,Nch) )
           dgap=0._long

           ! ... CALL NNFORCE IN THE SAME MESH
           if( .not.ALLOCATED(vNN)) ALLOCATE( vNN(Nkmesh,Nkmesh,Jmin:Jmax,Nch) )
           is_hf=.false.
           call NNPOT(is_hf,Jmin,Jmax,xk,v2N)

           vNN=v2N
           call WRITEVNN

           ! ... call NNFORCE FOR HARTREE-FOCK CALCULATION
           ! ... ONLY COMPUTES DIAGONAL ELEMENTS
           is_hf=.true.
           call NNPOT(is_hf,jminhf,jmaxhf,xkhf,v2HF)

           if(ikf == 1) then
              spe=xkmesh**2*htm(1)/2_long
              dgap=0._long
           endif

!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... LOOP OVER TEMPERATURE
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           do 200 itemp=1,nt

              t=t_i + real(itemp-1,long)*dt

              call ADJUST_MU( dgap )
              write(*,*) 'HF',xmu

              is_hf=.true.
              call N3POT_PNM(is_hf,jminhf,jmaxhf,xxkf,xkhf,wHF)
              vHF=v2HF+wHF/2._long
              call HARTREE_FOCK(vHF,sigma_hf)

              spe=xkmesh**2*htm(1)/2_long+sigma_hf

              call ADJUST_MU( dgap )
              write(*,*) spkf
              write(*,*) xmu


              is_hf=.false.
              call N3POT_PNM(is_hf,Jmin,Jmax,xxkf,xk,v3N)
              vNN=v3N
              call WRITEVNN

              vNN=v2N + v3N!/2_long
              call WRITEVNN

              vNN=vNN*fact
              CALL P_GAP

!              spe=xkmesh**2*htm(1)/2_long+sigma_hf

!              call ADJUST_MU( dgap )
!              write(*,*) spkf
!              write(*,*) xmu

200        enddo ! LOOP OVER TEMPERATURE

           ! ... TWO BODY EFFECTIVE INTERACTION

100     enddo ! LOOP OVER DENSITY

      END PROGRAM BCS_PAIRING
