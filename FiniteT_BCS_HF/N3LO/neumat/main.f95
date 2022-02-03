! ... ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     PROGRAM FOR BCS PAIRING IN NEUTRON MATTER
! ... ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      PROGRAM BCS_PAIRING

        USE precision_definition
        USE physical_constants
        USE NN_interaction
        USE mesh
        USE mesh_generator
        USE pairing
        USE interpolation
        USE KSPACE_MODULE
        USE read_write
        USE thermodynamical
        USE selfenergy

        IMPLICIT NONE
        REAL(long) :: fact,drho,rhoin
        INTEGER(ilong) ::ikf,itemp
        INTEGER(ilong) :: jc,ich
        CHARACTER(5) :: pwname
        LOGICAL :: cc
        INTEGER :: pw_size

        ! ... MY CONSTANTS 2012
        hbc=197.3269718_long
        hbc2=hbc*hbc
        hbc3=hbc*hbc*hbc
        pi=4._long*datan(1._long)
        pi2=pi*pi

        xmassn=939.56536_long
        xmassp=938.27203_long
        xmass=(xmassn+xmassp)/2._long

        htm(1)=hbc2/xmassn
        htm(2)=hbc2/xmassp
        htmm=(htm(1)+htm(2))/2._long
        zi=(0._long,1._long)

        Nch=8

        deg=2._long

        iz1=-1
        iz2=-1

        jminhf=0
        jmaxhf=5

!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... GET THE DENSITIES AND ASYMMETRIES
        write(6,*) 'Give me rhoin, incr de rho and nrho'
        !read(5,*) rhoin,drho,nxkf
        !write(6,*) rhoin,drho,nxkf
        read(5,*) xkfin,dxkf,nxkf
        write(6,*) xkfin,dxkf,nxkf

        write(6,*) 'Give me tin, incr temp and ntemp'
        read(5,*) t_i,dt,nt
        write(6,*) t_i,dt,nt


! ... READ MAXIMUM J
!        write(*,*) 'Jmax?'
!        read(*,*) jmin,jmax
!        write(*,*) jmin,jmax

        write(*,*) 'PWAVES? 1S0, 3PF2, etc'
        read (5, "(A)", SIZE=pw_size,ADVANCE='NO',EOR=10) p_waves
10      continue

        write(*,*) p_waves
        write(*,*) index(p_waves,'1S0')

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
                if( iunit1(jc,ich) /= 200 ) write(iunit1(jc,ich),130)
                if( iunit2(jc,ich) /= 300 ) write(iunit2(jc,ich),131)
                if( iunit3(jc,ich) /= 400 ) write(iunit3(jc,ich),131)
              endif
           enddo
        enddo


!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... LOOP OVER DENSITY
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do 100 ikf=1,nxkf

           xxkf=xkfin + real(ikf-1,long)*dxkf
           rho=deg/6._long/pi2*xxkf**3
!           rho=rhoin+real(ikf-1,long)*drho

           ! ... FERMI MOMENTUM - NEUTRONS
           xkfermi=(6_long*pi2*rho/deg)**(1d0/3d0)
           xxkf=xkfermi

           ef=htm(1)*xkfermi**2/2_long
           cutoff=30_long*xkfermi
           write(*,*) xkfermi

           call BUILDMESHES( xxkf )
           if(ikf==1) write(*,'(100i3)') N

           if( .not.ALLOCATED(spe) ) ALLOCATE( spe(Nkmesh) )
           if( .not.ALLOCATED(dgap) ) ALLOCATE( dgap(Nkmesh) )

           call NNPOT(Jmin,Jmax,xk,vNN)
           fact=hbc3
           vNN=vNN*fact

           call NNPOT(Jminhf,Jmaxhf,xkhf,vHF)
           vHF=vHF*fact

           if(ikf == 1) then
              spe=xkmesh**2*htm(1)/2_long
              dgap=0._long
           endif

           ! ... START LOOPS - DENSITY
           do 200 itemp=1,nt

              t=t_i + real(itemp-1,long)*dt

              call ADJUST_MU( dgap )
              write(*,*) spkf
              write(*,*) xmu

              write(*,*) 'HF'
              call HARTREE_FOCK(vHF,sigma_hf)
              write(*,*) 'end HF'
              spe=xkmesh**2*htm(1)/2_long+sigma_hf

              call ADJUST_MU( dgap )
              write(*,*) spkf
              write(*,*) xmu

              CALL P_GAP

200        enddo ! LOOP OVER TEMPERATURE
100     enddo ! LOOP OVER  DENSITY
      END PROGRAM BCS_PAIRING
