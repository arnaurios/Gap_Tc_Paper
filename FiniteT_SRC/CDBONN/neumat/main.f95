! ... ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     PROGRAM FOR BEYOND BCS PAIRING IN NEUTRON MATTER
!
! ... READS CHI DENOMINATOR FROM g2a_rhoxx FOLDER (wim_T0.dat)
! ... COMPUTES THE PAIRING GAP
! ... CODE IS SENSITIVE TO a) MESH AROUND KF and b) INTERPOLATION OF CHI
! ... ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      PROGRAM BEYOND_BCS_PAIRING

        USE precision_definition
        USE physical_constants
        USE NN_interaction
        USE mesh
        USE mesh_generator
        USE pairing
        USE interpolation
        USE KSPACE_MODULE
        USE read_write
        USE thermodynamics

        IMPLICIT NONE
        REAL (long) :: xf,xi
        REAL (long), DIMENSION(100) :: xrho
        REAL (long), DIMENSION(100) :: xtemp
        INTEGER(ilong) ::irho,itemp,ireg
        INTEGER(ilong) :: Ni,Nf,jc,ich,Nreg
        INTEGER(ilong), ALLOCATABLE, DIMENSION(:) :: N
        REAL(long), ALLOCATABLE, DIMENSION(:) :: xmsh
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

!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... GET THE DENSITIES AND ASYMMETRIES
        write(6,*) 'Give me rhoin, incr de rho and nrho'
        !        read(5,*) rhoin,drho,nrho
        read(5,*) nrho,(xrho(irho),irho=1,nrho)
        write(6,*) nrho,(xrho(irho),irho=1,nrho)

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

! ... READ XI? IF FALSE, USE FREE SP SPECTRA
        write(*,*) 'Read XI'
        read(*,*) xi_inp

! ... READ TEMPERATURE
        write(*,*) 'Read Temperature'
        read(5,*) ntemp,(xtemp(itemp),itemp=1,ntemp)
        write(*,*) ntemp,(xtemp(itemp),itemp=1,ntemp)

!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        iz1=-1
        iz2=-1

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
130     format(/,'#',3x,'kf [fm-1]',7x,'T [MeV]',7x,'Gap(kf) [MeV]')
131     format(/,'#',3x,'rho [fm-3]',6x,'T [MeV]',6x,'kf [fm-1]',7x,'k [fm-1]',7x,'Gap(k) [MeV]',7x,'auxgap')

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

!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! ... ALLOCATE ARRAYS
        Nreg=6 ! +1 on tangential
        ALLOCATE( N(Nreg+1),xmsh(Nreg+1) )

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

        Nkmesh=sum(N(1:Nreg+1))
        ALLOCATE( vNN(Nkmesh,Nkmesh,0:Jmax,Nch) )
        ALLOCATE( xkmesh(Nkmesh),wkmesh(Nkmesh) )
        ALLOCATE( xk(Nkmesh),wk(Nkmesh) )
        Nmsh=Nkmesh
        cmaxk=1e6

! ... FOR 2D DATA, INTERPOLATE TO kf x kout mesh
        Nkout=100
        ALLOCATE( xkout(Nkout),gkout(Nkout) )
        xi=real(0_long,long)
        xf=real(10_long,long)
        CALL LINEAR(xi,xf,Nkout,xkout,gkout)

        ! ... START LOOPS - DENSITY
        do 100 irho=1,nrho

           rho=xrho(irho) !rhoin + real(irho-1,long)*drho

! ... FERMI MOMENTUM - NEUTRONS
           xkfermi=(6_long*pi2*rho/deg)**(1d0/3d0)
           ef=htm(1)*xkfermi**2/2_long
           cutoff=30_long*xkfermi

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

           if(irho==1) write(*,'(100i3)') N

           ! ... CALL NNFORCE IN THE SAME MESH
           xk=xkmesh*hbc
           call NNPOT

           fact=hbc3
           vNN=vNN*fact
!           call WRITEVNN

           ! ... READ WIM_T0 OR WIM_TF

           ! ... READ WIM_T0 OR WIM_TF
           do itemp=1,ntemp
              temp=xtemp(itemp)

              if( .not. ALLOCATED(xid) ) ALLOCATE( xid(Nkmesh) )
              if( xi_inp ) call READXI
              if( .not.xi_inp ) xid=0_long

              ! ... READ WIM_T0 OR WIM_TF
              CALL P_GAP

           enddo ! LOOP OVER TEMPERATURE

100     enddo ! LOOP OVER DENSITY

   END PROGRAM BEYOND_BCS_PAIRING
