! ... NN interaction passed around
      MODULE NNN_INTERACTION

        USE precision_definition
        REAL(long), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: vNNN

        REAL (long), ALLOCATABLE, DIMENSION(:,:,:,:) :: v3N,wHF

        real(long), ALLOCATABLE, DIMENSION(:,:,:) :: q2
        real(long), ALLOCATABLE, DIMENSION(:,:) :: p2
        real(long), ALLOCATABLE, DIMENSION(:) :: angle, wangle
        real(long), ALLOCATABLE, DIMENSION(:) :: xk3,wk3,xnk,reg_int

        real(long) :: gA,fPI,c1,c3,c4,cD,cE,mN,mPI,lchi,xlmax
        INTEGER(ilong) :: Nangle,Nk3,nexp,ireg
        LOGICAL :: internalstep

        CONTAINS

!....cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... LECs & REGULATOR CHOICES FOR 3NF SUBROUTINE
!....cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          ! CONSTANTS HERE FOR NNLOpt AND NON-LOCAL REGULATOR AS IN
          ! PHYSICAL REVIEW C 89, 014319 (2014)
          SUBROUTINE LECS
            USE precision_definition
            USE physical_constants
            IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... Chiral constants
            gA=1.29d0
            fPI=92.4d0 !MeV
            mN=xmass !MeV
            mPI=138.04d0 !MeV
            lchi=700.d0 !MeV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... LECs
            c1=-0.81d0*1.E-3 !MeV^-1
            c3=-3.2d0*1.E-3 !MeV^-1

            c4=0.d0 !MeV^-1
            cD= 0d0
            cE= 0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! ... REGULATOR CUTOFF
            xlmax=500.d0 !MeV
            ! ... INTERNAL REGULATOR ireg=1 OR EXTERNAL REGULATOR ireg=2
            ireg=1
            ! ... nexp OF REGULATOR FUNCTION
            !           nexp=4
            nexp=2
            ! ... INTERNAL n(k) IS STEP FUNCTION (T) OR CORRELATED (F)
            internalstep=.false.

          END SUBROUTINE LECS

!....cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!....Density dependent effective NN interaction from                      c
!....N2LO three-body forces in Nuclear Matter.                            c
!....PRC 81, 024002 (2010)                                                c
!.... NEW VERSION JAN 2016 TO ACCOUNT FOR SEVERAL IMPROVEMENTS
!....cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          SUBROUTINE N3POT_PNM( ishf,jmn,jmx,xxkf,xk,vv )

            USE precision_definition
            USE physical_constants
            USE mesh_generator
            USE fermibose
            USE thermodynamical, ONLY : xmu,t,spe,momdis
            USE selfenergy, ONLY : dgap
            USE interpolation
            USE NN_INTERACTION, ONLY : iz1,iz2
            USE mesh, ONLY: xkmesh

            IMPLICIT NONE
            LOGICAL, INTENT(IN) :: ishf
            real(long), INTENT(IN) :: xxkf
            INTEGER(ilong), INTENT(IN) :: jmn,jmx
            REAL(long), ALLOCATABLE, DIMENSION(:),INTENT(IN) :: xk
            REAL(long), ALLOCATABLE, DIMENSION(:,:,:,:), INTENT(OUT) :: vv

            INTEGER(ilong) :: ik,jk,ian,Nk1,Nk2,jki,jkf
            INTEGER(ilong) :: Ni,Nf,i,ik3,Nk

            real(long) :: v1,v2,v3,reg,rho_int,rho
            real(long) :: G0,Gst,G2st,G1,G1st,G2,G3,ala,fac
            real(long) :: gamma0,gamma1,gamma2,gamma3
            real(long) :: xi,xf,xkf,p,pp2,qq2,q,xk1,xk2

            real(long), ALLOCATABLE, DIMENSION(:) :: spe3,dgp3

            real(long), ALLOCATABLE, DIMENSION(:) :: xaux,waux

            real(long), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: vNNN_aux

!            write(*,*) 'Loading the averaged NNN potential in NM...'

            ! ... CONSTANTS FOR DIFFERENT TERMS (all potential terms are in fm)
            call LECs

            v1=hbc3*hbc*gA**2*mN/(8.d0*pi*fPI**4)
            v2=hbc*gA**2*mN/(32.d0*pi*pi2*fPI**4)
            v3=hbc*gA**2*mN/(64.d0*pi2*pi*fPI**4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... MOMENTUM MESH FOR INTEGRALS OVER INTERNAL MOMENTA
            xkf=xxkf*hbc
            rho=deg/6._long/pi2*xxkf**3
            Nk=size( xk )

            if( internalstep ) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! ... INTERAL STEP FUNCTION INTEGRATION
               write(*,*) 'Step function averaged 3BF integration'

               Nk3=256
               ALLOCATE( xk3(Nk3),wk3(Nk3),xnk(Nk3),reg_int(Nk3) )

               ! ... gaussian set of points for momenta xk3 from 0 to kF
               xi=0.d0
               xf=xkf
               call GAUSS(xi,xf,Nk3,xk3,wk3)
               xnk=1d0

            else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               ! ... INTEGRAL OF CORRELATED MOMENTUM DISTRIBUTION
               ! ... NEED TO COMPUTE MOMENTUM DISTRIBUTION FROM ASF(k)
               write(*,*) 'n(k) averaged 3BF inegration'

               Nk1=32
               Nk2=32
               Nk3=2*Nk1+Nk2

               ALLOCATE( xk3(Nk3),wk3(Nk3) )

               Ni=1
               Nf=Nk1
               xi=0.d0
               xf=xkf
               call GAUSS(xi,xf,Nk1,xk3(Ni:Nf),wk3(Ni:Nf))

               Ni=Nf+1
               Nf=Ni+Nk2-1
               xi=xf
               xf=4.d0*xkf/2.d0
               call GAUSS(xi,xf,Nk2,xk3(Ni:Nf),wk3(Ni:Nf))

               Ni=Nf+1
               Nf=Ni+Nk1-1
               xi=xf
               xf=100.d0*xkf
               ALLOCATE( xaux(Nk1),waux(Nk1) )
               call GAUSS(0d0,1d0,Nk1,xaux,waux)
               fac= (xf-xi)/tan( pi/2d0*xaux(Nk1) )
               do i=1,Nk1
                  xk3(Ni+i-1) = xi + fac*tan( pi/2d0*xaux(i)  )
                  wk3(Ni+i-1) = waux(i)*fac*pi/2d0/(cos(pi/2d0*xaux(i))**2)
               enddo
               DEALLOCATE(xaux,waux)

               ! ... LOOP OVER KMESH
               ! ... interpolation of momentum distribution to mesh for integrals

               ALLOCATE( spe3(Nk3),dgp3(Nk3) )

               call LIN_INT(xkmesh*hbc, spe,xk3,spe3)
               call LIN_INT(xkmesh*hbc,dgap,xk3,dgp3)

               ALLOCATE( xnk(Nk3),reg_int(Nk3) )

               xnk = momdis(t,xmu,spe3,dgp3)

            endif

        ! ... TEST DENSITY VALUE FOR VALIDITY OF xnk
100         format(' TEST RHO IN 3NF ROUTINE rho=',f8.6,1x,2L)
            write(*,100) deg/2.d0/pi2*sum(xk3**2*xnk*wk3)/hbc3,ishf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! ... REGULATOR
            ! ... internal regulator function (x3lmax from modules)
            reg_int=1.d0 ! no internal regulator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! ... ANGULAR INTEGRATIONS
            Nangle=24
            ALLOCATE(angle(Nangle), wangle(Nangle))
            xi=-1d0
            xf= 1d0
            call GAUSS(xi,xf,Nangle,angle,wangle)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! ... ALLOCATE AUXILIARY POTENTIAL MATRIX
            ALLOCATE(vNNN_aux(Nk,Nk,Nangle,5,2,7))
            vNNN_aux=0.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! ... symmetric extrapolation for off-shell elements
            ! ... see PRC79,054331(2009)
            ALLOCATE(p2(Nk,Nk))
            FORALL(ik=1:Nk,jk=1:Nk)
               p2(ik,jk)=(xk(ik)**2+xk(jk)**2)/2.d0
            END FORALL

            ALLOCATE( q2(Nk,Nk,Nangle) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! ... BUILD FUNCTIONS
            ! ... ANGULAR INTEGRAL
            do ian=1,Nangle

               ! FIRST RELATIVE MOMENTUM INTEGRAL
               do ik=1,Nk
                  xk1=xk(ik)
                  ! SECOND RELATIVE MOMENTUM INTEGRAL
                  if(ishf) then
                     jki=ik
                     jkf=ik
                  else
                     jki=1
                     jkf=Nk
                  endif

                  do jk=jki,jkf
                     xk2=xk(jk)

! ... MOMENTUM TRANSFER SQUARED
                     q2(ik,jk,ian)=xk(ik)**2 + xk(jk)**2 &
                          - 2.d0*xk(ik)*xk(jk)*angle(ian)
                     qq2=q2(ik,jk,ian)
                     q=sqrt(q2(ik,jk,ian))

! ... TOTAL MOMENTUM
                     p=sqrt(p2(ik,jk))
                     pp2=p2(ik,jk)

! ... REGULATOR & DENSITY THAT APPEARS IN INTERACTION TERMS
                     if(ireg==1) then
                        reg = EXP( -2.d0*(p2(ik,jk)/xlmax**2)**nexp )
                        rho_int=rho
                     elseif(ireg==2) then
                        reg=1.d0
                        reg_int=exp( -2.d0*( (xk3**2/3.d0 + pp2 )/xlmax**2)**nexp )
                        rho_int=deg/2.d0/pi2*sum(xk3**2*xnk*wk3*reg_int)/hbc3
                        if(ian==50.and.ik==25.and.jk==25) write(*,100) rho_int
                     endif

                     call GAMMA( p,pp2,q,qq2, &
                          gamma0,gamma1,gamma2,gamma3,G0,Gst,G2st,G1,G1st,G2,G3 )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... TERMS OF INTERACTION
! ... FIRST TERM
                 ! Iso-vector tensor
                        vNNN_aux(ik,jk,ian,3,2,1) = reg*v1*rho_int &
                             *(2.d0*c1*mPI**2 + c3*q2(ik,jk,ian)) &
                             /(mPI**2 + q2(ik,jk,ian))**2
! ... SECOND TERM
                 ! Iso-vector tensor
                        vNNN_aux(ik,jk,ian,3,2,2) = reg*v2/( mPI**2+q2(ik,jk,ian) )* &
                             ( -4.d0*c1*mPI**2*(gamma0+gamma1) &
                             -(c3+c4)*(q2(ik,jk,ian)*(gamma0+2.d0*gamma1+gamma3) &
                             +4.d0*gamma2 ) &
                             +4.d0*c4*(pi2*rho_int*hbc3-mPI**2*gamma0))

! ... THIRD TERM ... A 1/3 TERM NEEDED IN NEUTRON MATTER
                   ! Iso-singlet central
                   ala=48d0/deg*pi2*rho_int*hbc3 &
                        -12d0*(2d0*mPI**2+q2(ik,jk,ian))*gamma0 &
                        -6d0*q2(ik,jk,ian)*gamma1 &
                        +3d0*(2d0*mPI**2+q2(ik,jk,ian))**2*G0

                   vNNN_aux(ik,jk,ian,1,1,3) = reg*v3*( &
                        -12d0*c1*mPI**2*(2d0*gamma0-(2d0*mPI**2+q2(ik,jk,ian))*G0) &
                        -c3*ala ) &
                        *(1d0/2d0+real(iz1+iz2)/12d0)

                   ! Iso-singlet SO
                   ala= 2d0*(gamma0+gamma1) &
                        -(2d0*mPI**2+q2(ik,jk,ian))*(G0+2d0*G1)

                   vNNN_aux(ik,jk,ian,4,1,3) = reg*v3*( &
                        -3d0*c3*( ala )  &
                        -12d0*c1*mPI**2*(G0 + 2d0*G1) ) &
                        *(1d0/2d0+real(iz1+iz2)/12d0)

                     enddo ! Loop over second momentum
                  enddo ! Loop over first momentum
               enddo ! Loop over angles

               ALLOCATE(vNNN(Nk,Nk,Nangle,5,2))
               vNNN=0.d0

! ... Sum of all density-dependent terms into two isospins
               vNNN(:,:,:,:,1)=vNNN_aux(:,:,:,:,1,3)+ &
                    vNNN_aux(:,:,:,:,1,5)+vNNN_aux(:,:,:,:,1,6)

               vNNN(:,:,:,:,2)=vNNN_aux(:,:,:,:,2,1)+ &
                    vNNN_aux(:,:,:,:,2,2)+vNNN_aux(:,:,:,:,2,3)+ &
                    vNNN_aux(:,:,:,:,2,4)+vNNN_aux(:,:,:,:,2,5)

               DEALLOCATE( xk3,wk3,xnk,reg_int )
               DEALLOCATE( vNNN_aux )

               if( any(isnan(vNNN)) ) then
                  write(*,*) 'NaN'
                  stop
               endif

               ! ... PROJECT INTO PARTIAL WAVES
               call NNNPOT_PW( ishf,jmn,jmx,xk,vv )

               DEALLOCATE( vNNN )

          END SUBROUTINE N3POT_PNM

!... **********************************************************
!...  Calculates partial wave matrix elements for the reduced *
!...  form of the 3-body part of the N2LO chiral interaction  *
!... **********************************************************
          SUBROUTINE NNNPOT_PW( ishf,jmn,jmx,xk,v3v )

            USE precision_definition
            USE physical_constants
            USE NN_INTERACTION, only : iz1,iz2,Nch

            IMPLICIT NONE
            LOGICAL, INTENT(IN) :: ishf
            REAL(long), ALLOCATABLE, DIMENSION(:),INTENT(IN) :: xk
            INTEGER(ilong), INTENT(IN) :: jmn,jmx
            REAL(long), ALLOCATABLE, DIMENSION(:,:,:,:), INTENT(OUT) :: v3v


            real(long), ALLOCATABLE, DIMENSION(:,:) :: pleg
            real(long), ALLOCATABLE, DIMENSION(:) :: aux1,aux2,paux1,paux2
            INTEGER(ilong) :: ik,jk,jj,ll,it,s,ich,itii,jki,jkf
            INTEGER :: ii,Nk

            real(long), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: uNNN
            real(long) :: xf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... PARTIAL WAVE EXPANSION
! ... allocate auxiliary variables
            ALLOCATE( aux1(Nangle), aux2(Nangle), paux1(Nangle), paux2(Nangle) )
            ALLOCATE( pleg(0:jmx+1,Nangle) )
            aux1=0.d0
            aux2=0.d0
            paux1=0.d0
            paux2=0.d0
            pleg=0.d0

            ! ... allocate potential matrix
            Nk=size(xk)
            if(.not. ALLOCATED(v3v)) ALLOCATE(v3v(Nk,Nk,jmn:jmx,Nch))    !fm

            ALLOCATE(uNNN(Nk,Nk,Nangle,5,0:1)) !fm
            uNNN=0.d0
            do it=0,1
               uNNN(:,:,:,:,it)=vNNN(:,:,:,:,1) &
                    +dble(4*it-3)*vNNN(:,:,:,:,2)
            enddo
! ... CALL LEGENDRE POLYNOMIALS AT ALL ANGLES
            ii=jmx+1
            call LPN(ii,angle,pleg)

            v3v=0.d0

            iz1=-1
            iz2=-1
        ! ... conditions on isospin
            if(iz1*iz2 == 1) then
               itii=1
            else
               itii=0
            endif

! ... momentum loop
            do ik = 1,Nk
               if(ishf) then
                  jki=ik
                  jkf=ik
               else
                  jki=1
                  jkf=Nk
               endif
               do jk=jki,jkf

              ! ... in the partial wave analysis, setting p2 this way resembles more the full non diagonal calculation of matrix elements (tested with Christian)
!        p2(ik,jk)=xk(ik)*xk(jk)

! ... filling channels
                  do it=itii,1
                     do ich=1,6
                        ! ... singlet
                        if(ich == 1) then
                           s=0
                           do jj=jmn,jmx
                              ll=jj
                              if(mod(s+ll+it,2) == 0) cycle
                              aux1 = uNNN(ik,jk,:,1,it)-3.d0*uNNN(ik,jk,:,2,it) &
                               -q2(ik,jk,:)*uNNN(ik,jk,:,3,it) &
                               +p2(ik,jk)**2*(angle**2-1.d0)*uNNN(ik,jk,:,5,it)
                              paux1= pleg(jj,:)
                              v3v(ik,jk,jj,ich)= 0.5d0*sum( aux1*paux1*wangle )
                           enddo
! ... triplet
                        elseif(ich == 2) then
                           s=1
                           do jj=max(1,jmn),jmx
                              ll=jj
                              if(mod(s+ll+it,2) == 0) cycle
                              aux1 = 2.d0*p2(ik,jk)*(uNNN(ik,jk,:,4,it) &
                                   -uNNN(ik,jk,:,3,it) &
                               +p2(ik,jk)*angle*uNNN(ik,jk,:,5,it) )
                              paux1 = pleg(jj+1,:)+pleg(jj-1,:)

                              aux2 =  uNNN(ik,jk,:,1,it)+uNNN(ik,jk,:,2,it) &
                                   +2.d0*p2(ik,jk)*(1.d0+angle)*uNNN(ik,jk,:,3,it) &
                                   -4.d0*p2(ik,jk)*angle*uNNN(ik,jk,:,4,it) &
                                   -p2(ik,jk)**2*(3.d0*angle**2+1.d0)*uNNN(ik,jk,:,5,it)
                              paux2 = pleg(jj,:)

                              v3v(ik,jk,jj,ich)= 0.5d0* &
                                   sum( (aux1*paux1 + aux2*paux2)*wangle )
                           enddo
                           ! ... coupled triplet V++
                        elseif(ich == 4)then
                           s=1
                           do jj=jmn,jmx
                              ll=jj+1
                              if(mod(s+ll+it,2) == 0) cycle
                              aux1 = 2.d0*p2(ik,jk)*(uNNN(ik,jk,:,4,it) &
                                   +1.d0/dble(2*jj+1) &
                                   *(uNNN(ik,jk,:,3,it)-p2(ik,jk)* &
                                   angle*uNNN(ik,jk,:,5,it)))
                              paux1 = pleg(jj,:)

                              aux2 = uNNN(ik,jk,:,1,it)+uNNN(ik,jk,:,2,it)+ &
                                   p2(ik,jk)* &
                                   ( p2(ik,jk)*(1-angle**2)*uNNN(ik,jk,:,5,it) &
                                   - 2.d0*angle*uNNN(ik,jk,:,4,it) &
                                   + 2.d0/dble(2*jj+1) &
                                   *(p2(ik,jk)*uNNN(ik,jk,:,5,it)-uNNN(ik,jk,:,3,it)))
                              paux2 = pleg(jj+1,:)

                              v3v(ik,jk,jj,ich)= 0.5d0* &
                                   sum( (aux1*paux1 + aux2*paux2)*wangle )
                              !                        if(jj == 0)then
                              !                           v3N(ik,jk,jj,2)=v3N(ik,jk,jj,ich)
                              !                           v3N(ik,jk,jj,4)=0.d0
                              !                        endif
                           enddo
! ... coupled triplet V--
                        elseif(ich == 3)then
                           s=1
                           do jj=max(1,jmn),jmx
                              ll=jj-1
                              if(mod(s+ll+it,2).eq.0)cycle
                              aux1 = 2.d0*p2(ik,jk)*(uNNN(ik,jk,:,4,it) &
                                   -1.d0/dble(2*jj+1) &
                                   *(uNNN(ik,jk,:,3,it)-p2(ik,jk)* &
                                   angle*uNNN(ik,jk,:,5,it)))
                              paux1 = pleg(jj,:)

                              aux2 = uNNN(ik,jk,:,1,it)+uNNN(ik,jk,:,2,it)+ &
                                   p2(ik,jk)* &
                                   ( p2(ik,jk)*(1-angle**2)*uNNN(ik,jk,:,5,it) &
                                   - 2.d0*angle*uNNN(ik,jk,:,4,it) &
                                   - 2.d0/dble(2*jj+1) &
                                   *(p2(ik,jk)*uNNN(ik,jk,:,5,it)-uNNN(ik,jk,:,3,it)))
                              paux2 = pleg(jj-1,:)

                              v3v(ik,jk,jj,ich)= 0.5d0* &
                                   sum( (aux1*paux1 + aux2*paux2)*wangle )
                           enddo
                           ! ... coupled triplet V-+
                        elseif(ich == 5)then
                           s=1
                           do jj=max(1,jmn),jmx
                              ll=jj-1
                              if(mod(s+ll+it,2) == 0)cycle
                              aux1 = uNNN(ik,jk,:,3,it)-p2(ik,jk)*uNNN(ik,jk,:,5,it)
                              paux1 = pleg(jj+1,:)

                              aux2= (dble(2*jj)-angle*dble(2*jj+1)) &
                                   *uNNN(ik,jk,:,3,it)+p2(ik,jk)*angle &
                                   *uNNN(ik,jk,:,5,it)
                              paux2=pleg(jj,:)

                              v3v(ik,jk,jj,ich)=(sqrt(dble(jj+1)))*p2(ik,jk)/  &
                                   (sqrt(dble(jj))*dble(2*jj+1))* &
                                   sum( (aux1*paux1 + aux2*paux2)*wangle )
                           enddo
                           ! ... coupled triplet V+-
                        elseif(ich == 6)then
                           s=1
                           do jj=max(1,jmn),jmx
                              ll=jj+1
                              if(mod(s+ll+it,2) == 0)cycle
                              aux1 = uNNN(ik,jk,:,3,it)-p2(ik,jk)*uNNN(ik,jk,:,5,it)
                              paux1 = pleg(jj+1,:)

                              aux2= (dble(2*jj)-angle*dble(2*jj+1)) &
                                   *uNNN(ik,jk,:,3,it)+p2(ik,jk)*angle &
                                   *uNNN(ik,jk,:,5,it)
                              paux2=pleg(jj,:)

                              v3v(ik,jk,jj,ich)=(sqrt(dble(jj+1)))*p2(ik,jk)/  &
                                   (sqrt(dble(jj))*dble(2*jj+1))* &
                                   sum( (aux1*paux1 + aux2*paux2)*wangle )
                           enddo
                        endif

                     enddo ! ... loop filling channels
                  enddo ! ... loop isospin
               enddo  ! ... loop second momentum
            enddo ! ... loop first momentum

            WHERE( abs(v3v) < 1e-30 ) v3v=0.d0

            ! ... TO HAVE MATRIX ELEMENTS IN UNITS OF MEV-2, as 2NF FORCE
            xf=4d0*pi/xmass/(2d0*pi**2)/hbc !2/m/pi
            v3v=v3v*xf

!            write(*,*) 'DONE WITH 3B MATRIX ELEMENTS'

            DEALLOCATE( aux1,aux2,paux1,paux2,pleg )
            DEALLOCATE( p2,q2 )

            DEALLOCATE( uNNN )
            DEALLOCATE( angle,wangle )

          END SUBROUTINE NNNPOT_PW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ...  Computes Gamma functions for density dependent 2B N2LO
    ! ... INPUT IS xk3,xkw3 and momentum distribution xnk
    ! ... INPUT IS xkxxx
    ! ... OUTPUT GAMMA0 TO GAMMA3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          SUBROUTINE GAMMA( p,pp2,q,qq2,gamma0,gamma1,gamma2,gamma3,G0,Gst,G2st,G1,G1st,G2,G3 )

            USE precision_definition
            USE physical_constants
            IMPLICIT NONE
            real(long), INTENT(IN) :: p,pp2,q,qq2
            real(long), INTENT(OUT) :: gamma0,gamma1,gamma2,gamma3
            real(long), INTENT(OUT) :: G0,Gst,G2st,G1,G1st,G2,G3
            real(long), ALLOCATABLE, DIMENSION(:) :: theta,aux,ak

            gamma0=0.d0
            gamma1=0.d0
            gamma2=0.d0
            gamma3=0.d0

            G0=0.d0
            Gst=0.d0
            G2st=0.d0
            G1=0.d0
            G1st=0.d0
            G2=0.d0
            G3=0.d0

            ! ... Allocate auxiliary variable for gamma integrals
            ALLOCATE(theta(Nk3),aux(Nk3),ak(Nk3))

            ! ... define auxiliary variable for gamma integrals
            theta=xk3*xnk &
                 *log( ((p+xk3)**2 + mPI**2) / ((p-xk3)**2+mPI**2) )

      ! ... Define gamma integrals
            gamma0=1.d0/2.d0/p*sum(wk3*reg_int*theta)

            gamma1=1.d0/2.d0/p**2*sum(wk3*reg_int*( &
                 2.d0*xk3**2*xnk - (p**2+xk3**2+mPI**2)/(2.d0*p)*theta &
                 ))

            gamma2=1.d0/4.d0/p*sum(wk3*reg_int*( &
                 xk3**2*xnk*(p**2+xk3**2+mPI**2)/p &
                 +(1.d0-( (p**2+xk3**2+mPI**2)/(2.d0*p*xk3))**2)*xk3**2*theta &
                 ))

            gamma3=1.d0/4.d0/p**3*sum(wk3*reg_int*( &
                 -3.d0*xk3**2/p*(p**2+xk3**2+mPI**2)*xnk &
                 +(3*((p**2+xk3**2+mPI**2)/(2.d0*p*xk3))**2-1.d0)*xk3**2*theta &
                 ))

            ! ... AUXILIARY INTEGRANDS
            ak=(mPI**2+(xk3+p)**2)*(mPI**2+(xk3-p)**2)
            aux= xk3/sqrt( ak+qq2*xk3**2 )* &
                 log( (q*xk3+sqrt(ak+qq2*xk3**2))/sqrt(ak) )

            G0=  2.d0/q*sum(wk3*reg_int* xnk*aux)
            Gst= 2.d0/q*sum(wk3*reg_int* xnk*aux*xk3**2)
            G2st=2.d0/q*sum(wk3*reg_int* xnk*aux*xk3**4)

            ! ... USE RECURSION RELATIONS
            G1=(gamma0-(mPI**2+pp2)*G0-Gst)/(4.d0*pp2-qq2)
            G1st=(3.d0*gamma2+pp2*gamma3  &
                 - (mPI**2+pp2)*Gst-G2st) / (4.d0*pp2-qq2)
            G2=(mPI**2+pp2)*G1+Gst+G1st
            G3=(gamma1/2.d0-2.d0*(mPI**2+pp2)*G1-2.d0*G1st-Gst) &
                 /(4.d0*pp2-qq2)

            DEALLOCATE(theta,aux,ak)
          END SUBROUTINE GAMMA

!       ========================================================
!       Purpose: This program computes the Legendre polynomials
!                Pn(x) and their derivatives Pn'(x) using
!                subroutine LPN
!       Input :  x --- Argument of Pn(x)
!                n --- Degree of Pn(x) ( n = 0,1,...)
!       Output:  PN(n) --- Pn(x)
!                PD(n) --- Pn'(x)
!       Example:    x = 0.5
!                  n          Pn(x)            Pn'(x)
!                ---------------------------------------
!                  0       1.00000000        .00000000
!                  1        .50000000       1.00000000
!                  2       -.12500000       1.50000000
!                  3       -.43750000        .37500000
!                  4       -.28906250      -1.56250000
!                  5        .08984375      -2.22656250
!       ========================================================

       SUBROUTINE LPN(N,X,PN)
!
!       ===============================================
!       Purpose: Compute Legendre polynomials Pn(x)
!                and their derivatives Pn'(x)
!       Input :  x --- Argument of Pn(x)
!                n --- Degree of Pn(x) ( n = 0,1,...)
!       Output:  PN(n) --- Pn(x)
!                PD(n) --- Pn'(x)
!       ===============================================
          IMPLICIT NONE
          REAL*8, INTENT(IN), DIMENSION(:) :: X
          INTEGER, INTENT(IN) :: N
          REAL*8, INTENT(OUT), DIMENSION(0:N,size(X)) :: PN

          REAL*8, DIMENSION(size(X)) :: P0,P1,PF
          INTEGER :: K
          REAL*8 :: XK

          PN(0,:)=1.0D0

          if(N > 0) then
             PN(1,:)=X
             P0(:)=1.0D0
             P1(:)=X
             DO 10 K=2,N
                XK=dble(K)
                PF=(2.0D0*XK-1.0D0)/XK*X*P1-(XK-1.0D0)/XK*P0
                PN(K,:)=PF
                P0=P1
                P1=PF
10           ENDDO
          endif
      END SUBROUTINE LPN


  END MODULE NNN_INTERACTION
