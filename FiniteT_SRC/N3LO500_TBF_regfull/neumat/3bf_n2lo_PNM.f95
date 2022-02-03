! ... NN interaction passed around
      MODULE THREE_N_INTERACTION
        USE precision_definition
        USE NN_interaction, ONLY : iz1,iz2,Jmin,Jmax,Nch
        USE mesh, ONLY : xk,Nkmesh

        REAL (long), ALLOCATABLE, DIMENSION(:,:,:,:) :: v3N
!        REAL(long), ALLOCATABLE, DIMENSION(:) :: gk,wgk,aux,ak
        real(long), ALLOCATABLE, DIMENSION(:) :: xk3,wk3,xnk,reg_int


        INTEGER(ilong) :: Nk
        REAL(long) :: kF,q,p
        REAL(long) :: gA,fPI,c1,c3,c4,cD,cE,mN,mPI,lchi,xlmax
        INTEGER(ilong) :: Nangle,Nk3,nexp,iregulator
        LOGICAL :: internalstep

      CONTAINS

!....cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... LECs & REGULATOR CHOICES FOR 3NF SUBROUTINE
!....cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SUBROUTINE LECs
          USE precision_definition
          USE physical_constants
          IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... Chiral constants
            gA=1.29_long
            fPI=92.4_long !MeV
            mN=xmass !MeV
            mPI=138.04_long !MeV
            lchi=700._long !MeV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... LECs
            c1=-0.81_long*1.E-3 !MeV^-1
            c3=-3.2_long*1.E-3 !MeV^-1

            c4= 0._long !MeV^-1
            cD= 0._long
            cE= 0._long

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! ... REGULATOR CUTOFF
            xlmax=500._long !MeV ! 2._long*
            ! ... INTERNAL REGULATOR ireg=2 OR EXTERNAL REGULATOR ireg=1
            iregulator=2
            ! ... nexp OF REGULATOR FUNCTION
            !           nexp=4
            nexp=2
            ! ... INTERNAL n(k) IS STEP FUNCTION (T) OR CORRELATED (F)
            internalstep=.true.

          END SUBROUTINE LECS
!....ccccccccccccccccccccccccccccccccccccccccccccccccc
!....Density dependent effective NN interaction from c
!....N2LO three-body forces in Nuclear Matter.       c
!....PRC 81, 024002 (2010)                           c
!....ccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE N3POT_PNM(rho)

        USE precision_definition
        USE physical_constants
        USE mesh_generator
        USE legendre_pols

        IMPLICIT NONE

        REAL(long), INTENT(IN) :: rho

        INTEGER(ilong) :: ik,jk,ian
        INTEGER(ilong) :: jj,ll,it,s,ich,itii
        INTEGER(ilong) :: Nangle,iterms

        REAL(long) :: v1,v2,v3,reg,rho_int
        REAL(long) :: G0,Gst,G2st,G1,G1st,G2,G3
        REAL(long) :: gamma0,gamma1,gamma2,gamma3
        REAL(long) :: xi,xf,pp2,qq2,ala

        REAL(long), ALLOCATABLE, DIMENSION(:,:) :: p2
        REAL(long), ALLOCATABLE, DIMENSION(:,:,:) :: q2

        REAL(long), ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: vNNN_aux
        REAL(long), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: vNNN,uNNN
        REAL(long), ALLOCATABLE, DIMENSION(:) :: angle, wangle

        REAL(long), ALLOCATABLE, DIMENSION(:,:) :: pleg
        REAL(long), ALLOCATABLE, DIMENSION(:) :: aux1,aux2,paux1,paux2,rgint

        write(*,*) 'Loading the averaged NNN potential in NM...'

        ! ... CONSTANTS FOR DIFFERENT TERMS (all potential terms are in fm)
        call LECs

        write(*,*) rho,mPi,c1,c3,mN,hbc
        write(*,*) iz1,iz2
! ... CONSTANTS FOR DIFFERENT TERMS (all potential terms are in fm)
        v1=hbc3*hbc*gA**2*mN/(8._long*pi*fPI**4)
        v2=hbc*gA**2*mN/(32._long*pi*pi2*fPI**4)
        v3=hbc*gA**2*mN/(64._long*pi*pi2*fPI**4)

! ... ANGULAR INTEGRATIONS
        Nk3=256
        ALLOCATE( xk3(Nk3),wk3(Nk3),xnk(Nk3),reg_int(Nk3),rgint(Nk3) )

        kf=(6._long*pi2/deg*rho)**(1d0/3d0)*hbc
        xi=0._long
        xf=kf
        call GAUSS(xi,xf,Nk3,xk3,wk3)
        xnk=1._long

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ... TEST DENSITY VALUE FOR VALIDITY OF xnk
100     format(' TEST RHO IN 3NF ROUTINE rho=',f8.6)
        write(*,100) deg/2._long/pi2*sum(xk3**2*xnk*wk3)/hbc3

        ! ... REGULATOR
        ! ... internal regulator function (x3lmax from modules)
        reg_int=1._long ! no internal regulator

! ... Allocating common matrices
        Nangle=24! 100
        ALLOCATE(vNNN(Nkmesh,Nkmesh,Nangle,5,2))
        ALLOCATE(angle(Nangle), wangle(Nangle))
        ALLOCATE(q2(Nkmesh,Nkmesh,Nangle))
        ALLOCATE(p2(Nkmesh,Nkmesh))
        xi=-1_long
        xf= 1_long
        call gauss(xi,xf,Nangle,angle,wangle)

! ... Allocating auxiliary matrices
        ALLOCATE(vNNN_aux(Nkmesh,Nkmesh,Nangle,5,2,6))

! ... Allocating auxiliary vectors
!        ALLOCATE(aux(Nk),aux1(Nk),aux2(Nk))

        vNNN_aux=0._long

        p2=0._long

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! ... symmetric extrapolation for off-shell elements
        FORALL(ik=1:Nkmesh,jk=1:Nkmesh)
           p2(ik,jk)=(xk(ik)**2+xk(jk)**2)/2._long
        END FORALL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... BUILD FUNCTIONS
        do ian=1,Nangle
           write(*,*) ian,Nangle
           do ik=1,Nkmesh
              do jk=1,Nkmesh

! ... momentum transfered
                 q2(ik,jk,ian)=xk(ik)**2 + xk(jk)**2 - &
                      2._long*xk(ik)*xk(jk)*angle(ian)
                 qq2=q2(ik,jk,ian)
                 q=sqrt(q2(ik,jk,ian))

! ... TOTAL MOMENTUM
                 p=sqrt(p2(ik,jk))
                 pp2=p2(ik,jk)
                
                 ! ... REGULATOR & DENSITY THAT APPEARS IN INTERACTION TERMS
                 if(iregulator==1) then
                    reg = EXP( -2._long*(pp2/xlmax**2)**nexp )
                    rho_int=rho
                 elseif(iregulator==2) then
                    reg=1._long
                    reg_int=1._long
                    rgint=-2._long*( (xk3**2/3._long + pp2 )/xlmax**2)**nexp
                    WHERE(abs(rgint) > 400d0 )
                      reg_int = 0.d0
                    ELSEWHERE
                      reg_int = exp( -2._long*( (xk3**2/3._long + pp2 )/xlmax**2)**nexp )
                    ENDWHERE
!                    reg_int =exp( -2._long*( (xk3**2/3._long + pp2 )/xlmax**2)**nexp )
                    rho_int=deg/2._long/pi2*sum(xk3**2*xnk*wk3*reg_int)/hbc3
                    if(ian==50.and.ik==25.and.jk==25) write(*,100) rho_int
                 endif
                 call GAMMA( p,pp2,q,qq2, &
                      gamma0,gamma1,gamma2,gamma3,G0,Gst,G2st,G1,G1st,G2,G3 )

                 ! ... TERMS OF INTERACTION
                 ! ... COMMENTS ARE CHANGES IN PNM
                 ! ... FIRST TERM
                 ! ... Term 1: 2 of rho + factor 2 in t1*t2 -> same as in SNM
                 vNNN_aux(ik,jk,ian,3,2,1) = reg*v1*rho_int* &
                      (2._long*c1*mPI**2 + c3*q2(ik,jk,ian)) &
                      /(mPI**2 + q2(ik,jk,ian))**2

                 ! ... SECOND TERM
                 ala=q2(ik,jk,ian)*(gamma0+2._long*gamma1+gamma3)+4._long*gamma2

                 ! Iso-singlet tensor
                 ! ... Term 2.1: 1/2 of t1*t2
                 vNNN_aux(ik,jk,ian,3,2,2) = reg*v2/(mPI**2+q2(ik,jk,ian))* &
                      ( -4._long*c1*mPI**2*(gamma0+gamma1) &
                        -c3*ala ) &
                      *(-real(iz1+iz2,long)/4d0)

                 ! Iso-vector tensor
                 ! ... Term 2.2: new term due to t1*t2
                 vNNN_aux(ik,jk,ian,3,1,2) = reg*v2/(mPI**2+q2(ik,jk,ian))* &
                      ( -4._long*c1*mPI**2*(gamma0+gamma1) &
                        -c3*ala ) &
                      *(0.5_long)

                 ! ... NOTE: THESE TWO ADD UP, SO ONE COULD TAKE (TERM 2.1) AND AVOID 1/2

                 ! Iso-singlet central
                 ! ... 1/3 in neutron matter
                 ala=  48._long/deg*pi2*rho_int*hbc3 &
                      -12._long*(2d0*mPI**2+q2(ik,jk,ian))*gamma0 &
                      -6._long*q2(ik,jk,ian)*gamma1 &
                      +3._long*(2d0*mPI**2+q2(ik,jk,ian))**2*G0

                 vNNN_aux(ik,jk,ian,1,1,3) = reg*v3*( &
                      -12d0*c1*mPI**2*(2d0*gamma0-(2d0*mPI**2+q2(ik,jk,ian))*G0) &
                      -c3*ala ) &
                      *(1d0/2d0+real(iz1+iz2,long)/12d0)


                 ! Iso-singlet SO
                   ala= 2._long*(gamma0+gamma1) &
                        -(2._long*mPI**2+q2(ik,jk,ian))*(G0+2d0*G1)

                 vNNN_aux(ik,jk,ian,4,1,3) = reg*v3*( &
                        -3d0*c3*( ala )  &
                        -12d0*c1*mPI**2*(G0 + 2d0*G1) ) &
                        *(1d0/2d0+real(iz1+iz2,long)/12d0)

              enddo ! Loop over second momentum
           enddo ! Loop over first momentum

        enddo ! Loop over angles

! ... choose density-dependent term or sum of all
        ! ... Sum of all density-dependent terms into two isospins
        do iterms=1,6
           vNNN(:,:,:,:,1)= vNNN(:,:,:,:,1) + vNNN_aux(:,:,:,:,1,iterms)
           vNNN(:,:,:,:,2)= vNNN(:,:,:,:,2) + vNNN_aux(:,:,:,:,2,iterms)
        enddo

        !        DEALLOCATE( aux,aux1,aux2,gk,wgk,vNNN_aux )
        DEALLOCATE( xk3,wk3,xnk,reg_int )
        DEALLOCATE( vNNN_aux )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ... PARTIAL WAVE EXPANSION
! ... allocate auxiliary variables
        ALLOCATE( aux1(Nangle), aux2(Nangle), paux1(Nangle), paux2(Nangle) )
        ALLOCATE( pleg(0:Jmax+1,Nangle) )

! ... allocate potential matrix
        ALLOCATE(uNNN(Nkmesh,Nkmesh,Nangle,5,0:1)) !fm
        if(.not. ALLOCATED(v3N)) ALLOCATE(v3N(Nkmesh,Nkmesh,0:Jmax,Nch))    !fm

      uNNN=0._long
      do it=0,1
         uNNN(:,:,:,:,it)=vNNN(:,:,:,:,1) &
              +real(4*it-3,long)*vNNN(:,:,:,:,2)
      enddo

      call LPN(Jmax+1,angle,pleg)

      v3N=0._long
! ... conditions on isospin
      if(iz1*iz2 == 1) then
         itii=1
      else
         itii=0
      endif

! ... momentum loop

      do ik = 1,Nkmesh
         do jk=1,Nkmesh

! ... filling channels
            do it=itii,1
               do ich=1,6
! ... singlet
                  if(ich == 1) then
                     s=0
                     do jj=Jmin,Jmax
                        ll=jj
                        if(mod(s+ll+it,2) == 0) cycle
                        aux1 = uNNN(ik,jk,:,1,it)-3._long*uNNN(ik,jk,:,2,it) &
                             -q2(ik,jk,:)*uNNN(ik,jk,:,3,it) &
                             +p2(ik,jk)**2*(angle**2-1._long)*uNNN(ik,jk,:,5,it)
                        paux1= pleg(jj,:)
                        v3N(ik,jk,jj,ich)= 0.5_long*sum( aux1*paux1*wangle )
                     enddo
! ... triplet
                  elseif(ich == 2) then
                     s=1
                     do jj=max(1,Jmin),Jmax
                        ll=jj
                        if(mod(s+ll+it,2) == 0) cycle
                        aux1 = 2._long*p2(ik,jk)*(uNNN(ik,jk,:,4,it) &
                             -uNNN(ik,jk,:,3,it)+p2(ik,jk)*angle* &
                             uNNN(ik,jk,:,5,it))
                        paux1 = pleg(jj+1,:)+pleg(jj-1,:)

                        aux2 =  uNNN(ik,jk,:,1,it)+uNNN(ik,jk,:,2,it) &
                             +2._long*p2(ik,jk)*(1._long+angle)*uNNN(ik,jk,:,3,it) &
                             -4._long*p2(ik,jk)*angle*uNNN(ik,jk,:,4,it) &
                             -p2(ik,jk)**2*(3._long*angle**2+1._long)*uNNN(ik,jk,:,5,it)
                        paux2 = pleg(jj,:)

                        v3N(ik,jk,jj,ich)= 0.5_long* &
                             sum( (aux1*paux1 + aux2*paux2)*wangle )
                     enddo
! ... coupled triplet V++
                  elseif(ich == 4)then
                     s=1
                     do jj=Jmin,Jmax
                        ll=jj+1
                        if(mod(s+ll+it,2) == 0) cycle
                        aux1 = 2._long*p2(ik,jk)*(uNNN(ik,jk,:,4,it) &
                             +1._long/real(2*jj+1,long) &
                             *(uNNN(ik,jk,:,3,it)-p2(ik,jk)* &
                             angle*uNNN(ik,jk,:,5,it)))
                        paux1 = pleg(jj,:)

                        aux2 = uNNN(ik,jk,:,1,it)+uNNN(ik,jk,:,2,it)+ &
                             p2(ik,jk)* &
                             ( p2(ik,jk)*(1-angle**2)*uNNN(ik,jk,:,5,it) &
                            - 2._long*angle*uNNN(ik,jk,:,4,it) &
                            + 2._long/real(2*jj+1,long) &
                             *(p2(ik,jk)*uNNN(ik,jk,:,5,it)-uNNN(ik,jk,:,3,it)))
                        paux2 = pleg(jj+1,:)

                        v3N(ik,jk,jj,ich)= 0.5_long* &
                             sum( (aux1*paux1 + aux2*paux2)*wangle )

                     enddo
! ... coupled triplet V--
                  elseif(ich == 3)then
                     s=1
                     do jj=max(1,Jmin),Jmax
                        ll=jj-1
                        if(mod(s+ll+it,2).eq.0) cycle
                        aux1 = 2._long*p2(ik,jk)*(uNNN(ik,jk,:,4,it) &
                             -1._long/real(2*jj+1,long) &
                             *(uNNN(ik,jk,:,3,it)-p2(ik,jk)* &
                             angle*uNNN(ik,jk,:,5,it)))
                        paux1 = pleg(jj,:)

                        aux2 = uNNN(ik,jk,:,1,it)+uNNN(ik,jk,:,2,it)+ &
                             p2(ik,jk)* &
                             ( p2(ik,jk)*(1-angle**2)*uNNN(ik,jk,:,5,it) &
                            - 2._long*angle*uNNN(ik,jk,:,4,it) &
                            - 2._long/real(2*jj+1,long) &
                             *(p2(ik,jk)*uNNN(ik,jk,:,5,it)-uNNN(ik,jk,:,3,it)))
                        paux2 = pleg(jj-1,:)

                        v3N(ik,jk,jj,ich)= 0.5_long* &
                             sum( (aux1*paux1 + aux2*paux2)*wangle )
                     enddo
! ... coupled triplet V-+
                  elseif(ich == 5)then
                     s=1
                     do jj=max(1,Jmin),Jmax
                        ll=jj-1
                        if(mod(s+ll+it,2) == 0)cycle
                        aux1 = uNNN(ik,jk,:,3,it)-p2(ik,jk)*uNNN(ik,jk,:,5,it)
                        paux1 = pleg(jj+1,:)

                        aux2= (real(2*jj,long)-angle*real(2*jj+1,long)) &
                             *uNNN(ik,jk,:,3,it)+p2(ik,jk)*angle &
                             *uNNN(ik,jk,:,5,it)
                        paux2=pleg(jj,:)

                        v3N(ik,jk,jj,ich)=(sqrt(real(jj+1,long)))*p2(ik,jk)/  &
                             (sqrt(real(jj,long))*real(2*jj+1,long))* &
                             sum( (aux1*paux1 + aux2*paux2)*wangle )
                     enddo
! ... coupled triplet V+-
                  elseif(ich == 6)then
                     s=1
                     do jj=max(1,Jmin),Jmax
                        ll=jj+1
                        if(mod(s+ll+it,2) == 0) cycle
                        !                        v3N(ik,jk,jj,ich)=v3N(jk,ik,jj,5)
                        v3N(jk,ik,jj,ich)=v3N(ik,jk,jj,5)
                     enddo
                  endif

               enddo ! ... loop filling channels
            enddo ! ... loop isospin
         enddo  ! ... loop second momentum
      enddo ! ... loop first momentum

      WHERE( abs(v3N) < 1e-30 ) v3N=0._long

! ... TO HAVE MATRIX ELEMENTS IN UNITS OF MEV-2, as 2NF FORCE
      xf=4_long*pi/xmass/(2_long*pi**2)/hbc !2/m/pi
      v3N=v3N*xf

      if( Jmin==0 ) then
         v3N(:,:,0,2)=v3N(:,:,0,4)
         v3N(:,:,0,4)=0._long
      endif

      write(*,*) 'DONE WITH 3B MATRIX ELEMENTS'

      DEALLOCATE( aux1,aux2,paux1,paux2,pleg )
      DEALLOCATE( uNNN,vNNN )
      DEALLOCATE( angle,wangle )

    END SUBROUTINE N3POT_PNM


    SUBROUTINE GAMMA( p,pp2,q,qq2,gamma0,gamma1,gamma2,gamma3,G0,Gst,G2st,G1,G1st,G2,G3 )

      USE physical_constants
      IMPLICIT NONE

      real(long), INTENT(IN) :: p,pp2,q,qq2
      real(long), INTENT(OUT) :: gamma0,gamma1,gamma2,gamma3
      real(long), INTENT(OUT) :: G0,Gst,G2st,G1,G1st,G2,G3
      real(long), ALLOCATABLE, DIMENSION(:) :: theta,aux,ak

      gamma0=0._long
      gamma1=0._long
      gamma2=0._long
      gamma3=0._long

      G0=0._long
      Gst=0._long
      G2st=0._long
      G1=0._long
      G1st=0._long
      G2=0._long
      G3=0._long

      ! ... Allocate auxiliary variable for gamma integrals
      ALLOCATE(theta(Nk3),aux(Nk3),ak(Nk3))

      ! ... define auxiliary variable for gamma integrals
      theta=xk3*xnk &
           *log( ((p+xk3)**2 + mPI**2) / ((p-xk3)**2+mPI**2) )


      ! ... Define gamma integrals
      gamma0=1._long/2._long/p*sum(wk3*reg_int*theta)

      gamma1=1._long/2._long/p**2*sum(wk3*reg_int*( &
           2._long*xk3**2*xnk - (p**2+xk3**2+mPI**2)/(2._long*p)*theta &
           ))

      gamma2=1._long/4._long/p*sum(wk3*reg_int*( &
           xk3**2*xnk*(p**2+xk3**2+mPI**2)/p &
           +(1._long-( (p**2+xk3**2+mPI**2)/(2._long*p*xk3))**2)*xk3**2*theta &
           ))

      gamma3=1._long/4._long/p**3*sum(wk3*reg_int*( &
           -3._long*xk3**2/p*(p**2+xk3**2+mPI**2)*xnk &
           +(3*((p**2+xk3**2+mPI**2)/(2._long*p*xk3))**2-1._long)*xk3**2*theta &
           ))


      ! ... AUXILIARY INTEGRANDS
      ak=(mPI**2+(xk3+p)**2)*(mPI**2+(xk3-p)**2)
      aux= xk3/sqrt( ak+qq2*xk3**2 )* &
           log( (q*xk3+sqrt(ak+qq2*xk3**2))/sqrt(ak) )

      G0=  2._long/q*sum(wk3*reg_int* xnk*aux)
      Gst= 2._long/q*sum(wk3*reg_int* xnk*aux*xk3**2)
      G2st=2._long/q*sum(wk3*reg_int* xnk*aux*xk3**4)

      ! ... USE RECURSION RELATIONS
      G1=(gamma0-(mPI**2+pp2)*G0-Gst)/(4._long*pp2-qq2)
      G1st=(3._long*gamma2+pp2*gamma3  &
           - (mPI**2+pp2)*Gst-G2st) / (4._long*pp2-qq2)
      G2=(mPI**2+pp2)*G1+Gst+G1st
      G3=(gamma1/2._long-2._long*(mPI**2+pp2)*G1-2._long*G1st-Gst) &
           /(4._long*pp2-qq2)

      DEALLOCATE(theta,aux,ak)

    END SUBROUTINE GAMMA

  END MODULE THREE_N_INTERACTION
