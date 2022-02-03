! ... cccccccccccccccccccccccccccccccccccccccccccccc
! ... MODULE TO TRANSFORM TO MOMENTUM SPACE AN
! ... INTERACTION THAT IS IN OPERATOR FORM IN REAL SPACE
! ... cccccccccccccccccccccccccccccccccccccccccccccc
      MODULE KSPACE_MODULE
        USE precision_definition
        INTEGER(ilong) :: Nrmax,l1,l2,maxpot
        INTEGER :: lpot,ss,tt

      CONTAINS

! ... cccccccccccccccccccccccccccccccccccccccccccccc
! ... DETERMINES WHICH EXPECTATION VALUES TO
! ... COMPUTE
! ... cccccccccccccccccccccccccccccccccccccccccccccc
          SUBROUTINE AVFACNAME(itz1,itz2,jc,lsj,pwname,coupled)
            USE precision_definition
            USE physical_constants

            IMPLICIT NONE
            INTEGER(ilong), INTENT(IN) :: itz1,itz2,jc,lsj
            CHARACTER(5), INTENT(OUT) :: pwname
            LOGICAL, INTENT(OUT) :: coupled

            INTEGER(ilong) :: s,isum,ittot,icheck
            CHARACTER(9) :: JPW,LPW
            CHARACTER(1) :: SPW
            CHARACTER(2) :: Ja,La,Lb

            if( jc > 7 ) then
               write(*,*) 'J>7 is too large'
               stop
            endif

            LPW="SPDFGHIJ"
            JPW="01234567"
            SPW="1"

            s=2
! ... J=0
            if(jc == 0) then
               coupled=.false.
! ..... 1S0
               if(lsj == 1) then
                  l1=0
                  l2=0
                  s=0
! .... 1P0
               elseif(lsj == 2) then
                  l1=1
                  l2=1
                  s=1
               endif
            else
! ... Singlet
               if(lsj==1) then
                  coupled=.false.
                  l1=jc
                  l2=jc
                  s=0
! ... Uncoupled triplet
               elseif(lsj==2) then
                  coupled=.false.
                  l1=jc
                  l2=jc
                  s=1
! ... Coupled triplet - V--
               elseif(lsj==3) then
                  coupled=.true.
                  l1=jc-1
                  l2=jc-1
                  s=1
! ... Coupled triplet - V++
               elseif(lsj==4) then
                  l1=jc+1
                  l2=jc+1
                  s=1
! ... Coupled triplet - V-+
               elseif(lsj==5) then
                  coupled=.true.
                  l1=jc-1
                  l2=jc+1
                  s=1
! ... Coupled triplet - V+-
               elseif(lsj==6) then
                  coupled=.true.
                  l1=jc+1
                  l2=jc-1
                  s=1
               endif
            endif

! ... ANTISIMMETRY CONSTRAINTS FOR T=1
            ittot=itz1+itz2
            if( abs(ittot) == 2) then
               isum=1+l1+s
               if(mod(isum,2)==0) s=2
            endif

! ... FIND ISOSPIN
            if(jc==0) then
               tt=1
            else
! ... If the channel is the singlet or the coupled triplet,
! ... the isospin is given by J+1 modulus 2
               tt=mod( mod(s+l1,2)+1,2 )

! ... SAFE CHECK
               icheck=mod(l1+s+tt,2)
               if(icheck==0) then
                  write(*,*) 'error in pws?'
                  stop
               endif
            endif

            ss=s

            if( ss /= 2) then
                  write(SPW,'(SS,i1)') 2*s+1
                  if(jc < 8) then
                     Ja=JPW(jc+1:jc+1)
                  else
                     write(Ja,'(SS,i2)') jc
                  endif

                  if(.not.coupled) then
                     if(l1<8) then
                        La=LPW(l1+1:l1+1)
                     else
                        write(La,'(SS,i2)') l1
                     endif
                     pwname=SPW//trim(La)//trim(Ja)
                  elseif(coupled) then

                     if(l1<8) then
                        La=LPW(jc:jc)
                        Lb=LPW(jc+2:jc+2)
                     else
                        write(La,'(SS,i2)') jc-1
                        write(Lb,'(SS,i2)') jc+1
                     endif
                     pwname=SPW//trim(La)//trim(Lb)//trim(Ja)
                  endif
               endif

          END SUBROUTINE AVFACNAME


        END MODULE KSPACE_MODULE
