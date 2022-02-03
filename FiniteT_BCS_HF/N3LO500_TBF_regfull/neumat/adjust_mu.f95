     ! ... THERMODYNAMICAL PROPERTIES OF THE SYSTEM
MODULE thermodynamical
  USE precision_definition
  REAL (long) :: rho,t,xmu
  
  INTEGER (ilong) :: nxkf,nt
  REAL (long) :: xxkf,xkfin,dxkf
  REAL (long) :: xkfermi,ef,spkf
  REAL (long) :: t_i,dt

  REAL (long), ALLOCATABLE, DIMENSION(:) :: spe
      
      CONTAINS
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... Computes the chemical potential by normalising the     c
! ... density                                                c
! ... THIS VERSION CORRECTS THE RANGE OF MU'S IF THE INITIAL c
! ... GUESS IS NOT GOOD ENOUGH                               c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE ADJUST_MU( gap )
        USE precision_definition
        USE physical_constants
        USE mesh
        USE interpolation
        USE fermibose
        USE mesh_generator
        
        IMPLICIT NONE
        REAL(long), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: gap
        REAL(long) :: xmurng,xmuin
        REAL(long), ALLOCATABLE, DIMENSION(:) :: xxmu,xden,xxx,rrr,xmui,xmuf,xxdh,xxdl
        INTEGER(ilong) :: Nmu

!        write(*,*) 'Finding the chemical potential...'

        Nmu=200
        xmurng=40d0

        ALLOCATE( xxmu(Nmu),xden(Nmu) )

! ... Range of MU's
! ... Lower limit is taken close to e(k_F)
        ALLOCATE( xmui(1),xxdh(1),xmuf(1),xxdl(1),xxx(1),rrr(1) )
        rrr=xxkf

        call LIN_INT(xkmesh,spe,rrr,xxx)
        spkf=xxx(1)

        xmui(1)=spkf-xmurng/2d0
        xmuf(1)=xmui(1)+xmurng

        go to 888
        call DENSITIES(xmui,xxdl,gap)
        call DENSITIES(xmuf,xxdh,gap)

! ... CHECK THAT THE INPUT RHO FALLS WITHIN THE REGION XMUI,XMUFI
      if(xxdl(1) .gt. rho) then
 886     continue 
         write(*,*) 'Initial den too high!'
         xmui(1)=xmui(1)-xmurng/2.d0
         xmuf(1)=xmui(1)+xmurng
         
         call DENSITIES(xmui,xxdl,gap)
         if(xxdl(1).gt.rho) go to 886

      elseif(xxdh(1) .lt. rho) then
 887     continue 
         write(*,*) 'Final den too low!'
         xmuf(1)=xmuf(1)+xmurng/2.d0
         xmui(1)=xmuf(1)-xmurng
         call DENSITIES(xmuf,xxdh,gap)
         if(xxdh(1).lt.rho) go to 887
      endif
888   continue
!      write(*,*) 'mu_in,mu_fi=',xmui,xmuf

! ... New mesh of chemical potentials
      call LINEAR(xmui(1),xmuf(1),Nmu,xxmu,xden)

      call DENSITIES(xxmu,xden,gap)

 900  format('Prev mu',2x,f18.9,/,'New  mu',2x,2(f18.9))
      xmuin=xmu
      rrr=rho

      call LIN_INT(xden,xxmu,rrr,xxx)
      xmu=xxx(1)
!      write(*,900) xmuin,xmu,rho
      xmui=xmu
      call DENSITIES(xmui,xxdl,gap)

!      write(*,*) xxdl
      
      if( isnan(xmu) ) then
         write(*,*) 'Program stopped NAN in mu'
         stop
      endif

      DEALLOCATE( rrr,xxx )
      
    END SUBROUTINE ADJUST_MU

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... Computes the densities from the spectral function c
! ... for a set of chemical potentials                  c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE DENSITIES(xxmu,xden,gap)
        USE precision_definition
        USE physical_constants
        USE mesh
        USE fermibose
        USE interpolation
        USE mesh_generator
        
        IMPLICIT NONE
        REAL(long), DIMENSION(:), INTENT(IN) :: xxmu
        REAL(long), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: gap
        REAL(long), DIMENSION(size(xxmu)), INTENT(OUT) :: xden

        REAL(long), ALLOCATABLE, DIMENSION(:) :: xi
        
        REAL(long) :: xsum,fac,xmmu
        INTEGER(ilong) :: Nmu,im
        REAL(long), ALLOCATABLE, DIMENSION(:,:) :: xint

        Nmu=size(xxmu)
        fac=deg*4._long*pi/(2._long*pi)**3!/hbc3

        ALLOCATE( xint(Nkmesh,Nmu) )
        ALLOCATE( xi(Nkmesh) )

        ! ... Loop over chemical potentials
        do 903 im=1,Nmu
           xmmu=xxmu(im)
           if( MAXVAL(abs(gap)) == 0._long ) then
              xint(:,im)=fermi( t,xmmu,spe )
           else
              xint(:,im)=momdis( t,xmmu,spe,gap )
           endif
903     enddo

! ... Integration of the momentum distribution
        do 905 im=1,Nmu
           xsum=sum( xint(:,im) * xkmesh**2 * wkmesh)
           xden(im)=fac*xsum
905     enddo        

        DEALLOCATE( xint,xi )

      END SUBROUTINE DENSITIES

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ... Computes the densities from the spectral function     c
! ... for a set of chemical potentials                      c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION MOMDIS(t,xmu,spe,gap)   
        USE precision_definition
        USE fermibose
        IMPLICIT NONE
        REAL(long) :: t,xmu
        REAL(long), INTENT(IN), ALLOCATABLE, DIMENSION(:) :: spe,gap
        REAL(long), DIMENSION(size(spe)) :: momdis
        REAL(long), DIMENSION(:), ALLOCATABLE :: xi

        ALLOCATE( xi((size(spe))) )

        if( MAXVAL(abs(gap)) == 0._long ) then
           momdis=fermi( t,xmu,spe )
        else
           xi=sqrt( (spe-xmu)**2 + gap**2 )
           momdis = (1._long - (spe-xmu)/xi*ttanh(t,xi) )/2._long
        endif
        
        DEALLOCATE(xi)
        
      END FUNCTION MOMDIS

    END MODULE thermodynamical
    
