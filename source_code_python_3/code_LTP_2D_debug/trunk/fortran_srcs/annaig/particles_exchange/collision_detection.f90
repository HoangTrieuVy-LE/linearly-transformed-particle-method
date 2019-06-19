!!====================================================================
!!
!!          Particle-particle detection
!!
!!====================================================================

subroutine COLLISION_DETECTION(NPART, &
                                XMIN, XMAX, &
                                YMIN, YMAX, &
                                ZMIN, ZMAX, &
                               XPART, &
                               YPART, &
                               ZPART, &
                              NCLOSE, &
                          NPCLOSEMAX, &
                                HOME, &
                                NEAR  )

!!====================================================================
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE
use PARAM_PHYS

implicit none


!!====================================================================
! ARRAYS STATEMENT
!---------------------------------------------------------------------
integer,                        intent(in) :: NPART
real(kind=8),                   intent(in) :: XMIN
real(kind=8),                   intent(in) :: XMAX
real(kind=8),                   intent(in) :: YMIN
real(kind=8),                   intent(in) :: YMAX
real(kind=8),                   intent(in) :: ZMIN
real(kind=8),                   intent(in) :: ZMAX
real(kind=8), dimension(NPART), intent(in) :: XPART
real(kind=8), dimension(NPART), intent(in) :: YPART
real(kind=8), dimension(NPART), intent(in) :: ZPART


!!- Number of particle close
integer,                     intent(out) :: NCLOSE
integer,                     intent(in ) :: NPCLOSEMAX
!!- index of home particle
integer, dimension(NPCLOSEMAX), intent(out) :: HOME

!!- index of near particle
integer, dimension(NPCLOSEMAX), intent(out) :: NEAR

!!- Boundary condition
!! in -x --> BC(1) 
!! in +x --> BC(2) 
!! in -y --> BC(3) 
!! in +y --> BC(4) 
!! in -z --> BC(5) 
!! in +z --> BC(6)
!!
!! if BC(X) = 0 : wall (or cpu)
!!          = 1 : Periodic
!integer, dimension(6) :: BC

!!--------------------------------------------------------------------
!!- Size of each box
real(kind=8) :: DXBOX, DYBOX, DZBOX

integer :: IBOX, JBOX, KBOX

!!- number of box
integer :: NBOX
integer, parameter :: NXBOX = 42
integer, parameter :: NYBOX = 42
integer, parameter :: NZBOX = 42

integer, dimension(NXBOX,NYBOX,NZBOX) :: NPBOX
integer, dimension(:,:,:,:), allocatable :: PARTBOX

!- Index
integer :: NP, I, J, K, IJK, I0, J0, K0, NP0, NP1, NP2
integer :: IP1, IM1, JP1, JM1, KP1, KM1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

HOME(:) = 0
NEAR(:) = 0
NCLOSE = 0


NBOX = NXBOX*NYBOX

!NCLOSE = max(20,int(3.*NPART/NBOX))


!write(*,*)'PARTBOX allocation'
!write(*,*)'PARTBOX NBOX = ',NBOX
!write(*,*)'PARTBOX NPART = ',NPART
!write(*,*)'PARTBOX 1.1*NPART/NBOX = ',1.1*NPART/NBOX
!write(*,*)'PARTBOX NCLOSE = ',NCLOSE


!allocate(PARTBOX(NXBOX,NYBOX,NZBOX,NCLOSE))
!NCLOSE = 0


DXBOX = (XMAX - XMIN)/real(NXBOX)
DYBOX = (YMAX - YMIN)/real(NYBOX)


!write(900+MYID,*)XMIN,XMAX,DXBOX
write(900+MYID,*)YMIN,YMAX,DYBOX
write(900+MYID,*)ZMIN,ZMAX,DZBOX


NPBOX(:,:,:) = 0

!!--------------------------------------------------------------------
!!- Number of particle in each box
!!--------------------------------------------------------------------

do NP = 1, NPART

!!- index of the box containing the particle
 I0 = int((XPART(NP) - XMIN)/DXBOX) + 1
 J0 = int((YPART(NP) - YMIN)/DYBOX) + 1 


 NPBOX(I0,J0) = NPBOX(I0,J0) + 1

end do


NCLOSE = maxval(NPBOX(:,:,:))

allocate(PARTBOX(NXBOX,NYBOX,NCLOSE))


!!write(*,*)'Id=',MYID,' PARTBOX allocated with:',NCLOSE



NCLOSE = 0

NPBOX(:,:,:) = 0


do NP = 1, NPART

!!- index of the box containing the particle
 I0 = int((XPART(NP) - XMIN)/DXBOX) + 1
 J0 = int((YPART(NP) - YMIN)/DYBOX) + 1 

 NPBOX(I0,J0) = NPBOX(I0,J0) + 1

 PARTBOX(I0,J0,NPBOX(I0,J0)) = NP

end do



!!--------------------------------------------------------------------
!!- List of particles in each box
!!--------------------------------------------------------------------
NP0 = 1
do IJK = 1, NBOX

 J0 = int((IJK-1-NXBOX*NYBOX*(K0-1))/NXBOX) + 1
 I0 = IJK - NXBOX*(J0-1) - NXBOX*NYBOX*(K0-1)

!!--------------------------------------------------------------------
!!- Neighboring boxes
!!--------------------------------------------------------------------
!!- x-direction
 IP1 = I0 + 1
 IM1 = I0 - 1

 if(IP1 > NXBOX) IP1 = 1
 if(IM1 < 1    ) IM1 = NXBOX

!!- y-direction
 JP1 = J0 + 1
 JM1 = J0 - 1

 if(JP1 > NYBOX) JP1 = 1
 if(JM1 < 1    ) JM1 = NYBOX


 if(NPBOX(I0,J0) /= 0) then


!write(*,*)'Box:',I0,J0
!write(*,*)'contient:',NPBOX(I0,J0,K0),'part.'
!write(*,*)(PARTBOX(I0,J0,K0,NP),NP=1,NPBOX(I0,J0,K0))


 do NP1 = 1, NPBOX(I0,J0)

!!--------------------------------------------------------------------
!!- Box IM1, JM1, K0 - box 1
!!--------------------------------------------------------------------
  if(NPBOX(IM1,JM1) /= 0)  then
   do NP2 = 1, NPBOX(IM1,JM1)
    HOME(NP0) = PARTBOX(I0 ,J0 ,NP1)
    NEAR(NP0) = PARTBOX(IM1,JM1,NP2)
    NP0 = NP0 + 1
   end do
  end if

!!--------------------------------------------------------------------
!!- Box I0,JM1, K0 - box 2
!!--------------------------------------------------------------------
  if(NPBOX(I0,JM1) /= 0)  then
   do NP2 = 1, NPBOX(I0,JM1)
    HOME(NP0) = PARTBOX(I0 ,J0 ,NP1)
    NEAR(NP0) = PARTBOX(I0 ,JM1,NP2)
    NP0 = NP0 + 1
   end do
  end if
!!--------------------------------------------------------------------
!!- Box IP1, JM1, K0 - box 3
!!--------------------------------------------------------------------
  if(NPBOX(IP1,JM1) /= 0)  then 
   do NP2 = 1, NPBOX(IP1,JM1)
    HOME(NP0) = PARTBOX(I0 ,J0 ,NP1) 
    NEAR(NP0) = PARTBOX(IP1,JM1,NP2)
    NP0 = NP0 + 1
   end do
  end if

!!--------------------------------------------------------------------
!!- Box IM1, J0, K0 - box 4
!!--------------------------------------------------------------------
  if(NPBOX(IM1,J0) /= 0)  then 
   do NP2 = 1, NPBOX(IM1,J0)
    HOME(NP0) = PARTBOX(I0 ,J0 ,NP1) 
    NEAR(NP0) = PARTBOX(IM1,J0 ,NP2)
    NP0 = NP0 + 1
   end do
  end if

!!--------------------------------------------------------------------
!!- Box I0, J0, K0 - box 5
!!--------------------------------------------------------------------
  do NP2 = NP1+1, NPBOX(I0,J0)
   HOME(NP0) = PARTBOX(I0 ,J0 ,NP1) 
   NEAR(NP0) = PARTBOX(I0 ,J0 ,NP2)
   NP0 = NP0 + 1
  end do

!!--------------------------------------------------------------------
!!- Box IP1, J0, K0 - box 6
!!--------------------------------------------------------------------
  if(NPBOX(IP1,J0) /= 0)  then 
   do NP2 = 1, NPBOX(IP1,J0)
    HOME(NP0) = PARTBOX(I0 ,J0,NP1) 
    NEAR(NP0) = PARTBOX(IP1,J0,NP2)
    NP0 = NP0 + 1
   end do
  end if

!!--------------------------------------------------------------------
!!- Box IM1, JP1, K0 - box 7
!!--------------------------------------------------------------------
  if(NPBOX(IM1,JP1) /= 0)  then 
   do NP2 = 1, NPBOX(IM1,JP1)
    HOME(NP0) = PARTBOX(I0 ,J0 ,NP1) 
    NEAR(NP0) = PARTBOX(IM1,JP1,NP2)
    NP0 = NP0 + 1
   end do
  end if

!!--------------------------------------------------------------------
!!- Box I0,JP1, K0 - box 8 
!!--------------------------------------------------------------------
  if(NPBOX(I0,JP1,K0) /= 0)  then 
   do NP2 = 1, NPBOX(I0,JP1)
    HOME(NP0) = PARTBOX(I0 ,J0 ,NP1) 
    NEAR(NP0) = PARTBOX(I0 ,JP1,NP2)
    NP0 = NP0 + 1
   end do
  end if

!!--------------------------------------------------------------------
!!- Box IP1, JP1, K0 - box 9
!!--------------------------------------------------------------------
  if(NPBOX(IP1,JP1) /= 0)  then 
   do NP2 = 1, NPBOX(IP1,JP1)
    HOME(NP0) = PARTBOX(I0 ,J0 ,NP1) 
    NEAR(NP0) = PARTBOX(IP1,JP1,NP2)
    NP0 = NP0 + 1
   end do
  end if

 end do



!do NP=1,NP0-1
!write(*,*)NP,HOME(NP),NEAR(NP)
!end do
!pause

 end if !- if  NPBOX(I0,J0) /= 0

end do



NCLOSE = NP0 -1


!write(*,*)'NCLOSE=',NCLOSE
!do NP=1,NCLOSE
!!write(*,*)NP,HOME(NP),NEAR(NP)
!end do
!pause

deallocate(PARTBOX)




!!--------------------------------------------------------------------
10100 format (A,I2.2,A)
10101 format (A,I2.2,A,A)
10102 format (A,I2.2,A,I8)
10300 format (A,A)
10303 format (A,A,A)


end subroutine COLLISION_DETECTION
