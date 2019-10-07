!!--------------------------------------------------------------------
!!
!!--------------------------------------------------------------------

module PARTICLE_PARALLEL

use DNS_DIM

implicit none

type PARTTYPE
 real(kind=8) :: XP
 real(kind=8) :: YP
 real(kind=8) :: ZP
 real(kind=8) :: UP
 real(kind=8) :: VP
 real(kind=8) :: WP


 real(kind=8) :: UP_NM1
 real(kind=8) :: VP_NM1
 real(kind=8) :: WP_NM1

 real(kind=8) :: UP_NM2
 real(kind=8) :: VP_NM2
 real(kind=8) :: WP_NM2

 real(kind=8) :: UFAP
 real(kind=8) :: VFAP
 real(kind=8) :: WFAP

 real(kind=8) :: UFAP_NM1
 real(kind=8) :: VFAP_NM1
 real(kind=8) :: WFAP_NM1

 real(kind=8) :: UFAP_NM2
 real(kind=8) :: VFAP_NM2
 real(kind=8) :: WFAP_NM2

!!- Stored particle velocity for Lagrangian correlation
 real(kind=8),dimension(NBLGRMAX) :: UPT0
 real(kind=8),dimension(NBLGRMAX) :: VPT0
 real(kind=8),dimension(NBLGRMAX) :: WPT0

!!- Stored fluid velocity for Lagrangian correlation
 real(kind=8),dimension(NBLGRMAX) :: UFAPT0
 real(kind=8),dimension(NBLGRMAX) :: VFAPT0
 real(kind=8),dimension(NBLGRMAX) :: WFAPT0

!!- Stored Scalar
 real(kind=8) :: TP
 real(kind=8) :: TP_NM1
 real(kind=8) :: TP_NM2
 real(kind=8) :: TFAP
 real(kind=8) :: TFAP_NM1
 real(kind=8) :: TFAP_NM2
 
!!- Stored particle temperature for Lagrangian correlation
 real(kind=8),dimension(NBLGRMAX) :: TPT0
 
!!- Stored fluid temperature for Lagrangian correlation
 real(kind=8),dimension(NBLGRMAX) :: TFAPT0

!!- Color of the particles
 real(kind=8) :: COLOR

!!====================================================
!!- Specific variable for investigating the subgrid
!!  turbulent scales
 real(kind=8) :: DUFAP
 real(kind=8) :: DVFAP
 real(kind=8) :: DWFAP
 real(kind=8),dimension(NBLGRMAX) :: DUFAPT0
 real(kind=8),dimension(NBLGRMAX) :: DVFAPT0
 real(kind=8),dimension(NBLGRMAX) :: DWFAPT0
!!====================================================


! integer :: PROC_ID
! integer :: ID 
end type PARTTYPE

!integer :: NDOUBLE = 27
!integer :: NDOUBLE = 35
!integer :: NDOUBLE = 36
!integer :: NDOUBLE = 42
integer :: NDOUBLE = 28 + NBLGRMAX*11


!integer :: NINTEGER = 2
integer :: NINTEGER = 0

!! MPI TYPE FOR PARTICLE EXCHANGE
integer :: MPI_PARTICLETYPE
integer :: TAG


!- Number of particle exchanged
integer, dimension(:), allocatable :: NBR_EXCHANGE

!- Particle structure
type(PARTTYPE), dimension(:,:), allocatable :: PART




contains

  subroutine MPI_PART_TYPE

    integer, dimension(0:1) ::  oldtypes, blockcounts, offsets(0:1)
!    integer(kind=MPI_OFFSET_KIND) :: extent
    integer :: extent

!   Setup description of the NDOUBLE MPI_DOUBLE_PRECISION fields in PARTTYPE  
    offsets(0) = 0 
    oldtypes(0) = MPI_DOUBLE_PRECISION
    blockcounts(0) = NDOUBLE

!  Setup description of the NINTEGER MPI_INTEGER fields in PARTYPE
!  Need to first figure offset by getting size of MPI_DOUBLE_PRECISION 
    call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent, ierr) 
    offsets(1) = NDOUBLE * extent 
    oldtypes(1) = MPI_INTEGER 
    blockcounts(1) = NINTEGER 

!  Now define structured type and commit it  
    call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, MPI_PARTICLETYPE, ierr) 
    call MPI_TYPE_COMMIT(MPI_PARTICLETYPE, ierr) 


  end subroutine MPI_PART_TYPE


end module PARTICLE_PARALLEL


