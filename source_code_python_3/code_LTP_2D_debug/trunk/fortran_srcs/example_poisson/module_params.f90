!module_params.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE  type_params
  IMPLICIT NONE

!! Parameter initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Number of points on each dimension
  INTEGER                                   :: ntx,nty

  !Local grid boundary coordinates (global indexes)
  INTEGER                                   :: sx, ex, sy, ey

  !Number of time step iteration
  INTEGER, PARAMETER                        :: it_max=100000

  integer, parameter                   :: dp = kind(1.d0)
  !F90 EPSILON function
  REAL(kind=dp), PARAMETER                   :: eps=EPSILON(1._8)

END MODULE  type_params
