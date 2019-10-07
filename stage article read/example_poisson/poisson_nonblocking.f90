!**********************************************************************
!   poisson_nonblocking.f90 - Resolution of the Poisson's equation using Jacobi
!   on the [0,1]x[0,1] domain by a finite difference method
!   with a Jacobi solver.
!   Delta u = f(x,y)= 2*(x*x-x+y*y -y)
!   u equal 0 on the boudaries
!   Exact solution : u= x*y*(x-1)*(y-1)
!
!   The u value is :
!    coef(1) = (0.5*hx*hx*hy*hy)/(hx*hx+hy*hy)
!    coef(2) = 1./(hx*hx)
!    coef(3) = 1./(hy*hy)
!
!    u(i,j)(n+1)= coef(1) * ( coef(2)*(u(i+1,j)+u(i-1,j)) &
!               + coef(3)*(u(i,j+1)+u(i,j-1)) - f(i,j))
!
!   In this version, we give the total interior number of points along x (ntx) 
!   and y (nty).
! 
!   hx is the step along x, hy along y
!    hx = 1./(ntx+1)
!    hy = 1./(nty+1)
!
!   On each proess :
!   1) decomposition of the domain
!   2) knowing his 4 neighbors
!   3) exchange the points on the interfaces (non_blocking version)
!   4) calculate
!   5) recomposition of the u matrix in a output file data.dat
!
!   Author          : Isabelle DUPAYS (CNRS/IDRIS - France)
!                     <Isabelle.Dupays@idris.fr>
!   Creation        : April 2012
!****************************************************************************

PROGRAM poisson
  USE TYPE_PARAMS
  USE PARALLEL
  USE CALCUL_POISSON

  IMPLICIT NONE

  !u and u_new solution at the n and n+1 iteration
  REAL(kind=dp),ALLOCATABLE,DIMENSION(:, :) :: u, u_new
  !Exact solution
  REAL(kind=dp),ALLOCATABLE,DIMENSION(:, :) :: u_exact
  !Number of iterations
  INTEGER                                  :: it
  !Convergence
  REAL(kind=dp)                             :: diffnorm, dwork
  !Time measurement
  REAL(kind=dp)                             :: t1, t2
  !Convergence test
  LOGICAL                                  :: convergence

  !****************************************************************************
  !MPI init
  CALL env_init

  !Creation of the 2D cartesian topology
  CALL topology_init

  !Index for each domain
  CALL domain_boundaries

  !Initialization of the f, u, u_new and u_exact
  CALL initialization(u, u_new, u_exact)

  !Neighbors 
  CALL domain_neighbours

  !Creation of the derived datatypes type_line and type_column
  CALL derived_datatypes

  !Time scheme
  it = 0
  convergence = .FALSE.

  !Time measure in the loop
  t1 = MPI_WTIME()

  DO WHILE ((.NOT. convergence) .AND. (it < it_max))

     it = it +1

     u(sx:ex,sy:ey) = u_new(sx:ex,sy:ey)

     !Exchange the points at the interface for u at the n iteration
     CALL non_blocking_communications(u)

     !Calculation  u at the n+1 iteration
     CALL computation(u,  u_new)

     !Calculation of the global error
     diffnorm = global_error (u, u_new)

     !Stop if we obtained the machine precision using F90 funtion EPSILON
     convergence = (diffnorm < eps)

     !Print for the process 0 diffnorm
     IF ((rank == 0) .AND. (MOD(it,100) == 0)) THEN
        PRINT *, 'Iteration ',it, ' global_error = ', diffnorm
     END IF

  END DO

  !Time measure at the end of the loop
  t2 = MPI_WTIME()

  IF (rank ==  0) THEN
     !Print convergence time for the process 0
     PRINT *, 'Convergence after ', it, ' iterations in  ', t2 - t1, ' seconds '

     !Compare the calculate and the exact solution on the process 0
     CALL output_results(u, u_exact)
  END IF
 
  !Write the results u(sx:ex,sy:ey) 
  !on each process
  CALL mpi_write(u)

  !MPI finalize
  CALL env_finalize

END PROGRAM poisson
