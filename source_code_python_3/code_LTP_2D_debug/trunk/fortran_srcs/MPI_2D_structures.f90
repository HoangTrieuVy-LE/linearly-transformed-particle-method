!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> MPI 2D structures 
!--------------------------------------------------------------------------- 


MODULE mpi_2d_structures_modf90
	USE mod_particle_2D_modf90
	USE calculsfor_var_modf90
	use calculsfor_ini_modf90
	USE data_launch
	USE Jacobi_method
	
	IMPLICIT NONE
	
	! DECALARATIONS
	
	! NB_NEIGHBOURS is already declared in pack
	
	INTEGER, DIMENSION(NB_NEIGHBOURS)                     :: NEIGHBOUR
	INTEGER, PARAMETER                                    :: FJ=1,BJ=2,FK=3,BK=4 ! neighbours : Forward_J, Backward_J, ... 
  	INTEGER, PARAMETER                                    :: FJBK=5,BJFK=6,BJBK=7,FJFK=8 ! Common edge neighbour
		
	! IND_leave indicates the identity(in ALL_PARTICLES) of particles will leave the present block
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: IND_leave
	
	! IND_recv_overlap indicates the identity of particles receving from neighbours
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: IND_recv_overlap
	
	! IND_inside indicates the identity of particles inside the present block
	INTEGER, DIMENSION(:), ALLOCATABLE                    :: IND_inside
	
	! IND_overlap indicates the identity of overlap particles in each direction
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: IND_overlap
	
	! COUNTER_inside indicates the number of particles inside the present block
	INTEGER                                               :: COUNTER_inside 
	
	! COUNTER_recv_overlap indicates the number of overlap particles receiving in each direction 

	INTEGER                                               :: COUNTER_recv_overlap 
	
	! COUNTER_leave indicates the number of particles leaving in each direction
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: COUNTER_leave
	
	! COUNTER_overlap indicates the number of overlap particles in each direction
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: COUNTER_overlap
	
	! ALL_PARTICLES stock all particle informations and should be broadcasted after being updated all particles displacement.
	TYPE(PARTTYPE), DIMENSION(:), ALLOCATABLE             :: ALL_PARTICLES
	
	! These following variable is the return of deformation matrix calculation
	DOUBLE PRECISION                                      :: Norm_inf_Dm1
    INTEGER(8)                                            :: indice_max_norm_Dm1
    INTEGER                                               :: Npart_block	
    
    ! These 5 following variables is about to calculate the velocity_field and deformation matrix in block loop
    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE        :: Dout
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)         :: Xpart_block
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)           :: Mblock
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)         :: Dblock
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)         :: velocity_field
	
	CONTAINS
	
		SUBROUTINE neighbouring
			NEIGHBOUR(FJ) = rank_up 
			NEIGHBOUR(BJ) = rank_down
			NEIGHBOUR(FK) = rank_right
			NEIGHBOUR(BK) = rank_left
			
			NEIGHBOUR(FJBK) = rank_up_left
			NEIGHBOUR(BJFK) = rank_down_right
			NEIGHBOUR(BJBK) = rank_down_left
			NEIGHBOUR(FJFK) = rank_up_right
			
	
		END SUBROUTINE neighbouring
		
		SUBROUTINE initiation_table
			ALLOCATE(IND_inside(number_of_particles))
			ALLOCATE(IND_leave(number_of_particles,NB_NEIGHBOURS))
			ALLOCATE(IND_overlap(number_of_particles,NB_NEIGHBOURS))
			
			ALLOCATE(IND_recv_overlap(number_of_particles,NB_NEIGHBOURS))
			

			ALLOCATE(COUNTER_overlap(NB_NEIGHBOURS,0:nb_proc-1))
			ALLOCATE(COUNTER_leave(NB_NEIGHBOURS,0:nb_proc-1))
			
			ALLOCATE(ALL_PARTICLES(number_of_particles))
			
			ALLOCATE(Xpart_block(2,number_of_particles))
			ALLOCATE(Mblock(number_of_particles))
			ALLOCATE(Dblock(4,number_of_particles))
			ALLOCATE(velocity_field(2,number_of_particles))
			ALLOCATE(Dout(4,number_of_particles))		
		END SUBROUTINE initiation_table
	
		
		SUBROUTINE parttype_convert
			DO i=1,number_of_particles
				ALL_PARTICLES(i)%Xp = Xparticle_read(1,i)
				ALL_PARTICLES(i)%Yp = Xparticle_read(2,i)
				ALL_PARTICLES(i)%Upx = velocity_read(1,i)
				ALL_PARTICLES(i)%Upy = velocity_read(2,i)
				ALL_PARTICLES(i)%Dp1 = D_read(1,i)
				ALL_PARTICLES(i)%Dp2 = D_read(2,i)
				ALL_PARTICLES(i)%Dp3 = D_read(3,i)
				ALL_PARTICLES(i)%Dp4 = D_read(4,i)
				ALL_PARTICLES(i)%ID = i
				ALL_PARTICLES(i)%mass = mass_read(i) 
			END DO
		
		END SUBROUTINE parttype_convert
		

		!----------------------------------------!
		! OVERLAP CRITERION FOR JUST ONE PARTICLE!
		!----------------------------------------!
			
		SUBROUTINE overlap_criterion(particle_k,overlap_check,pure_inside_check,pointu1,pointu2,pointu3,pointu4)
		! Check a particle if it satifies overlap_criterion
		! Retrieve the eigenvalues and eigenvectors to determine the direction movement of particles, 	then decide its status
		
		! Input: k-th_particle in "PARTYPE" type 
		! Output: logical
		
		
			IMPLICIT NONE
			TYPE(PARTTYPE), INTENT(in)                              :: particle_k
			LOGICAL, INTENT(out)                                    :: overlap_check
			LOGICAL, INTENT(out)                                    :: pure_inside_check
			LOGICAL                                                 :: partial_inside_check			
			! Initial idea
			INTEGER, PARAMETER                                      :: n=2
			DOUBLE PRECISION                                        :: a(n,n), x(n,n)
			DOUBLE PRECISION, PARAMETER                             :: abserr=1.0e-09
			INTEGER i, j


			double precision, dimension(n)             	            :: eigenvector_1
			double precision, dimension(n)                          :: eigenvector_2
			DOUBLE PRECISION                                        :: eigenvalue_1
			DOUBLE PRECISION                                        :: eigenvalue_2
			DOUBLE PRECISION, dimension(n), intent(out)             :: pointu1, pointu2,pointu3,pointu4
			DOUBLE PRECISION, dimension(n)                          :: pointu5,pointu6,pointu7,pointu8

			DOUBLE PRECISION		   :: D11,D12,D21,D22
			
			LOGICAL :: pointu1inside,pointu2inside,pointu3inside,pointu4inside, &
			 		  pointu5inside,pointu6inside,pointu7inside,pointu8inside
			
			D11 = particle_k%Dp1
			D12 = particle_k%Dp2
			D21 = particle_k%Dp3
			D22 = particle_k%Dp4

			! matrix D
  			a(1,1) = D11
  			a(1,2) = D12
  			a(2,1) = D21
  			a(2,2) = D22


			! print a header and the original matrix
			call Jacobi(a,x,abserr,n)

			eigenvector_1 = x(1,:) 
			eigenvector_2 = x(2,:)
			eigenvalue_1  = a(1,1)
			eigenvalue_2  = a(2,2)	
				
			pointu1(1) =  eigenvector_1(1) + eigenvector_2(1) + particle_k%Xp
			pointu1(2) =  eigenvector_1(2) + eigenvector_2(2) + particle_k%Yp
			pointu2(1) = -eigenvector_1(1)- eigenvector_2(1) + particle_k%Xp
			pointu2(2) = -eigenvector_1(2) - eigenvector_2(2) + particle_k%Yp
			pointu3(1) = -eigenvector_1(1)+ eigenvector_2(1) + particle_k%Xp
			pointu3(2) = -eigenvector_1(2) + eigenvector_2(2) + particle_k%Yp
			pointu4(1) =  eigenvector_1(1) - eigenvector_2(1) + particle_k%Xp
			pointu4(2) =  eigenvector_1(2) - eigenvector_2(2) + particle_k%Yp
			
			
			pointu5(1) =  eigenvector_1(1) + particle_k%Xp
			pointu5(2) =  eigenvector_1(2) + particle_k%Yp
			pointu6(1) = -eigenvector_1(1) + particle_k%Xp
			pointu6(2) = -eigenvector_1(2) + particle_k%Yp
			pointu7(1) =  eigenvector_2(1) + particle_k%Xp
			pointu7(2) =  eigenvector_2(2) + particle_k%Yp
			pointu8(1) = -eigenvector_2(1) + particle_k%Xp
			pointu8(2) = -eigenvector_2(2) + particle_k%Yp
			
			
!			print*,pointu1
!			print*,a
!			print*,'EIGENNN', eigenvalue_1,eigenvalue_2
!			print*,'EIGNENN', eigenvector_1,eigenvector_2

			pointu1inside = (pointu1(1)-start_x>0) &
				.and.(pointu1(1)-end_x<0)          &
				.and.(pointu1(2)-start_y>0)        &
				.and.(pointu2(2)-end_y<0)
			pointu2inside = (pointu2(1)-start_x>0) &
				.and.(pointu2(1)-end_x<0)          &
				.and.(pointu2(2)-start_y>0)        &
				.and.(pointu2(2)-end_y<0)
			pointu3inside = (pointu3(1)-start_x>0) &
				.and.(pointu3(1)-end_x<0)          &
				.and.(pointu3(2)-start_y>0)        &
				.and.(pointu3(2)-end_y<0)
			pointu4inside = (pointu4(1)-start_x>0) &
				.and.(pointu4(1)-end_x<0)          &
				.and.(pointu4(2)-start_y>0)        &
				.and.(pointu4(2)-end_y<0)
!			pointu5inside = (pointu5(1)-start_x>0) &
!				.and.(pointu5(1)-end_x<0)          &
!				.and.(pointu5(2)-start_y>0)        &
!				.and.(pointu5(2)-end_y<0)
!			pointu6inside = (pointu6(1)-start_x>0) &
!				.and.(pointu6(1)-end_x<0)          &
!				.and.(pointu6(2)-start_y>0)        &
!				.and.(pointu6(2)-end_y<0)
!			pointu7inside = (pointu7(1)-start_x>0) &
!				.and.(pointu7(1)-end_x<0)          &
!				.and.(pointu7(2)-start_y>0)        &
!				.and.(pointu7(2)-end_y<0)
!			pointu8inside = (pointu8(1)-start_x>0) &
!				.and.(pointu8(1)-end_x<0)          &
!				.and.(pointu8(2)-start_y>0)        &
!				.and.(pointu8(2)-end_y<0)
			
			! WE HAVE 3 CASE HERE, ONE IS PARTICLE IS ABSOLUTELY INSIDE ("PURE INSIDE") 
			! THE BLOCK HAVING THE RANK OF THE EXECUTING PROCESS OR ANOTHER BLOCK  
			! ANOTHER CASE IS BARYCENTER OF PARTICLE IS INSIDE THE BLOCK AND IT IS 
			! AN OVERLAP PARTICLE
			
			if(pointu1inside.and.pointu2inside.and.pointu3inside.and.pointu4inside)then 
				pure_inside_check = .true. ! This particle is inside block
			else 
				pure_inside_check = .false.	! This particle is maybe an overlap particle
			end if
			
			if ((particle_k%Xp>start_x .and. particle_k%Xp<end_x) &
		  .and. (particle_k%Yp>start_y .and. particle_k%Yp<end_y)) then
		  		partial_inside_check = .true.
		  	else
		  		partial_inside_check = .false.
		  	end if
			
			if ((partial_inside_check).and.(.not. pure_inside_check).and.(pointu1inside .or. &
										pointu2inside .or. &
										pointu3inside .or. &
										pointu4inside &
!										pointu5inside .or. &
!										pointu6inside .or. &
!										pointu7inside .or. &
!										pointu8inside &
										)) then
				overlap_check = .true.
			else
				overlap_check = .false.
			end if 

		END SUBROUTINE overlap_criterion
	
		
		
		SUBROUTINE particle_distribution
			IMPLICIT NONE
			logical                           :: local_overlapcheck, local_insidecheck
			DOUBLE PRECISION, dimension(2)    :: pointu1, pointu2,pointu3,pointu4

			integer :: ID
			
			COUNTER_inside            = 0
			COUNTER_leave(:,rank)     = 0
			COUNTER_overlap(:,rank)   = 0

	
			!for every particle coordinates in Xparticle, assign it rightly to the right domain
			DO i =1,number_of_particles
				CALL overlap_criterion(ALL_PARTICLES(i),local_overlapcheck,local_insidecheck,pointu1,pointu2,pointu3,pointu4)
				
				ALL_PARTICLES(i)%pointu1 = pointu1
				ALL_PARTICLES(i)%pointu2 = pointu2
				ALL_PARTICLES(i)%pointu3 = pointu3
				ALL_PARTICLES(i)%pointu4 = pointu4
				

				if (local_insidecheck) then
					COUNTER_inside = COUNTER_inside + 1
					IND_inside(COUNTER_inside) = i
				end if
				
				if (local_overlapcheck) then
					call particle_screening(ALL_PARTICLES(i))
				end if
			END DO
			
			! MPI_BCAST update COUNTER_overlap and COUNTER_leave
    		
			if(rank /= 0) then
      			call MPI_SEND(COUNTER_overlap(1,rank),NB_NEIGHBOURS,MPI_INTEGER,0,101,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1
       				call MPI_RECV(COUNTER_overlap(1,ID),NB_NEIGHBOURS,MPI_INTEGER,ID,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			end do
    		end if
    		call MPI_BCAST(COUNTER_overlap,NB_NEIGHBOURS*nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)
    		
    		if(rank /= 0) then
      			call MPI_SEND(COUNTER_leave(1,rank),NB_NEIGHBOURS,MPI_INTEGER,0,1001,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1
       				call MPI_RECV(COUNTER_leave(1,ID),NB_NEIGHBOURS,MPI_INTEGER,ID,1001,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			end do
    		end if
    		call MPI_BCAST(COUNTER_leave,NB_NEIGHBOURS*nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)
			
		END SUBROUTINE particle_distribution
		
		
		SUBROUTINE particle_screening(particle_k)
			! For an overlap table of a block, we will precise in which neighbour block, these overlap
			! particles access.
			
			! Idea, we make a loop in overlap table, and condition for an overlap particle access in
			! a neighbour block is JUST ONE of 4 pointu locate inside of the neigbour.
			type(PARTTYPE), INTENT(in)     :: particle_k
			double precision, dimension(2) :: p1,p2,p3,p4
		
				p1 = particle_k%pointu1
				p2 = particle_k%pointu2
				p3 = particle_k%pointu3
				p4 = particle_k%pointu4
				
				
				if (rank_up /= -1) then
					if(((p1(2) > start_y+block_step_y .and. p1(2) < end_y+block_step_y) & 
					.and. (p1(1) > start_x .and. p1(1) < end_x)) .or. &
					((p2(2) > start_y+block_step_y .and. p2(2) < end_y+block_step_y) & 
					.and. (p2(1) > start_x .and. p2(1) < end_x)) .or. & 
					((p3(2) > start_y+block_step_y .and. p3(2) < end_y+block_step_y) & 
					.and. (p3(1) > start_x .and. p3(1) < end_x)) .or. &
					((p4(2) > start_y+block_step_y .and. p4(2) < end_y+block_step_y) & 
					.and. (p4(1) > start_x .and. p4(1) < end_x)))&
					then
						COUNTER_overlap(FJ,rank) = COUNTER_overlap(FJ,rank) + 1
						IND_overlap(COUNTER_overlap(FJ,rank),FJ) =  particle_k%ID
						
					end if
				end if
			
				if (rank_down /= -1) then
					if (( (p1(2) > start_y-block_step_y .and. p1(2) < end_y-block_step_y) &
					.and. (p1(1) > start_x .and. p1(1) < end_x) ) .or.&
					 ( (p2(2) > start_y-block_step_y .and. p2(2) < end_y-block_step_y) &
					.and. (p2(1) > start_x .and. p2(1) < end_x) ) .or.&
					 ( (p3(2) > start_y-block_step_y .and. p3(2) < end_y-block_step_y) &
					.and. (p3(1) > start_x .and. p3(1) < end_x) ) .or.&
					 ( (p4(2) > start_y-block_step_y .and. p4(2) < end_y-block_step_y) &
					.and. (p4(1) > start_x .and. p4(1) < end_x) )) &
					then
						COUNTER_overlap(BJ,rank) = COUNTER_overlap(BJ,rank) + 1
						IND_overlap(COUNTER_overlap(BJ,rank),BJ) =  particle_k%ID
					end if
					
				end if
				if (rank_left /= -1) then
					if (( (p1(2) > start_y .and. p1(2) < end_y) .and. (p1(1) > start_x-block_step_x .and. p1(1) < end_x-block_step_x) ) .or. &
					 ( (p2(2) > start_y .and. p2(2) < end_y) .and. (p2(1) > start_x-block_step_x .and. p2(1) < end_x-block_step_x) ) .or. &
					 ( (p3(2) > start_y .and. p3(2) < end_y) .and. (p3(1) > start_x-block_step_x .and. p3(1) < end_x-block_step_x) ) .or. &
					 ( (p4(2) > start_y .and. p4(2) < end_y) .and. (p4(1) > start_x-block_step_x .and. p4(1) < end_x-block_step_x) )) &
					 
					then
						COUNTER_overlap(BK,rank) = COUNTER_overlap(BK,rank) + 1
						IND_overlap(COUNTER_overlap(BK,rank),BK) =  particle_k%ID
!						if (rank==2) then 
!							print*,COUNTER_overlap(BK,rank)
!						end if
						
					end if
				
				end if
				if (rank_right /= -1) then
					if (( (p1(2) > start_y .and. p1(2) < end_y) & 
					.and. (p1(1) > start_x+block_step_x .and. p1(1) < end_x+block_step_x) ).or. &
					( (p2(2) > start_y .and. p2(2) < end_y) & 
					.and. (p2(1) > start_x+block_step_x .and. p2(1) < end_x+block_step_x) ).or. &
					( (p3(2) > start_y .and. p3(2) < end_y) & 
					.and. (p3(1) > start_x+block_step_x .and. p3(1) < end_x+block_step_x) ).or. &
					( (p4(2) > start_y .and. p4(2) < end_y) & 
					.and. (p4(1) > start_x+block_step_x .and. p4(1) < end_x+block_step_x) )) &
					
					 then
					 	COUNTER_overlap(FK,rank) = COUNTER_overlap(FK,rank) + 1
						IND_overlap(COUNTER_overlap(FK,rank),FK) =  particle_k%ID
					end if
			
				end if
				if (rank_up_left /= -1) then
					if (( (p1(2) > start_y+block_step_y .and. p1(2) < end_y+block_step_y) &
			       .and. (p1(1) > start_x-block_step_x .and. p1(1) < end_x-block_step_x) ).or.&
			       	( (p2(2) > start_y+block_step_y .and. p2(2) < end_y+block_step_y) &
			       .and. (p2(1) > start_x-block_step_x .and. p2(1) < end_x-block_step_x) ).or.&
			       ( (p3(2) > start_y+block_step_y .and. p3(2) < end_y+block_step_y) &
			       .and. (p3(1) > start_x-block_step_x .and. p3(1) < end_x-block_step_x) ).or.&
			       ( (p4(2) > start_y+block_step_y .and. p4(2) < end_y+block_step_y) &
			       .and. (p4(1) > start_x-block_step_x .and. p4(1) < end_x-block_step_x) ))&
			       
			        then
			        	COUNTER_overlap(FJBK,rank) = COUNTER_overlap(FJBK,rank) + 1
						IND_overlap(COUNTER_overlap(FJBK,rank),FJBK) =  particle_k%ID
					end if
					
				
				end if
				if (rank_up_right /= -1) then
					if (( (p1(2) > start_y+block_step_y .and. p1(2) < end_y+block_step_y) &
			       .and. (p1(1) > start_x+block_step_x .and. p1(1) < end_x+block_step_x) ).or. &
			       ( (p2(2) > start_y+block_step_y .and. p2(2) < end_y+block_step_y) &
			       .and. (p2(1) > start_x+block_step_x .and. p2(1) < end_x+block_step_x) ).or. &
			       ( (p3(2) > start_y+block_step_y .and. p3(2) < end_y+block_step_y) &
			       .and. (p3(1) > start_x+block_step_x .and. p3(1) < end_x+block_step_x) ).or. &
			       ( (p4(2) > start_y+block_step_y .and. p4(2) < end_y+block_step_y) &
			       .and. (p4(1) > start_x+block_step_x .and. p4(1) < end_x+block_step_x) ))&
			        then
			        	COUNTER_overlap(FJFK,rank) = COUNTER_overlap(FJFK,rank) + 1
						IND_overlap(COUNTER_overlap(FJFK,rank),FJFK) =  particle_k%ID
			        end if
				
				end if
				if (rank_down_left /= -1) then 
					if (( (p1(2) > start_y-block_step_y .and. p1(2) < end_y-block_step_y) &
			       .and. (p1(1) > start_x-block_step_x .and. p1(1) < end_x-block_step_x) ).or. &
			       ( (p2(2) > start_y-block_step_y .and. p2(2) < end_y-block_step_y) &
			       .and. (p2(1) > start_x-block_step_x .and. p2(1) < end_x-block_step_x) ).or. &
			       ( (p3(2) > start_y-block_step_y .and. p3(2) < end_y-block_step_y) &
			       .and. (p3(1) > start_x-block_step_x .and. p3(1) < end_x-block_step_x) ).or. &
			       ( (p4(2) > start_y-block_step_y .and. p4(2) < end_y-block_step_y) &
			       .and. (p4(1) > start_x-block_step_x .and. p4(1) < end_x-block_step_x) )) &
			        then
						COUNTER_overlap(BJBK,rank) = COUNTER_overlap(BJBK,rank) + 1
						IND_overlap(COUNTER_overlap(BJBK,rank),BJBK) =  particle_k%ID
					end if
				end if
				if (rank_down_right /= -1) then
					if ( ((p1(2) > start_y-block_step_y .and. p1(2) < end_y-block_step_y) &
			       .and. (p1(1) > start_x+block_step_x .and. p1(1) < end_x+block_step_x) ).or.&
			       ( (p2(2) > start_y-block_step_y .and. p2(2) < end_y-block_step_y) &
			       .and. (p2(1) > start_x+block_step_x .and. p2(1) < end_x+block_step_x) ).or.&
			       ( (p3(2) > start_y-block_step_y .and. p3(2) < end_y-block_step_y) &
			       .and. (p3(1) > start_x+block_step_x .and. p3(1) < end_x+block_step_x) ).or.&
			       ( (p4(2) > start_y-block_step_y .and. p4(2) < end_y-block_step_y) &
			       .and. (p4(1) > start_x+block_step_x .and. p4(1) < end_x+block_step_x) ))&
			       
			        then 
						COUNTER_overlap(BJFK,rank) = COUNTER_overlap(BJFK,rank) + 1
						IND_overlap(COUNTER_overlap(BJFK,rank),BJFK) =  particle_k%ID
					end if
				end if
			
		END SUBROUTINE
		
		
		SUBROUTINE send_overlap_particle
		! RECEIVE overlap particles from all possible neighbour blocks
		! After particle_screening, we will know in every blocks, which overlap particle
		! will overlap which block	
		
		! Idea, we concatenate all possible overlap particles from all possible neighbour
		! blocks into IND_overlap
		
		! Firstly, send its informations about "overlap_direction" to all neighbour blocks
		! Then receive informations from its neighbour blocks.
		
		! MPI.SEND
		! MPI.RECV
		integer,dimension(MPI_STATUS_SIZE) :: status
		integer :: i,neighloop
		integer, dimension(NB_NEIGHBOURS) :: OPP
		
		COUNTER_recv_overlap = 0
		
		OPP = (/2,1,4,3,6,5,8,7/)		

			DO neighloop = 1,8

				DO i=1,COUNTER_overlap(neighloop,rank)
					if (NEIGHBOUR(neighloop)<0) then
						cycle
					end if
					call MPI_SEND(IND_overlap(i,neighloop),1,MPI_INTEGER,NEIGHBOUR(neighloop),101,MPI_COMM_WORLD,code)
					
				END DO
			
				DO i= 1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
					if (NEIGHBOUR(neighloop)<0) then
						cycle
					end if
					call MPI_RECV(IND_recv_overlap(i,NEIGHBOUR(neighloop)),1,MPI_INTEGER,NEIGHBOUR(neighloop),101, & 
					MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
				END DO
				if (NEIGHBOUR(neighloop)<0) then
						cycle
				else	
					COUNTER_recv_overlap = COUNTER_recv_overlap + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
				end if
				
			END DO	
  			call MPI_BARRIER(MPI_COMM_WORLD,code)
			
		END SUBROUTINE send_overlap_particle
		
		
		
		SUBROUTINE block_loop_on_block_global_table	
		IMPLICIT NONE
		integer, dimension(NB_NEIGHBOURS) :: OPP
		integer :: i,j,neighloop,Ncum1,Ncum2,Ncum3,Ncum4,Ncum5
		
		Ncum1 = 0
		Ncum2 = 0
		Ncum3 = 0
		Ncum4 = 0
		Ncum5 = 0
		
		OPP  = (/2,1,4,3,6,5,8,7/)	
		! UPDATE GLOBAL_TABLE
		! This part's mission is about to include not only IND_inside and IND_overlap
		! but also all possible overlap particles from all possible neighbours 

		! We will have "Npart_block" became: COUNTER_inside + inside sum of COUNTER_overlap in 
		! all direction + COUNTER_recv
		! Where COUNTER_recv = sum of COUNTER_overlap from neigbour !TODO COMPLICATED
			
!			print*,'rank',rank,'COUNTER_recv_overlap',COUNTER_recv_overlap
			
			Npart_block = COUNTER_recv_overlap + COUNTER_inside + SUM(COUNTER_overlap(:,rank))
!			print*,'rank',rank,'Npart_block = ',Npart_block	

		DO i= 1,COUNTER_inside
			Xpart_block(1,i) = ALL_PARTICLES(IND_inside(i))%Xp
			Xpart_block(2,i) = ALL_PARTICLES(IND_inside(i))%Yp
		END DO
		Ncum1 = Ncum1 + COUNTER_inside
		
		DO neighloop=1,8
			if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
			DO j=1,COUNTER_overlap(neighloop,rank)
				Xpart_block(1,j+Ncum1) = ALL_PARTICLES(IND_overlap(j,neighloop))%Xp
				Xpart_block(2,j+Ncum1) = ALL_PARTICLES(IND_overlap(j,neighloop))%Yp
			END DO
			Ncum1 = Ncum1 + COUNTER_overlap(neighloop,rank)
!			if (rank==0) then
!				print*, Ncum
!			end if
			DO j=1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
				Xpart_block(1,j+Ncum1) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Xp
				Xpart_block(1,j+Ncum1) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Yp
			END DO
			Ncum1 = Ncum1 + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
		
		END DO
		
		DO i= 1,COUNTER_inside
			Mblock(i) = ALL_PARTICLES(IND_inside(i))%mass
		END DO
		Ncum2 = Ncum2 + COUNTER_inside
		
		DO neighloop=1,8
			if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
			DO j=1,COUNTER_overlap(neighloop,rank)
				Mblock(j+Ncum2) = ALL_PARTICLES(IND_overlap(j,neighloop))%mass
			END DO
			Ncum2 = Ncum2 + COUNTER_overlap(neighloop,rank)
			
			DO j=1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
				Mblock(j+Ncum2) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%mass
			END DO
			Ncum2 = Ncum2 + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
		END DO
	
!		
		DO i= 1,COUNTER_inside
			Dblock(1,i) = ALL_PARTICLES(IND_inside(i))%Dp1
			Dblock(2,i) = ALL_PARTICLES(IND_inside(i))%Dp2
			Dblock(3,i) = ALL_PARTICLES(IND_inside(i))%Dp3
			Dblock(4,i) = ALL_PARTICLES(IND_inside(i))%Dp4
		END DO
		
		Ncum3 = Ncum3 + COUNTER_inside
!		
		DO neighloop=1,8
			if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
			DO j=1,COUNTER_overlap(neighloop,rank)
				Dblock(1,j+Ncum3) = ALL_PARTICLES(IND_overlap(j,neighloop))%Dp1
				Dblock(2,j+Ncum3) = ALL_PARTICLES(IND_overlap(j,neighloop))%Dp2
				Dblock(3,j+Ncum3) = ALL_PARTICLES(IND_overlap(j,neighloop))%Dp3
				Dblock(4,j+Ncum3) = ALL_PARTICLES(IND_overlap(j,neighloop))%Dp4
			END DO
			Ncum3 = Ncum3 + COUNTER_overlap(neighloop,rank)
			
			DO j=1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
				Dblock(1,j+Ncum3) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Dp1
				Dblock(2,j+Ncum3) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Dp2
				Dblock(3,j+Ncum3) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Dp3
				Dblock(4,j+Ncum3) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Dp4
			END DO
			Ncum3 = Ncum3 + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
		END DO
		
!		
		DO i= 1,COUNTER_inside
			velocity_field(1,i) = ALL_PARTICLES(IND_inside(i))%Upx
			velocity_field(2,i) = ALL_PARTICLES(IND_inside(i))%Upy
		END DO
		Ncum4 = Ncum4 + COUNTER_inside
		DO neighloop=1,8
			if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
			DO j=1,COUNTER_overlap(neighloop,rank)
				velocity_field(1,j+Ncum4) = ALL_PARTICLES(IND_overlap(j,neighloop))%Upx
				velocity_field(2,j+Ncum4) = ALL_PARTICLES(IND_overlap(j,neighloop))%Upy
			END DO
			
			Ncum4 = Ncum4 + COUNTER_overlap(neighloop,rank)
			
			DO j=1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
				velocity_field(1,j+Ncum4) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Upx
				velocity_field(2,j+Ncum4) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Upy
			END DO
			Ncum4 = Ncum4 + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
		END DO
!		
		
		if(rank==0) then
!		print*, Dblock(2,1:1500)
!		print*,'Mblock(Npart_block)',Xpart_block(1,1:Npart_block+1)
!		print*, 'rank',rank,'velocity',velocity_field(1,1:10)
!		print*,'Before',velocity_field(:,1:5)
!		print*,'Before',Xpart_block(:,1:10)
!		print*, Mblock(1:10)
!		print*,'hx', hx, 'dt', dt	
		end if
	
		call diffusion_field_ltp(Xpart_block,Mblock,Dblock,hx,hy,velocity_field,Npart_block)
!		call diffusion_field_ltp(Xparticle_read,mass_read,D_read,hx,hy,velocity_read,number_of_particles)
		


		
!		call update_d_diffusion(Xpart_block, Dblock, Mblock, Tini, dt, time_scheme, hx, hy, &
!     indice_max_norm_Dm1, Norm_inf_Dm1, Dout, Npart_block)
!     	call update_d_diffusion(Xparticle_read, D_read, mass_read, Tini, dt, time_scheme, hx, hy, &
!      indice_max_norm_Dm1, Norm_inf_Dm1, Dout, number_of_particles)
     
		DO i=1,Npart_block
			Xpart_block(1,i) = Xpart_block(1,i) - velocity_field(1,i)*dt
			Xpart_block(2,i) = Xpart_block(2,i) - velocity_field(2,i)*dt
		END DO
		
		if(rank==0) then
!		print*, 'After', Xpart_block(1,100:110)
!     	print*,'rank',rank,'After', Dout(1,1:Npart_block)
!		print*,'After',velocity_field(:,1:5)
     	end if

	
		END SUBROUTINE block_loop_on_block_global_table	
	
	
	
	
		SUBROUTINE update_displacement_on_just_IND_inside_and_IND_overlap
		integer :: neighloop, Ncum1, Ncum2, Ncum3, Ncum4, Ncum5
		integer :: i,j,ID
		integer,dimension(MPI_STATUS_SIZE) :: status_msg
		integer,dimension(1:nb_proc-1)   ::address
		address(:) = 0
		Ncum1 = 0
		Ncum2 = 0
		Ncum3 = 0
		Ncum4 = 0
		Ncum5 = 0

		call MPI_PART_TYPE
		
!		if (rank==0) then

!		print*,rank,'Before', IND_inside(1:COUNTER_inside)
!		print*, 'Before',ALL_PARTICLES(IND_overlap(1:COUNTER_overlap(1,rank),1))%Xp
			DO i=1,COUNTER_inside
				ALL_PARTICLES(IND_inside(i))%Xp = Xpart_block(1,i)
				ALL_PARTICLES(IND_inside(i))%Yp = Xpart_block(2,i)
				if (rank /=0) then
				address(rank) = IND_inside(i)
				call MPI_SEND(address(rank),1,MPI_INTEGER,0,10,MPI_COMM_WORLD,code)
				call MPI_SEND(ALL_PARTICLES(address(rank)),1,MPI_PARTICLETYPE,0,9,MPI_COMM_WORLD,code)
				else
					DO ID =1,nb_proc-1
					call MPI_RECV(address(ID),1,MPI_INTEGER,ID,10,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
					call MPI_RECV(ALL_PARTICLES(address(ID)),1,MPI_PARTICLETYPE,ID,9,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
					END DO
				end if
			END DO

			Ncum1 = Ncum1 + COUNTER_inside
			
			address(:) = 0
			
			DO neighloop=1,8
				if (NEIGHBOUR(neighloop)<0) then
					cycle
				end if
				DO j=1,COUNTER_overlap(neighloop,rank)
					ALL_PARTICLES(IND_overlap(j,neighloop))%Xp = Xpart_block(1,j+Ncum1)
					ALL_PARTICLES(IND_overlap(j,neighloop))%Yp = Xpart_block(2,j+Ncum1)
							
				if (rank /=0) then
				address(rank) = IND_overlap(j,neighloop)
				call MPI_SEND(address(rank),1,MPI_INTEGER,0,110,MPI_COMM_WORLD,code)
				call MPI_SEND(ALL_PARTICLES(address(rank)),1,MPI_PARTICLETYPE,0,99,MPI_COMM_WORLD,code)
				else
					DO ID =1,nb_proc-1
					call MPI_RECV(address(ID),1,MPI_INTEGER,ID,110,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
					call MPI_RECV(ALL_PARTICLES(address(ID)),1,MPI_PARTICLETYPE,ID,99,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
					END DO
				end if
				END DO
				Ncum1 = Ncum1 + COUNTER_overlap(neighloop,rank)
			END DO
		
			
		END SUBROUTINE update_displacement_on_just_IND_inside_and_IND_overlap
		
!		SUBROUTINE update_deformation_matrix_on_just_IND_inside_and_IND_overlap
!			DO i=1,COUNTER_inside
!				ALL_PARTICLES(IND_inside(i))%Dp1 = Dout(1,i) 
!				ALL_PARTICLES(IND_inside(i))%Dp2 = Dout(2,i)
!				ALL_PARTICLES(IND_inside(i))%Dp3 = Dout(3,i)
!				ALL_PARTICLES(IND_inside(i))%Dp4 = Dout(4,i) 
!			END DO
!			DO i=1,COUNTER_overlap
!				ALL_PARTICLES(IND_overlap(i))%Dp1 = Dout(1,i+COUNTER_inside) 
!				ALL_PARTICLES(IND_overlap(i))%Dp2 = Dout(2,i+COUNTER_inside)
!				ALL_PARTICLES(IND_overlap(i))%Dp3 = Dout(3,i+COUNTER_inside)
!				ALL_PARTICLES(IND_overlap(i))%Dp4 = Dout(4,i+COUNTER_inside)
!			END DO
!		END SUBROUTINE update_deformation_matrix_on_just_IND_inside_and_IND_overlap
!		
!		SUBROUTINE update_table_block
!			IMPLICIT NONE
!			LOGICAL :: overlap_check, pure_inside_check
!			double precision, dimension(2) :: pointu1,pointu2,pointu3,pointu4

!			DO i=1,COUNTER_inside
!				call overlap_criterion(ALL_PARTICLES(IND_inside(i)),overlap_check,pure_inside_check,pointu1,pointu2,pointu3,pointu4)
!							
!				if (pure_inside_check) then
!					cycle
!					
!				else
!					if(overlap_check) then
!						COUNTER_overlap = COUNTER_overlap + 1
!						
!						IND_inside(i:(COUNTER_inside-1)) = IND_inside((i+1):COUNTER_inside)
!						COUNTER_inside = COUNTER_inside - 1
!					else
!						IND_leave(COUNTER_leave+1) = ALL_PARTICLES(IND_inside(i))%ID
!						COUNTER_leave = COUNTER_leave + 1
!					end if
!				end if
!			END DO		
!			
!			DO i=1, COUNTER_overlap
!				call overlap_criterion(ALL_PARTICLES(IND_overlap(i)),overlap_check,pure_inside_check,pointu1,pointu2,pointu3,pointu4)
!				
!				if (overlap_check) then
!					cycle
!				else
!					if(pure_inside_check) then
!						IND_inside(COUNTER_inside+1) = ALL_PARTICLES(IND_overlap(i))%ID
!						COUNTER_inside = COUNTER_inside + 1
!						
!						IND_overlap(i:(COUNTER_overlap-1)) = IND_overlap((i+1):COUNTER_overlap)
!						COUNTER_overlap = COUNTER_overlap - 1
!					else 
!						IND_leave(COUNTER_leave+1) = ALL_PARTICLES(IND_overlap(i))%ID
!						COUNTER_leave = COUNTER_leave + 1
!					end if 
!				end if
!			END DO
!			
!				
!		END SUBROUTINE update_table_block
!		
!		SUBROUTINE update_all_particle_information
!			DO i=1,COUNTER_inside
!!				print*, 'Xparticle_read(1,IND_inside(i))',Xparticle_read(1,IND_inside(i))
!!				print*, 'ALL_PARTICLES(IND_inside(i))%Xp',ALL_PARTICLES(IND_inside(i))%Xp
!!				print*, "Before",Xparticle_read(1,IND_inside(i))
!				Xparticle_read(1,IND_inside(i)) = ALL_PARTICLES(IND_inside(i))%Xp
!				Xparticle_read(2,IND_inside(i)) = ALL_PARTICLES(IND_inside(i))%Yp
!!				print*, "After", Xparticle_read(1,IND_inside(i))
!				D_read(1,i) = ALL_PARTICLES(IND_inside(i))%Dp1
!				D_read(2,i) = ALL_PARTICLES(IND_inside(i))%Dp2
!				D_read(3,i) = ALL_PARTICLES(IND_inside(i))%Dp3 
!				D_read(4,i) = ALL_PARTICLES(IND_inside(i))%Dp4
!!				print*,'AAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
!!				print*, 'Xparticle_read(1,IND_inside(i))',Xparticle_read(1,IND_inside(i))
!!				print*, 'ALL_PARTICLES(IND_inside(i))%Xp',ALL_PARTICLES(IND_inside(i))%Xp
!			END DO
!			
!			DO i=1,COUNTER_overlap
!				Xparticle_read(1,IND_overlap(i)) = ALL_PARTICLES(IND_overlap(i))%Xp
!				Xparticle_read(2,IND_overlap(i)) = ALL_PARTICLES(IND_overlap(i))%Yp
!				
!				D_read(1,i) = ALL_PARTICLES(IND_overlap(i))%Dp1
!				D_read(2,i) = ALL_PARTICLES(IND_overlap(i))%Dp2
!				D_read(3,i) = ALL_PARTICLES(IND_overlap(i))%Dp3 
!				D_read(4,i) = ALL_PARTICLES(IND_overlap(i))%Dp4
!			END DO
!			
!			
!			write(2) Xparticle_read

!			write(3) D_read
!			
!			
!			close(2)
!			close(3)

!		END SUBROUTINE update_all_particle_information


		SUBROUTINE dealloc_all_particle_table
			DEALLOCATE(IND_leave)
			DEALLOCATE(IND_inside)
			DEALLOCATE(IND_overlap)
			DEALLOCATE(ALL_PARTICLES)
			
			DEALLOCATE(COUNTER_leave)
			DEALLOCATE(COUNTER_overlap)
			DEALLOCATE(IND_recv_overlap)
		
			DEALLOCATE(Dout)
			DEALLOCATE(Xpart_block)
			DEALLOCATE(Mblock)
			DEALLOCATE(Dblock)
			DEALLOCATE(velocity_field)	
		
		END SUBROUTINE dealloc_all_particle_table

END MODULE mpi_2d_structures_modf90
