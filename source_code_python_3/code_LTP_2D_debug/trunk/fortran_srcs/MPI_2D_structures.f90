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
	DOUBLE PRECISION                                      :: Norm_inf_Dm1
    INTEGER(8)                                            :: indice_max_norm_Dm1
    DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE        :: Dout


	INTEGER                                               :: Npart_block	
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)         :: Xpart_block
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)           :: Mblock
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)         :: Dblock
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)         :: velocity_field
		

	! TODO derived particle dataypes

	INTEGER, DIMENSION(:), ALLOCATABLE                    :: IND_LEAVE
	INTEGER, DIMENSION(:), ALLOCATABLE                    :: IND_INSIDE
	INTEGER, DIMENSION(:), ALLOCATABLE                    :: IND_OVERLAP
	INTEGER, DIMENSION(:), ALLOCATABLE                    :: overlap_north, & 
															 overlap_south, & 
															 overlap_west, & 
															 overlap_east, & 
															 overlap_north_east, & 
															 overlap_north_west, & 
															 overlap_south_east, & 
															 overlap_south_west
	
	TYPE(PARTTYPE), DIMENSION(:), ALLOCATABLE :: ALL_PARTICLES
	
	
	INTEGER                                               :: COUNTER_inside, COUNTER_overlap
	INTEGER                                               :: COUNTER_leave !, DIMENSION(:,:), ALLOCATABLE
	
	INTEGER                                               :: COUNTER_overlap_north, & 
															 COUNTER_overlap_south, & 
															 COUNTER_overlap_east, & 
															 COUNTER_overlap_west, & 
															 COUNTER_overlap_north_east, & 
															 COUNTER_overlap_north_west, & 
															 COUNTER_overlap_south_east, & 
															 COUNTER_overlap_south_west
	CONTAINS
		
		SUBROUTINE initiation_table
			ALLOCATE(IND_LEAVE(number_of_particles))
			ALLOCATE(IND_INSIDE(number_of_particles))
			ALLOCATE(IND_OVERLAP(number_of_particles))
			ALLOCATE(ALL_PARTICLES(number_of_particles))

			ALLOCATE(overlap_north(number_of_particles))
			ALLOCATE(overlap_south(number_of_particles))
			ALLOCATE(overlap_west(number_of_particles))
			ALLOCATE(overlap_east(number_of_particles))
			ALLOCATE(overlap_north_east(number_of_particles))
			ALLOCATE(overlap_north_west(number_of_particles))
			ALLOCATE(overlap_south_east(number_of_particles))
			ALLOCATE(overlap_south_west(number_of_particles))
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
		
		! OVERLAP CRITERION FOR JUST ONE PARTICLE
		
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
			logical :: overlapcheck_1, insidecheck_1
			INTEGER :: position_inside, position_overlap
			DOUBLE PRECISION, dimension(2)             :: pointu1, pointu2,pointu3,pointu4
			
			
			position_overlap = 1
			position_inside = 1
			COUNTER_inside = 0
			COUNTER_overlap = 0
			COUNTER_leave = 0 
			
			!for every particle coordinates in Xparticle, assign it rightly to the right domain
			DO i =1,number_of_particles
				CALL overlap_criterion(ALL_PARTICLES(i),overlapcheck_1,insidecheck_1,pointu1,pointu2,pointu3,pointu4)
				
				ALL_PARTICLES(i)%pointu1 = pointu1
				ALL_PARTICLES(i)%pointu2 = pointu2
				ALL_PARTICLES(i)%pointu3 = pointu3
				ALL_PARTICLES(i)%pointu4 = pointu4
				

				if (insidecheck_1) then
					IND_INSIDE(position_inside) = i
					position_inside = position_inside +1
					COUNTER_inside = COUNTER_inside + 1

				end if
				
				if (overlapcheck_1) then
					IND_OVERLAP(position_overlap) = i
					position_overlap = position_overlap +1
					COUNTER_overlap = COUNTER_overlap + 1
				end if
			END DO
		END SUBROUTINE particle_distribution
		
		
		SUBROUTINE particle_screening
			! For an overlap table of a block, we will precise in which neighbour block, these overlap
			! particles access.
			
			! Idea, we make a loop in overlap table, and condition for an overlap particle access in
			! a neighbour block is JUST ONE of 4 pointu locate inside of the neigbour.
						
			double precision, dimension(2) :: p1,p2,p3,p4
			
				COUNTER_overlap_north = 0
				COUNTER_overlap_south = 0	
				COUNTER_overlap_east   = 0
 				COUNTER_overlap_west       = 0
				COUNTER_overlap_north_east  = 0
				COUNTER_overlap_north_west = 0
				COUNTER_overlap_south_east  = 0
				COUNTER_overlap_south_west = 0
							
			DO i = 1, COUNTER_overlap
				p1 = ALL_PARTICLES(IND_OVERLAP(i))%pointu1
				p2 = ALL_PARTICLES(IND_OVERLAP(i))%pointu2
				p3 = ALL_PARTICLES(IND_OVERLAP(i))%pointu3
				p4 = ALL_PARTICLES(IND_OVERLAP(i))%pointu4

		
				
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
						overlap_north(COUNTER_overlap_north) = ALL_PARTICLES(IND_OVERLAP(i))%ID
						COUNTER_overlap_north = COUNTER_overlap_north + 1 
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
						overlap_south(COUNTER_overlap_south) = ALL_PARTICLES(IND_OVERLAP(i))%ID
						COUNTER_overlap_south = COUNTER_overlap_south + 1
					end if
				
				end if
				if (rank_left /= -1) then 
					if (( (p1(2) > start_y .and. p1(2) < end_y) .and. (p1(1) > start_x-block_step_x .and. p1(1) < end_y-block_step_x) ) .or. &
					 ( (p2(2) > start_y .and. p2(2) < end_y) .and. (p2(1) > start_x-block_step_x .and. p2(1) < end_y-block_step_x) ) .or. &
					 ( (p3(2) > start_y .and. p3(2) < end_y) .and. (p3(1) > start_x-block_step_x .and. p3(1) < end_y-block_step_x) ) .or. &
					 ( (p4(2) > start_y .and. p4(2) < end_y) .and. (p4(1) > start_x-block_step_x .and. p4(1) < end_y-block_step_x) )) &
					 
					then
						overlap_west(COUNTER_overlap_west) = ALL_PARTICLES(IND_OVERLAP(i))%ID
						COUNTER_overlap_west = COUNTER_overlap_west + 1
					
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
					 	overlap_east(COUNTER_overlap_east) = ALL_PARTICLES(IND_OVERLAP(i))%ID
						COUNTER_overlap_east = COUNTER_overlap_east + 1
					
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
			        	overlap_north_west(COUNTER_overlap_east) = ALL_PARTICLES(IND_OVERLAP(i))%ID
						COUNTER_overlap_north_west = COUNTER_overlap_north_west + 1
					
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
			        	overlap_north_east(COUNTER_overlap_east) = ALL_PARTICLES(IND_OVERLAP(i))%ID
						COUNTER_overlap_north_east = COUNTER_overlap_north_east + 1
			        
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
						overlap_south_west(COUNTER_overlap_east) = ALL_PARTICLES(IND_OVERLAP(i))%ID
						COUNTER_overlap_south_west = COUNTER_overlap_south_west + 1
					
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
						overlap_south_east(COUNTER_overlap_east) = ALL_PARTICLES(IND_OVERLAP(i))%ID
						COUNTER_overlap_south_east = COUNTER_overlap_south_east + 1
					
					end if
				end if
			
			END DO
		END SUBROUTINE
		
		SUBROUTINE update_overlap_table
		! RECEIVE overlap particles from all possible neighbour blocks
		! TODO after particle_screening, we will know in every blocks, which overlap particle
		! will overlap which block	
		
		! Idea, we concatenate all possible overlap particles from all possible neighbour
		! blocks into IND_OVERLAP
		
		! Firstly, send its informations about "overlap_direction" to all neighbour blocks
		! Then receive informations from its neighbour blocks.
		
		! MPI.SEND
		! MPI.RECV
		
		
		
		END SUBROUTINE update_overlap_table
		
		
		SUBROUTINE block_loop	
		
		IMPLICIT NONE
		! Input: X_block_k (k is rank of processor),M_block_k,D_old_block_k,hx,hy, N_part_block_k	
		! Output: U_block_k, Df2py_block_k
		Npart_block = COUNTER_inside + COUNTER_overlap

		DO i= 1,COUNTER_inside
			Xpart_block(1,i) = ALL_PARTICLES(IND_INSIDE(i))%Xp
			Xpart_block(2,i) = ALL_PARTICLES(IND_INSIDE(i))%Yp
		END DO
		DO i=1,COUNTER_overlap
			Xpart_block(1,i+COUNTER_inside) = ALL_PARTICLES(IND_OVERLAP(i))%Xp
			Xpart_block(2,i+COUNTER_inside) = ALL_PARTICLES(IND_OVERLAP(i))%Yp
		END DO
		
		
		
		DO i= 1,COUNTER_inside
			Mblock(i) = ALL_PARTICLES(IND_INSIDE(i))%mass
		END DO
		DO i=1,COUNTER_overlap
			Mblock(i+COUNTER_inside) = ALL_PARTICLES(IND_OVERLAP(i))%mass
		END DO

		
		DO i= 1,COUNTER_inside
			Dblock(1,i) = ALL_PARTICLES(IND_INSIDE(i))%Dp1
			Dblock(2,i) = ALL_PARTICLES(IND_INSIDE(i))%Dp2
			Dblock(3,i) = ALL_PARTICLES(IND_INSIDE(i))%Dp3
			Dblock(4,i) = ALL_PARTICLES(IND_INSIDE(i))%Dp4
		END DO
		DO i=1,COUNTER_overlap
			Dblock(1,i+COUNTER_inside) = ALL_PARTICLES(IND_OVERLAP(i))%Dp1
			Dblock(2,i+COUNTER_inside) = ALL_PARTICLES(IND_OVERLAP(i))%Dp2
			Dblock(3,i+COUNTER_inside) = ALL_PARTICLES(IND_OVERLAP(i))%Dp3
			Dblock(4,i+COUNTER_inside) = ALL_PARTICLES(IND_OVERLAP(i))%Dp4
		END DO
		
		
		DO i= 1,COUNTER_inside
			velocity_field(1,i) = ALL_PARTICLES(IND_INSIDE(i))%Upx
			velocity_field(2,i) = ALL_PARTICLES(IND_INSIDE(i))%Upy
		END DO
		DO i=1,COUNTER_overlap
			velocity_field(1,i+COUNTER_inside) = ALL_PARTICLES(IND_OVERLAP(i))%Upx
			velocity_field(2,i+COUNTER_inside) = ALL_PARTICLES(IND_OVERLAP(i))%Upy
		END DO
		
!		print*, shape(velocity_field)
!		print*,velocity_field(:,1:10)
!		print*,Xpart_block(:,1:10)
!		print*, Mblock(1:10)
!		print*, hx
		! Velocity calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!HUGE PROBLEM
!Need to initiate all the value in new library calculsfortran_rec,init,2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		call diffusion_field_ltp(Xpart_block,Mblock,Dblock,hx,hy,velocity_field,Npart_block)
!		call update_d_diffusion(Xpart_block, Dblock, Mblock, Tini, dt, time_scheme, hx, hy, &
!     indice_max_norm_Dm1, Norm_inf_Dm1, Dout, Npart_block)
     	
		END SUBROUTINE block_loop	
		
		SUBROUTINE update_displacement
			DO i=1,COUNTER_inside
				ALL_PARTICLES(IND_INSIDE(i))%Xp = ALL_PARTICLES(IND_INSIDE(i))%Xp - velocity_field(1,i)*dt
				ALL_PARTICLES(IND_INSIDE(i))%Yp = ALL_PARTICLES(IND_INSIDE(i))%Yp - velocity_field(2,i)*dt
			END DO
			DO i=1,COUNTER_overlap
				ALL_PARTICLES(IND_OVERLAP(i))%Xp = ALL_PARTICLES(IND_OVERLAP(i))%Xp - velocity_field(1,i+COUNTER_inside)*dt
				ALL_PARTICLES(IND_OVERLAP(i))%Yp = ALL_PARTICLES(IND_OVERLAP(i))%Yp - velocity_field(2,i+COUNTER_inside)*dt
			END DO

		END SUBROUTINE update_displacement
		
		
		SUBROUTINE update_deformation_matrix
			DO i=1,COUNTER_inside
				ALL_PARTICLES(IND_INSIDE(i))%Dp1 = Dout(1,i) 
				ALL_PARTICLES(IND_INSIDE(i))%Dp2 = Dout(2,i)
				ALL_PARTICLES(IND_INSIDE(i))%Dp3 = Dout(3,i)
				ALL_PARTICLES(IND_INSIDE(i))%Dp4 = Dout(4,i) 
			END DO
			DO i=1,COUNTER_overlap
				ALL_PARTICLES(IND_OVERLAP(i))%Dp1 = Dout(1,i+COUNTER_inside) 
				ALL_PARTICLES(IND_OVERLAP(i))%Dp2 = Dout(2,i+COUNTER_inside)
				ALL_PARTICLES(IND_OVERLAP(i))%Dp3 = Dout(3,i+COUNTER_inside)
				ALL_PARTICLES(IND_OVERLAP(i))%Dp4 = Dout(4,i+COUNTER_inside)
			END DO
		END SUBROUTINE update_deformation_matrix
		
		
		SUBROUTINE update_table_block
			IMPLICIT NONE
			LOGICAL :: overlap_check, pure_inside_check
			double precision, dimension(2) :: pointu1,pointu2,pointu3,pointu4

			DO i=1,COUNTER_inside
				call overlap_criterion(ALL_PARTICLES(IND_INSIDE(i)),overlap_check,pure_inside_check,pointu1,pointu2,pointu3,pointu4)
							
				if (pure_inside_check) then
					cycle
					
				else
					if(overlap_check) then
						COUNTER_overlap = COUNTER_overlap + 1
						
						IND_INSIDE(i:(COUNTER_inside-1)) = IND_INSIDE((i+1):COUNTER_inside)
						COUNTER_inside = COUNTER_inside - 1
					else
						IND_LEAVE(COUNTER_leave+1) = ALL_PARTICLES(IND_INSIDE(i))%ID
						COUNTER_leave = COUNTER_leave + 1
					end if
				end if
			END DO		
			
			DO i=1, COUNTER_overlap
				call overlap_criterion(ALL_PARTICLES(IND_OVERLAP(i)),overlap_check,pure_inside_check,pointu1,pointu2,pointu3,pointu4)
				
				if (overlap_check) then
					cycle
				else
					if(pure_inside_check) then
						IND_INSIDE(COUNTER_inside+1) = ALL_PARTICLES(IND_OVERLAP(i))%ID
						COUNTER_inside = COUNTER_inside + 1
						
						IND_OVERLAP(i:(COUNTER_overlap-1)) = IND_OVERLAP((i+1):COUNTER_overlap)
						COUNTER_overlap = COUNTER_overlap - 1
					else 
						IND_LEAVE(COUNTER_leave+1) = ALL_PARTICLES(IND_OVERLAP(i))%ID
						COUNTER_leave = COUNTER_leave + 1
					end if 
				end if
			END DO
			
				
		END SUBROUTINE update_table_block
		
		SUBROUTINE update_particle_information
			DO i=1,COUNTER_inside
!				print*, 'Xparticle_read(1,IND_INSIDE(i))',Xparticle_read(1,IND_INSIDE(i))
!				print*, 'ALL_PARTICLES(IND_INSIDE(i))%Xp',ALL_PARTICLES(IND_INSIDE(i))%Xp
				Xparticle_read(1,IND_INSIDE(i)) = ALL_PARTICLES(IND_INSIDE(i))%Xp
				Xparticle_read(2,IND_INSIDE(i)) = ALL_PARTICLES(IND_INSIDE(i))%Yp
				
				D_read(1,i) = ALL_PARTICLES(IND_INSIDE(i))%Dp1
				D_read(2,i) = ALL_PARTICLES(IND_INSIDE(i))%Dp2
				D_read(3,i) = ALL_PARTICLES(IND_INSIDE(i))%Dp3 
				D_read(4,i) = ALL_PARTICLES(IND_INSIDE(i))%Dp4
!				print*,'AAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
!				print*, 'Xparticle_read(1,IND_INSIDE(i))',Xparticle_read(1,IND_INSIDE(i))
!				print*, 'ALL_PARTICLES(IND_INSIDE(i))%Xp',ALL_PARTICLES(IND_INSIDE(i))%Xp
			END DO
			
			DO i=1,COUNTER_overlap
				Xparticle_read(1,IND_OVERLAP(i)) = ALL_PARTICLES(IND_OVERLAP(i))%Xp
				Xparticle_read(2,IND_OVERLAP(i)) = ALL_PARTICLES(IND_OVERLAP(i))%Yp
				
				D_read(1,i) = ALL_PARTICLES(IND_OVERLAP(i))%Dp1
				D_read(2,i) = ALL_PARTICLES(IND_OVERLAP(i))%Dp2
				D_read(3,i) = ALL_PARTICLES(IND_OVERLAP(i))%Dp3 
				D_read(4,i) = ALL_PARTICLES(IND_OVERLAP(i))%Dp4
			END DO
			
			
			write(2) Xparticle_read

			write(3) D_read
			
			
			close(2)
			close(3)

		END SUBROUTINE update_particle_information


		SUBROUTINE dealloc_all_particle_table
			DEALLOCATE(IND_LEAVE)
			DEALLOCATE(IND_INSIDE)
			DEALLOCATE(IND_OVERLAP)
			DEALLOCATE(ALL_PARTICLES)
			DEALLOCATE(overlap_north)
			DEALLOCATE(overlap_south)
			DEALLOCATE(overlap_west)
			DEALLOCATE(overlap_east)
			DEALLOCATE(overlap_north_east)
			DEALLOCATE(overlap_north_west)
			DEALLOCATE(overlap_south_east)
			DEALLOCATE(overlap_south_west)
			
			
			
			DEALLOCATE(Dout)
			DEALLOCATE(Xpart_block)
			DEALLOCATE(Mblock)
			DEALLOCATE(Dblock)
			DEALLOCATE(velocity_field)	
		
		END SUBROUTINE dealloc_all_particle_table

END MODULE mpi_2d_structures_modf90
