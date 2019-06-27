!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> MPI 2D structures 
!--------------------------------------------------------------------------- 


MODULE mpi_2d_structures_modf90
	USE mod_particle_2D_modf90
	USE PACKMPI
	USE data_launch
	USE Jacobi_method
	USE calculsfor_rec_modf90
	IMPLICIT NONE
	
	! DECALARATIONS
	! TODO

	INTEGER                                               :: Npart_block	
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)         :: Xpart_block
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)           :: Mblock
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)         :: Dblock
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)         :: velocity_field
		

	! TODO derived particle dataypes
	TYPE(PARTTYPE), DIMENSION(:), ALLOCATABLE :: inside_particle_table
	TYPE(PARTTYPE), DIMENSION(:), ALLOCATABLE :: overlap_particle_table
	TYPE(PARTTYPE), DIMENSION(:), ALLOCATABLE :: ALL_PARTICLES
	
	INTEGER                    ::     COUNTER_inside, COUNTER_overlap
	CONTAINS
		
		SUBROUTINE initiation_table
			ALLOCATE(inside_particle_table(number_of_particles))
			ALLOCATE(overlap_particle_table(number_of_particles))
			ALLOCATE(ALL_PARTICLES(number_of_particles))
		END SUBROUTINE initiation_table
		
		SUBROUTINE neighbour_counter
			nb_neighbours = 0
			if (rank_left .ge. 0) then 
				nb_neighbours = nb_neighbours+ 1
			end if
			if (rank_right .ge. 0) then 
			nb_neighbours =nb_neighbours+ 1
			end if
			if (rank_up .ge. 0)then
				nb_neighbours =nb_neighbours + 1
			end if
			if (rank_down .ge. 0)then
				nb_neighbours = nb_neighbours+ 1
			end if
			if (rank_up_left .ge. 0)then
				nb_neighbours = nb_neighbours+1
			end if
			if (rank_up_right .ge. 0)then
				nb_neighbours =nb_neighbours+ 1
			end if
			if (rank_down_left .ge. 0)then 
				nb_neighbours =nb_neighbours+ 1
			end if
			if (rank_down_right  .ge. 0) then 
				nb_neighbours =nb_neighbours+ 1
			end if		
		END SUBROUTINE neighbour_counter
		
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
				! TODO ALL_PARTICLES(i)%Mp =
				ALL_PARTICLES(i)%ID = i
				ALL_PARTICLES(i)%Rhop = density_read(i)  			
			END DO
		
		END SUBROUTINE parttype_convert
		
		! OVERLAP CRITERION FOR JUST ONE PARTICLE
		
		SUBROUTINE overlap_criterion(particle_k,overlap_check,inside_check)
		! Check a particle if it satifies overlap_criterion
		! Retrieve the eigenvalues and eigenvectors to determine the direction movement of particles, 	then decide its status
		
		! Input: k-th_particle in "PARTYPE" type 
		! Output: logical
		
		
			IMPLICIT NONE
			TYPE(PARTTYPE), INTENT(in) :: particle_k
			LOGICAL, INTENT(out)       :: overlap_check
			LOGICAL, INTENT(out)       :: inside_check			
			! Initial idea
			INTEGER, PARAMETER                                      :: n=2
			DOUBLE PRECISION                                        :: a(n,n), x(n,n)
			DOUBLE PRECISION, PARAMETER                             :: abserr=1.0e-09
			INTEGER i, j


			double precision, dimension(n)             	            :: eigenvector_1
			double precision, dimension(n)                          :: eigenvector_2
			DOUBLE PRECISION                                        :: eigenvalue_1
			DOUBLE PRECISION                                        :: eigenvalue_2
			DOUBLE PRECISION, dimension(n)                          :: pointu1, pointu2,pointu3,pointu4,pointu5,pointu6,pointu7,pointu8

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
				
			pointu1(1) = eigenvector_1(1) + eigenvector_2(1) + particle_k%Xp
			pointu1(2) = eigenvector_1(2) + eigenvector_2(2) + particle_k%Yp
			pointu2(1) = -eigenvector_1(1)- eigenvector_2(1) + particle_k%Xp
			pointu2(2) = -eigenvector_1(2) - eigenvector_2(2) + particle_k%Yp
			pointu3(1) = -eigenvector_1(1)+ eigenvector_2(1) + particle_k%Xp
			pointu3(2) = -eigenvector_1(2) + eigenvector_2(2) + particle_k%Yp
			pointu4(1) = eigenvector_1(1) - eigenvector_2(1) + particle_k%Xp
			pointu4(2) = eigenvector_1(2) - eigenvector_2(2) + particle_k%Yp
			
			
			pointu5(1) = eigenvector_1(1) + particle_k%Xp
			pointu5(2) = eigenvector_1(2) + particle_k%Yp
			pointu6(1) = -eigenvector_1(1) + particle_k%Xp
			pointu6(2) = -eigenvector_1(2) + particle_k%Yp
			pointu7(1) = eigenvector_2(1) + particle_k%Xp
			pointu7(2) = eigenvector_2(2) + particle_k%Yp
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
			pointu5inside = (pointu5(1)-start_x>0) &
				.and.(pointu5(1)-end_x<0)          &
				.and.(pointu5(2)-start_y>0)        &
				.and.(pointu5(2)-end_y<0)
			pointu6inside = (pointu6(1)-start_x>0) &
				.and.(pointu6(1)-end_x<0)          &
				.and.(pointu6(2)-start_y>0)        &
				.and.(pointu6(2)-end_y<0)
			pointu7inside = (pointu7(1)-start_x>0) &
				.and.(pointu7(1)-end_x<0)          &
				.and.(pointu7(2)-start_y>0)        &
				.and.(pointu7(2)-end_y<0)
			pointu8inside = (pointu8(1)-start_x>0) &
				.and.(pointu8(1)-end_x<0)          &
				.and.(pointu8(2)-start_y>0)        &
				.and.(pointu8(2)-end_y<0)
			
			
			
			if(pointu1inside.and.pointu2inside.and.pointu3inside.and.pointu4inside)then 
				inside_check = .true. ! This particle is inside block
			else 
				inside_check = .false.	! This particle is maybe an overlap particle
			end if
			
			if ((.not. inside_check).and.(pointu1inside .or. &
										pointu2inside .or. &
										pointu3inside .or. &
										pointu4inside .or. &
										pointu5inside .or. &
										pointu6inside .or. &
										pointu7inside .or. &
										pointu8inside )) then
				overlap_check = .true.
			else
				overlap_check = .false.
			end if 

		END SUBROUTINE overlap_criterion
	
		SUBROUTINE particle_screening()
			!table_block_particle
			!type(PARTYPE), dimension(:,:), allocatable, intent(in)  :: table_block_particle
			!type(PARTYPE), dimension(:,:), allocatable, intent(out) :: table_particle_inside
			!type(PARTYPE), dimension(:,:), allocatable, intent(out) :: table_particle_overlap
			! Loop on ALL_PARTICLES
			! Check overlap_criterion
			! overlapcheck false: append to table_particle_inside
			! overlapcheck true : TODO
		END SUBROUTINE
		
		SUBROUTINE particle_distribution
			IMPLICIT NONE
			logical :: overlapcheck_1, insidecheck_1
			INTEGER :: position_inside, position_overlap
			
			position_overlap = 1
			position_inside = 1
			COUNTER_inside = 0
			COUNTER_overlap = 0

			!for every particle coordinates in Xparticle, assign it rightly to the right domain
			DO i =1,number_of_particles
				CALL overlap_criterion(ALL_PARTICLES(i),overlapcheck_1,insidecheck_1)
				!print*, overlapcheck_1
				
!				if(((ALL_PARTICLES(i)%Xp > start_x ) &
!				.and. (ALL_PARTICLES(i)%Xp <= end_x)) &
!				.and. ((ALL_PARTICLES(i)%Yp > start_y )&
!				.and.(ALL_PARTICLES(i)%Yp <= end_y)) &
!				.and. (.not. overlapcheck_1) & 
!				) 
				if (insidecheck_1) then
					inside_particle_table(position_inside) = ALL_PARTICLES(i)
					inside_particle_table(i)%PROC_ID = rank
					COUNTER_inside = COUNTER_inside + 1
					position_inside = position_inside +1
				end if
				
				if (overlapcheck_1) then
					overlap_particle_table(position_overlap) = ALL_PARTICLES(i)
					COUNTER_overlap = COUNTER_overlap + 1
					position_overlap = position_overlap +1
				end if
				
			END DO
		END SUBROUTINE particle_distribution
		
		
		SUBROUTINE block_loop	
		
		IMPLICIT NONE
		! TODO
		! Input: X_block_k (k is rank of processor),M_block_k,D_old_block_k,hx,hy, N_part_block_k	
		! Output: U_block_k, Df2py_block_k
		
		
		
		Npart_block = COUNTER_inside + COUNTER_overlap
		
		ALLOCATE(Xpart_block(2,Npart_block))
		ALLOCATE(Mblock(Npart_block))
		ALLOCATE(Dblock(4,Npart_block))
		ALLOCATE(velocity_field(2,Npart_block))
		
		
		DO i= 1,COUNTER_inside
			Xpart_block(1,i) = inside_particle_table(i)%Xp
			Xpart_block(2,i) = inside_particle_table(i)%Yp
		END DO
		DO i=1,COUNTER_overlap
			Xpart_block(1,i+COUNTER_inside) = overlap_particle_table(i)%Xp
			Xpart_block(2,i+COUNTER_inside) = overlap_particle_table(i)%Yp
		END DO
		
		
		
		DO i= 1,COUNTER_inside
			Mblock(i) = inside_particle_table(i)%Mp
		END DO
		DO i=1,COUNTER_overlap
			Mblock(i+COUNTER_inside) = overlap_particle_table(i)%Mp
		END DO
		
		
		
		DO i= 1,COUNTER_inside
			Dblock(1,i) = inside_particle_table(i)%Dp1
			Dblock(2,i) = inside_particle_table(i)%Dp2
			Dblock(3,i) = inside_particle_table(i)%Dp3
			Dblock(4,i) = inside_particle_table(i)%Dp4
		END DO
		DO i=1,COUNTER_overlap
			Dblock(1,i+COUNTER_inside) = overlap_particle_table(i)%Dp1
			Dblock(2,i+COUNTER_inside) = overlap_particle_table(i)%Dp2
			Dblock(3,i+COUNTER_inside) = overlap_particle_table(i)%Dp3
			Dblock(4,i+COUNTER_inside) = overlap_particle_table(i)%Dp4
		END DO
		
		
		DO i= 1,COUNTER_inside
			velocity_field(1,i) = inside_particle_table(i)%Upx
			velocity_field(2,i) = inside_particle_table(i)%Upy
		END DO
		DO i=1,COUNTER_overlap
			velocity_field(1,i+COUNTER_inside) = overlap_particle_table(i)%Upx
			velocity_field(2,i+COUNTER_inside) = overlap_particle_table(i)%Upy
		END DO
		
		! Velocity calculation
		call diffusion_field_ltp(Xpart_block,Mblock,Dblock,hx_remap,hy_remap,velocity_field,Npart_block)
	
		END SUBROUTINE block_loop	
		
		SUBROUTINE update_displacement
			DO i=1,COUNTER_inside
				inside_particle_table(i)%Xp = inside_particle_table(i)%Xp + velocity_field(1,i)*dt
				inside_particle_table(i)%Yp = inside_particle_table(i)%Yp + velocity_field(2,i)*dt
			END DO
			DO i=1,COUNTER_overlap
				overlap_particle_table(i)%Xp = overlap_particle_table(i)%Xp + velocity_field(1,i+COUNTER_inside)*dt
				overlap_particle_table(i)%Yp = overlap_particle_table(i)%Yp + velocity_field(2,i+COUNTER_inside)*dt
			END DO
		END SUBROUTINE update_displacement
		
		SUBROUTINE update_table
			
		END SUBROUTINE update_table



END MODULE mpi_2d_structures_modf90
