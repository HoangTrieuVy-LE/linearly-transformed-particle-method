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
	IMPLICIT NONE
	
	! DECALARATIONS
	! TODO
	
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: IND_OVERLAP
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: IND_INSIDE
	
	! TODO derived particle dataypes
	TYPE(PARTTYPE), DIMENSION(:), ALLOCATABLE :: inside_particle_table
	TYPE(PARTTYPE), DIMENSION(:,:), ALLOCATABLE :: overlap_particle_table
	TYPE(PARTTYPE), DIMENSION(:), ALLOCATABLE :: ALL_PARTICLES
	
	INTEGER                    ::     COUNTER
	CONTAINS
		
		SUBROUTINE initiation_table
			ALLOCATE(inside_particle_table(number_of_particles))
			ALLOCATE(overlap_particle_table(number_of_particles,nb_neighbours))
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
		
		SUBROUTINE particle_distribution
			COUNTER = 0
			!for every particle coordinates in Xparticle, assign it rightly to the right domain
			DO i =1,number_of_particles
				if(((ALL_PARTICLES(i)%Xp .gt. start_x ).and. (ALL_PARTICLES(i)%Xp .le. end_x))&
				.and. ((ALL_PARTICLES(i)%Yp .gt. start_y )&
				.and.(ALL_PARTICLES(i)%Yp .le. end_y))) THEN
					inside_particle_table(i) = ALL_PARTICLES(i)
					inside_particle_table(i)%PROC_ID = rank
					COUNTER = COUNTER + 1
				end if
				
			END DO
		END SUBROUTINE particle_distribution
		
	
		SUBROUTINE initiation_overlap_table
		! Create overlap table for each process 
		
		END SUBROUTINE initiation_overlap_table


		SUBROUTINE block_loop	
		! TODO

		! Input: X_casier_k (k is rank of processor),M_casier_k,D_old_casier_k,hx,hy, N_part_casier_k	
		! Output: U_casier_k, Df2py_casier_k	
	
		END SUBROUTINE block_loop	
		
		SUBROUTINE overlap_criterion(particle_k,overlap_check)
		! Check a particle if it satifies overlap_criterion
		! Retrieve the eigenvalues and eigenvectors to determine the direction movement of particles, 	then decide it leaves the block or not
		! Input: k-th_particle in "PARTYPE" type 
		! Output: logical
			IMPLICIT NONE
			TYPE(PARTTYPE), INTENT(in) :: particle_k
			LOGICAL, INTENT(out)       :: overlap_check
			
			! Initial idea
			INTEGER, PARAMETER                                      :: n=2
			DOUBLE PRECISION                                        :: a(n,n), x(n,n)
			DOUBLE PRECISION, PARAMETER                             :: abserr=1.0e-09
			INTEGER i, j


			double precision, dimension(n)             	            :: eigenvector_1
			double precision, dimension(n)                          :: eigenvector_2
			DOUBLE PRECISION                                        :: eigenvalue_1
			DOUBLE PRECISION                                        :: eigenvalue_2
			DOUBLE PRECISION, dimension(n)                          :: pointu1, pointu2,pointu3,pointu4

			DOUBLE PRECISION		   :: D11,D12,D21,D22
			
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
			pointu1(1) = eigenvector_1(1) + eigenvector_2(1)
			pointu1(2) = eigenvector_1(2) + eigenvector_2(2)
			pointu2(1) = -eigenvector_1(1)- eigenvector_2(1)
			pointu2(2) = eigenvector_1(2) + eigenvector_2(2)
			pointu3(1) = -eigenvector_1(1)- eigenvector_2(1)
			pointu3(2) = -eigenvector_1(2) - eigenvector_2(2)
			pointu4(1) = eigenvector_1(1) + eigenvector_2(1)
			pointu4(2) = -eigenvector_1(2) - eigenvector_2(2) 
			
			
			if ((pointu1(1)-start_x>0) &
				.and.(pointu1(1)-end_x<0) &
				.and.(pointu1(2)-start_y>0) &
				.and.(pointu1(2)-end_y<0)) then 
				overlap_check = .false.
			else 
				overlap_check = .true.
			end if
			if ((pointu2(1)-start_x>0) &
				.and.(pointu2(1)-end_x<0) &
				.and.(pointu2(2)-start_y>0) &
				.and.(pointu2(2)-end_y<0)) then 
				overlap_check = .false.
			else 
				overlap_check = .true.
			end if
			if ((pointu3(1)-start_x>0) &
				.and.(pointu3(1)-end_x<0) &
				.and.(pointu3(2)-start_y>0) &
				.and.(pointu3(2)-end_y<0)) then 
				overlap_check = .false.
			else 
				overlap_check = .true.
			end if
			if ((pointu4(1)-start_x>0) &
				.and.(pointu4(1)-end_x<0) &
				.and.(pointu4(2)-start_y>0) &
				.and.(pointu4(2)-end_y<0)) then 
				overlap_check = .false.
			else 
				overlap_check = .true.
			end if
			
		END SUBROUTINE overlap_criterion

		SUBROUTINE particle_screening()
			!table_block_particle
			!type(PARTYPE), dimension(:,:), allocatable, intent(in)  :: table_block_particle
			!type(PARTYPE), dimension(:,:), allocatable, intent(out) :: table_particle_inside
			!type(PARTYPE), dimension(:,:), allocatable, intent(out) :: table_particle_overlap
			! Loop on table_block_particles
			! Check overlap_criterion
			! Yes: append to table_particle_inside
			! No : append to table_particle_overlap
		END SUBROUTINE
		
		
		


END MODULE mpi_2d_structures_modf90
