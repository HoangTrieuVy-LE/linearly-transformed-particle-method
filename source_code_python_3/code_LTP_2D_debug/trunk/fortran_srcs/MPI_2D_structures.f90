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
	
	
	
	! IND_overlap indicates the identity of danger particles in each direction
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: IND_danger
	! IND_recv_danger indicates the identity of particles receving from neighbours
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: IND_recv_danger
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: IND_recv_leave
		
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
	
	! COUNTER_recv_danger indicates the number of overlap particles receiving in each direction 

	INTEGER                                               :: COUNTER_recv_danger
	
	! COUNTER_recv_overlap indicates the number of overlap particles receiving in each direction 

	INTEGER                                               :: COUNTER_recv_overlap
	INTEGER                                               :: COUNTER_recv_leave
	
	! COUNTER_leave indicates the number of particles leaving in each direction
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: COUNTER_leave
	
	! COUNTER_overlap indicates the number of overlap particles in each direction
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: COUNTER_overlap
	! COUNTER_overlap indicates the number of danger particles in each direction
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: COUNTER_danger
	
	
	
	
	
	
	! ALL_PARTICLES stock all particle informations and should be broadcasted after being updated all particles displacement.
	TYPE(PARTTYPE), DIMENSION(:), ALLOCATABLE             :: ALL_PARTICLES
	
	
	
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
			ALLOCATE(IND_danger(number_of_particles,NB_NEIGHBOURS))	
			
			ALLOCATE(IND_recv_overlap(number_of_particles,NB_NEIGHBOURS))
			ALLOCATE(IND_recv_danger(number_of_particles,NB_NEIGHBOURS))
			ALLOCATE(IND_recv_leave(number_of_particles,NB_NEIGHBOURS))
			
			ALLOCATE(COUNTER_danger(NB_NEIGHBOURS,0:nb_proc-1))
			ALLOCATE(COUNTER_overlap(NB_NEIGHBOURS,0:nb_proc-1))
			ALLOCATE(COUNTER_leave(NB_NEIGHBOURS,0:nb_proc-1))
			
			ALLOCATE(ALL_PARTICLES(number_of_particles))
			
					
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
			
		SUBROUTINE overlap_criterion(particle_k,overlap_check,totally_inside_check,danger_check, & 
		pointu1,pointu2,pointu3,pointu4,axe)
		! Check a particle if it satifies overlap_criterion
		! Retrieve the eigenvalues and eigenvectors to determine the direction movement of particles, 	then decide its status
		
		! Input: k-th_particle in "PARTYPE" type 
		! Output: logical
		
		
			IMPLICIT NONE
			TYPE(PARTTYPE), INTENT(in)                              :: particle_k
			LOGICAL, INTENT(out)                                    :: overlap_check
			LOGICAL, INTENT(out)                                    :: totally_inside_check
			LOGICAL, INTENT(out)                                    :: danger_check
			DOUBLE PRECISION, INTENT(out)                           :: axe
			
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

			DOUBLE PRECISION		   :: D11,D12,D21,D22, d1,d2,d3,d4
			
			
			LOGICAL :: pointu1inside,pointu2inside,pointu3inside,pointu4inside
			
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
			
			axe = max(abs(eigenvalue_1),abs(eigenvalue_2))
				
			pointu1(1) =  eigenvector_1(1) + eigenvector_2(1) + particle_k%Xp
			pointu1(2) =  eigenvector_1(2) + eigenvector_2(2) + particle_k%Yp
			pointu2(1) = -eigenvector_1(1)- eigenvector_2(1) + particle_k%Xp
			pointu2(2) = -eigenvector_1(2) - eigenvector_2(2) + particle_k%Yp
			pointu3(1) = -eigenvector_1(1)+ eigenvector_2(1) + particle_k%Xp
			pointu3(2) = -eigenvector_1(2) + eigenvector_2(2) + particle_k%Yp
			pointu4(1) =  eigenvector_1(1) - eigenvector_2(1) + particle_k%Xp
			pointu4(2) =  eigenvector_1(2) - eigenvector_2(2) + particle_k%Yp
			
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
			
			! WE HAVE 3 CASE HERE, ONE IS PARTICLE IS ABSOLUTELY INSIDE ("totally INSIDE") 
			! THE BLOCK HAVING THE RANK OF THE EXECUTING PROCESS OR ANOTHER BLOCK  
			! ANOTHER CASE IS BARYCENTER OF PARTICLE IS INSIDE THE BLOCK AND IT IS 
			! AN OVERLAP PARTICLE
		  
		  if ((particle_k%Xp>start_x .and. particle_k%Xp<end_x) &
		  .and. (particle_k%Yp>start_y .and. particle_k%Yp<end_y)) then
		  
		  	if(pointu1inside.and.pointu2inside.and.pointu3inside.and.pointu4inside)then 
		  		totally_inside_check = .true. ! This particle is inside block
				overlap_check = .false.
				d1 = particle_k%Xp - start_x
				d2 = -particle_k%Xp + end_x
				d3 = particle_k%Yp - start_y
				d4 = -particle_k%Yp + end_y
		  		if (d1<2*axe .or. d2<2*axe .or. d3<2*axe .or. d4 <2*axe) then
					danger_check = .true.
				else
					danger_check = .false.
				end if
			else
				totally_inside_check = .false.	! This particle is maybe an overlap particle
				overlap_check = .true.
				danger_check  = .false.	
		  	end if
		  else 
				overlap_check = .false.
				totally_inside_check = .false.
				danger_check = .false.
				
		  end if

		END SUBROUTINE overlap_criterion
	
		
		
		SUBROUTINE particle_distribution
			IMPLICIT NONE
			logical                           :: local_overlapcheck, local_insidecheck,local_dangercheck,local_leave_check
			DOUBLE PRECISION, dimension(2)    :: pointu1, pointu2,pointu3,pointu4
			DOUBLE PRECISION                  :: axe
			integer :: ID
			
			COUNTER_inside            = 0
			COUNTER_leave(:,rank)     = 0
			COUNTER_overlap(:,rank)   = 0
			COUNTER_danger(:,rank)    = 0

			local_leave_check = .false.
			!for every particle coordinates in Xparticle, assign it rightly to the right domain
			DO i =1,number_of_particles
				CALL overlap_criterion(ALL_PARTICLES(i),local_overlapcheck,local_insidecheck,local_dangercheck & 
				,pointu1,pointu2,pointu3,pointu4,axe)
				
				ALL_PARTICLES(i)%pointu1 = pointu1
				ALL_PARTICLES(i)%pointu2 = pointu2
				ALL_PARTICLES(i)%pointu3 = pointu3
				ALL_PARTICLES(i)%pointu4 = pointu4
				ALL_PARTICLES(i)%axe = axe

				if (local_insidecheck) then
					COUNTER_inside = COUNTER_inside + 1
					IND_inside(COUNTER_inside) = i	
					call particle_screening(ALL_PARTICLES(i),local_overlapcheck,local_dangercheck,local_leave_check)
				else 
					call particle_screening(ALL_PARTICLES(i),local_overlapcheck,local_dangercheck,local_leave_check)
				end if
			END DO
    		
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
    		
    		
    		if(rank /= 0) then
      			call MPI_SEND(COUNTER_danger(1,rank),NB_NEIGHBOURS,MPI_INTEGER,0,111,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1
       				call MPI_RECV(COUNTER_danger(1,ID),NB_NEIGHBOURS,MPI_INTEGER,ID,111,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			end do
    		end if
    		call MPI_BCAST(COUNTER_danger,NB_NEIGHBOURS*nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)
    		
			
		END SUBROUTINE particle_distribution
		
		
		SUBROUTINE particle_screening(particle_k,overlap,danger,leave)
			! For an overlap table of a block, we will precise in which neighbour block, these overlap
			! particles access.
			
			! Idea, we make a loop in overlap table, and condition for an overlap particle access in
			! a neighbour block is JUST ONE of 4 pointu locate inside of the neigbour.
			type(PARTTYPE), INTENT(in)     :: particle_k
			double precision, dimension(2) :: p1,p2,p3,p4
			double precision               :: d1,d2,d3,d4,  axe_k, xk,yk
			logical, intent(in)            :: overlap, danger, leave
			logical                        :: in_up,in_down,in_left,in_right, & 
			in_up_left,in_up_right,in_down_left,in_down_right
			
			
			if (leave) then
				
				xk = particle_k%Xp
				yk = particle_k%Yp
			
				in_up =((yk > start_y+block_step_y .and. yk < end_y+block_step_y) & 
					.and. (xk > start_x .and. xk < end_x))
					
				in_down = ((yk > start_y-block_step_y .and. yk < end_y-block_step_y) &
					.and. (xk > start_x .and. xk < end_x))
				
				in_left = ((yk > start_y .and. yk < end_y) .and. (xk > start_x-block_step_x .and. p1(1) < end_x-block_step_x))
				
				in_right = ((yk > start_y .and. yk < end_y) & 
					.and. (xk > start_x+block_step_x .and. xk < end_x+block_step_x)) 
				
				in_up_left = ((yk > start_y+block_step_y .and. yk < end_y+block_step_y) &
			       .and. (xk > start_x-block_step_x .and. xk < end_x-block_step_x))
			       
				in_up_right = ((yk > start_y+block_step_y .and. yk < end_y+block_step_y) &
			       .and. (xk > start_x+block_step_x .and. xk < end_x+block_step_x))
				
				in_down_left = ((yk > start_y-block_step_y .and. yk < end_y-block_step_y) &
			       .and. (xk > start_x-block_step_x .and. xk < end_x-block_step_x))
				
				in_down_right = ((yk > start_y-block_step_y .and. yk < end_y-block_step_y) &
			       .and. (xk > start_x+block_step_x .and. xk < end_x+block_step_x))
			    
			    if (rank_up /= -1) then
			    	if(in_up) then
			    		COUNTER_leave(FJ,rank) = COUNTER_leave(FJ,rank) + 1
			    		IND_leave(COUNTER_leave(FJ,rank),FJ) = particle_k%ID
			    	end if
			    end if
			    if (rank_down /= -1) then
			    	if(in_down) then
			    		COUNTER_leave(BJ,rank) = COUNTER_leave(BJ,rank) + 1
			    		IND_leave(COUNTER_leave(BJ,rank),BJ) = particle_k%ID
			    	end if
			    end if
			    if (rank_left /= -1) then
			    	if(in_left) then
			    		COUNTER_leave(FK,rank) = COUNTER_leave(FK,rank) + 1
			    		IND_leave(COUNTER_leave(FK,rank),FK) = particle_k%ID
			    	end if
			    end if
			    if (rank_right /= -1) then
			    	if(in_right) then
			    		COUNTER_leave(BK,rank) = COUNTER_leave(BK,rank) + 1
			    		IND_leave(COUNTER_leave(BK,rank),BK) = particle_k%ID
			    	end if
			    end if
			    if (rank_up_left /= -1) then
			    	if(in_up_left) then
			    		COUNTER_leave(FJBK,rank) =COUNTER_leave(FJBK,rank) + 1
			    		IND_leave(COUNTER_leave(FJBK,rank),FJBK) = particle_k%ID
			    	end if
			    end if
			    
			    if (rank_up_right /= -1) then
			    	if(in_up_right) then
			    		COUNTER_leave(FJFK,rank) = COUNTER_leave(FJFK,rank) +1
			    		IND_leave(COUNTER_leave(FJFK,rank),FJFK) = particle_k%ID
			    	end if
			    end if
			    
			    if (rank_down_left /= -1) then
			    	if(in_down_left) then
			    		COUNTER_leave(BJBK,rank) = COUNTER_leave(BJBK,rank) + 1
			    		IND_leave(COUNTER_leave(BJBK,rank),rank) = particle_k%ID
			    	end if 
			    end if
			    
				if (rank_down_right /= -1) then
					if(in_down_right) then
						COUNTER_leave(BJFK,rank) = COUNTER_leave(BJFK,rank) + 1
						IND_leave(COUNTER_leave(BJFK,rank),rank) = particle_k%ID
					end if
				end if
				
			end if

			if(danger) then
				
				d1 = particle_k%Xp - start_x
				d2 = -particle_k%Xp + end_x
				d3 = particle_k%Yp - start_y
				d4 = -particle_k%Yp + end_y
				axe_k = particle_k%axe
				
				
				if (rank_up /= -1) then
				if (d4<2*axe_k .and. (.not. d1<2*axe_k) .and. (.not.d2<2*axe_k)) then
						COUNTER_danger(FJ,rank) = COUNTER_danger(FJ,rank) + 1
						IND_danger(COUNTER_danger(FJ,rank),FJ) =  particle_k%ID
					end if
				end if
				
				if (rank_down /= -1) then
				if (d3<2*axe_k.and. (.not. d1<2*axe_k) .and. (.not.d2<2*axe_k)) then
						COUNTER_danger(BJ,rank) = COUNTER_danger(BJ,rank) + 1
						IND_danger(COUNTER_danger(BJ,rank),BJ) =  particle_k%ID
					end if
				end if
				
				if (rank_left /= -1) then
					if(d1<2*axe_k.and. (.not. d3<2*axe_k) .and. (.not.d4<2*axe_k)) then
						COUNTER_danger(BK,rank) = COUNTER_danger(BK,rank) + 1
						IND_danger(COUNTER_danger(BK,rank),BK) =  particle_k%ID
					end if
				end if
				
				if (rank_right /= -1) then	
					if(d2<2*axe_k.and. (.not. d3<2*axe_k) .and. (.not.d4<2*axe_k)) then
						COUNTER_danger(FK,rank) = COUNTER_danger(FK,rank) + 1
						IND_danger(COUNTER_danger(FK,rank),FK) =  particle_k%ID
					end if
				end if
				
				if (rank_up_left /= -1) then		
					if(d4<2*axe_k .and.d1<2*axe_k) then
						COUNTER_danger(FJBK,rank) = COUNTER_danger(FJBK,rank) + 1
						IND_danger(COUNTER_danger(FJBK,rank),FJBK) =  particle_k%ID
					end if
				end if
				
				if (rank_up_right /= -1) then	    
			        if (d4<2*axe_k .and. d2<2*axe_k) then
			        	COUNTER_danger(FJFK,rank) = COUNTER_danger(FJFK,rank) + 1
						IND_danger(COUNTER_danger(FJFK,rank),FJFK) =  particle_k%ID
			        end if
			    end if
			    
				if (rank_down_left /= -1) then	
					if (d3<2*axe_k .and. d1<2*axe_k) then
						COUNTER_danger(BJBK,rank) = COUNTER_danger(BJBK,rank) + 1
						IND_danger(COUNTER_danger(BJBK,rank),BJBK) =  particle_k%ID
					end if
				end if
					
				if (rank_down_right /= -1) then
					if (d3<2*axe_k .and. d2<2*axe_k) then
						COUNTER_danger(BJFK,rank) = COUNTER_danger(BJFK,rank) + 1
						IND_danger(COUNTER_danger(BJFK,rank),BJFK) =  particle_k%ID
					end if
				end if
			end if
				
			if(overlap) then
				
				p1 = particle_k%pointu1
				p2 = particle_k%pointu2
				p3 = particle_k%pointu3
				p4 = particle_k%pointu4
				
				
				in_up =((p1(2) > start_y+block_step_y .and. p1(2) < end_y+block_step_y) & 
					.and. (p1(1) > start_x .and. p1(1) < end_x)) .or. &
					((p2(2) > start_y+block_step_y .and. p2(2) < end_y+block_step_y) & 
					.and. (p2(1) > start_x .and. p2(1) < end_x)) .or. & 
					((p3(2) > start_y+block_step_y .and. p3(2) < end_y+block_step_y) & 
					.and. (p3(1) > start_x .and. p3(1) < end_x)) .or. &
					((p4(2) > start_y+block_step_y .and. p4(2) < end_y+block_step_y) & 
					.and. (p4(1) > start_x .and. p4(1) < end_x))
					
				in_down = (( (p1(2) > start_y-block_step_y .and. p1(2) < end_y-block_step_y) &
					.and. (p1(1) > start_x .and. p1(1) < end_x) ) .or.&
					 ( (p2(2) > start_y-block_step_y .and. p2(2) < end_y-block_step_y) &
					.and. (p2(1) > start_x .and. p2(1) < end_x) ) .or.&
					 ( (p3(2) > start_y-block_step_y .and. p3(2) < end_y-block_step_y) &
					.and. (p3(1) > start_x .and. p3(1) < end_x) ) .or.&
					 ( (p4(2) > start_y-block_step_y .and. p4(2) < end_y-block_step_y) &
					.and. (p4(1) > start_x .and. p4(1) < end_x) ))
				
				in_left = ( (p1(2) > start_y .and. p1(2) < end_y) .and. (p1(1) > start_x-block_step_x .and. p1(1) < end_x-block_step_x) ) .or. & 
				( (p2(2) > start_y .and. p2(2) < end_y) .and. (p2(1) > start_x-block_step_x .and. p2(1) < end_x-block_step_x) ) .or. & 
				( (p3(2) > start_y .and. p3(2) < end_y) .and. (p3(1) > start_x-block_step_x .and. p3(1) < end_x-block_step_x) ) .or. & 
				( (p4(2) > start_y .and. p4(2) < end_y) .and. (p4(1) > start_x-block_step_x .and. p4(1) < end_x-block_step_x) )
				
				in_right = (( (p1(2) > start_y .and. p1(2) < end_y) & 
					.and. (p1(1) > start_x+block_step_x .and. p1(1) < end_x+block_step_x) ).or. &
					( (p2(2) > start_y .and. p2(2) < end_y) & 
					.and. (p2(1) > start_x+block_step_x .and. p2(1) < end_x+block_step_x) ).or. &
					( (p3(2) > start_y .and. p3(2) < end_y) & 
					.and. (p3(1) > start_x+block_step_x .and. p3(1) < end_x+block_step_x) ).or. &
					( (p4(2) > start_y .and. p4(2) < end_y) & 
					.and. (p4(1) > start_x+block_step_x .and. p4(1) < end_x+block_step_x) )) 
				
				in_up_left = (( (p1(2) > start_y+block_step_y .and. p1(2) < end_y+block_step_y) &
			       .and. (p1(1) > start_x-block_step_x .and. p1(1) < end_x-block_step_x) ).or.&
			       	( (p2(2) > start_y+block_step_y .and. p2(2) < end_y+block_step_y) &
			       .and. (p2(1) > start_x-block_step_x .and. p2(1) < end_x-block_step_x) ).or.&
			       ( (p3(2) > start_y+block_step_y .and. p3(2) < end_y+block_step_y) &
			       .and. (p3(1) > start_x-block_step_x .and. p3(1) < end_x-block_step_x) ).or.&
			       ( (p4(2) > start_y+block_step_y .and. p4(2) < end_y+block_step_y) &
			       .and. (p4(1) > start_x-block_step_x .and. p4(1) < end_x-block_step_x) ))
			       
				in_up_right = (( (p1(2) > start_y+block_step_y .and. p1(2) < end_y+block_step_y) &
			       .and. (p1(1) > start_x+block_step_x .and. p1(1) < end_x+block_step_x) ).or. &
			       ( (p2(2) > start_y+block_step_y .and. p2(2) < end_y+block_step_y) &
			       .and. (p2(1) > start_x+block_step_x .and. p2(1) < end_x+block_step_x) ).or. &
			       ( (p3(2) > start_y+block_step_y .and. p3(2) < end_y+block_step_y) &
			       .and. (p3(1) > start_x+block_step_x .and. p3(1) < end_x+block_step_x) ).or. &
			       ( (p4(2) > start_y+block_step_y .and. p4(2) < end_y+block_step_y) &
			       .and. (p4(1) > start_x+block_step_x .and. p4(1) < end_x+block_step_x) ))
				
				in_down_left = (( (p1(2) > start_y-block_step_y .and. p1(2) < end_y-block_step_y) &
			       .and. (p1(1) > start_x-block_step_x .and. p1(1) < end_x-block_step_x) ).or. &
			       ( (p2(2) > start_y-block_step_y .and. p2(2) < end_y-block_step_y) &
			       .and. (p2(1) > start_x-block_step_x .and. p2(1) < end_x-block_step_x) ).or. &
			       ( (p3(2) > start_y-block_step_y .and. p3(2) < end_y-block_step_y) &
			       .and. (p3(1) > start_x-block_step_x .and. p3(1) < end_x-block_step_x) ).or. &
			       ( (p4(2) > start_y-block_step_y .and. p4(2) < end_y-block_step_y) &
			       .and. (p4(1) > start_x-block_step_x .and. p4(1) < end_x-block_step_x) ))
				
				in_down_right = ((p1(2) > start_y-block_step_y .and. p1(2) < end_y-block_step_y) &
			       .and. (p1(1) > start_x+block_step_x .and. p1(1) < end_x+block_step_x) ).or.&
			       ( (p2(2) > start_y-block_step_y .and. p2(2) < end_y-block_step_y) &
			       .and. (p2(1) > start_x+block_step_x .and. p2(1) < end_x+block_step_x) ).or.&
			       ( (p3(2) > start_y-block_step_y .and. p3(2) < end_y-block_step_y) &
			       .and. (p3(1) > start_x+block_step_x .and. p3(1) < end_x+block_step_x) ).or.&
			       ( (p4(2) > start_y-block_step_y .and. p4(2) < end_y-block_step_y) &
			       .and. (p4(1) > start_x+block_step_x .and. p4(1) < end_x+block_step_x) )
			       
			       
				if (rank_up /= -1) then
					if(in_up .and. (.not. in_up_right).and.(.not. in_up_left)) then
						COUNTER_overlap(FJ,rank) = COUNTER_overlap(FJ,rank) + 1
						IND_overlap(COUNTER_overlap(FJ,rank),FJ) =  particle_k%ID		
					end if	
					
				end if
			
				if (rank_down /= -1) then
					if (in_down .and. (.not. in_down_left) .and. (.not. in_down_right) ) then
						COUNTER_overlap(BJ,rank) = COUNTER_overlap(BJ,rank) + 1
						IND_overlap(COUNTER_overlap(BJ,rank),BJ) =  particle_k%ID
					end if
					
					
				end if
				if (rank_left /= -1) then
					if (in_left .and.(.not. in_up_left).and.(.not. in_down_left)) then
						COUNTER_overlap(BK,rank) = COUNTER_overlap(BK,rank) + 1
						IND_overlap(COUNTER_overlap(BK,rank),BK) =  particle_k%ID
					end if
				
				
				end if
				if (rank_right /= -1) then
					if (in_right.and.(.not. in_up_right).and.(.not. in_down_right)) then
					 	COUNTER_overlap(FK,rank) = COUNTER_overlap(FK,rank) + 1
						IND_overlap(COUNTER_overlap(FK,rank),FK) =  particle_k%ID
					end if
					
				end if
				if (rank_up_left /= -1) then
					if (in_up_left) then
			        	COUNTER_overlap(FJBK,rank) = COUNTER_overlap(FJBK,rank) + 1
						IND_overlap(COUNTER_overlap(FJBK,rank),FJBK) =  particle_k%ID
					end if
			
				end if
				if (rank_up_right /= -1) then
					if (in_up_right) then
			        	COUNTER_overlap(FJFK,rank) = COUNTER_overlap(FJFK,rank) + 1
						IND_overlap(COUNTER_overlap(FJFK,rank),FJFK) =  particle_k%ID
			        end if
			    
				end if
				if (rank_down_left /= -1) then 
					if (in_down_left) then
						COUNTER_overlap(BJBK,rank) = COUNTER_overlap(BJBK,rank) + 1
						IND_overlap(COUNTER_overlap(BJBK,rank),BJBK) =  particle_k%ID
					end if
				
				end if
				if (rank_down_right /= -1) then
					if (in_down_right) then 
						COUNTER_overlap(BJFK,rank) = COUNTER_overlap(BJFK,rank) + 1
						IND_overlap(COUNTER_overlap(BJFK,rank),BJFK) =  particle_k%ID
					end if
					
				end if
			
			end if 
		
		END SUBROUTINE
		
		
		SUBROUTINE send_overlap_and_danger_particle
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
		COUNTER_recv_danger  = 0
		
		OPP = (/2,1,4,3,6,5,8,7/)		
!==============================================================================
			DO neighloop = 1,8
				DO i=1,COUNTER_overlap(neighloop,rank)
					IF (NEIGHBOUR(neighloop)<0) THEN
						cycle
					end if
					call MPI_SEND(IND_overlap(i,neighloop),1,MPI_INTEGER,NEIGHBOUR(neighloop),101,MPI_COMM_WORLD,code)
				END DO
			
				DO i= 1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
					IF (NEIGHBOUR(neighloop)<0) THEN
						CYCLE
					END IF
					call MPI_RECV(IND_recv_overlap(i,NEIGHBOUR(neighloop)),1,MPI_INTEGER,NEIGHBOUR(neighloop),101, & 
					MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
				END DO
				
				IF (NEIGHBOUR(neighloop)<0) THEN
						cycle
				ELSE	
					COUNTER_recv_overlap = COUNTER_recv_overlap + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
					
				END IF

				DO i=1,COUNTER_danger(neighloop,rank)
					IF (NEIGHBOUR(neighloop)<0) THEN
						cycle
					end if
					call MPI_SEND(IND_danger(i,neighloop),1,MPI_INTEGER,NEIGHBOUR(neighloop),1011,MPI_COMM_WORLD,code)
					
				END DO
			
				DO i= 1,COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
					IF (NEIGHBOUR(neighloop)<0) THEN
						CYCLE
					END IF
					call MPI_RECV(IND_recv_danger(i,NEIGHBOUR(neighloop)),1,MPI_INTEGER,NEIGHBOUR(neighloop),1011, & 
					MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
				END DO
				
				IF (NEIGHBOUR(neighloop)<0) THEN
						cycle
				ELSE	
					COUNTER_recv_danger = COUNTER_recv_danger + COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
					
				END IF
				
			END DO
!==============================================================================
  			call MPI_BARRIER(MPI_COMM_WORLD,code)
			
		END SUBROUTINE send_overlap_and_danger_particle
		
		
		
		SUBROUTINE block_loop_on_block_global_table	
		IMPLICIT NONE
		integer, dimension(NB_NEIGHBOURS) :: OPP
		integer :: i,j,neighloop,Ncum1,Ncum2,Ncum3,Ncum4,Ncum5,Ncum6
		
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
			
		
		Ncum1 = 0
		Ncum2 = 0
		Ncum3 = 0
		Ncum4 = 0
		Ncum5 = 0
		Ncum6 = 0
		
		OPP  = (/2,1,4,3,6,5,8,7/)	
		! UPDATE GLOBAL_TABLE
		! This part's mission is about to include not only IND_inside and IND_overlap
		! but also all possible overlap particles from all possible neighbours 

		! We will have "Npart_block" became: COUNTER_inside + inside sum of COUNTER_overlap in 
		! all direction + COUNTER_recv + COUNTER_danger
		
		! Where COUNTER_danger is the number of particle inside of neighbour block but can 
		! effect on the particles inside the present block
		
		! Where COUNTER_recv_overlap = sum of COUNTER_overlap from neigbour 

		Npart_block = COUNTER_recv_overlap + COUNTER_inside + SUM(COUNTER_overlap(:,rank)) + COUNTER_recv_danger
		
			ALLOCATE(Xpart_block(2,Npart_block))
			ALLOCATE(Mblock(Npart_block))
			ALLOCATE(Dblock(4,Npart_block))
			ALLOCATE(velocity_field(2,Npart_block))
			ALLOCATE(Dout(4,Npart_block))
!==============================================================================
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
		END DO
	
		DO neighloop=1,8
		if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
		DO j=1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
				Xpart_block(1,j+Ncum1) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Xp
				Xpart_block(1,j+Ncum1) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Yp
			END DO
			Ncum1 = Ncum1 + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
		
		END DO
		
		DO neighloop=1,8
		if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
		DO j=1,COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
				Xpart_block(1,j+Ncum1) = ALL_PARTICLES(IND_recv_danger(j,NEIGHBOUR(neighloop)))%Xp
				Xpart_block(1,j+Ncum1) = ALL_PARTICLES(IND_recv_danger(j,NEIGHBOUR(neighloop)))%Yp
			END DO
			Ncum1 = Ncum1 + COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
		END DO

!==============================================================================		
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
		END DO
		DO neighloop=1,8
		if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
		DO j=1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
				Mblock(j+Ncum2) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%mass
			END DO
			Ncum2 = Ncum2 + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))	
		END DO
		DO neighloop=1,8
		if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
		DO j=1,COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
				Mblock(j+Ncum2) = ALL_PARTICLES(IND_recv_danger(j,NEIGHBOUR(neighloop)))%mass
			END DO
			Ncum2 = Ncum2 + COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))	
		END DO
!==============================================================================				
		DO i= 1,COUNTER_inside
			Dblock(1,i) = ALL_PARTICLES(IND_inside(i))%Dp1
			Dblock(2,i) = ALL_PARTICLES(IND_inside(i))%Dp2
			Dblock(3,i) = ALL_PARTICLES(IND_inside(i))%Dp3
			Dblock(4,i) = ALL_PARTICLES(IND_inside(i))%Dp4
		END DO	
		Ncum3 = Ncum3 + COUNTER_inside	
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
		END DO
		DO neighloop=1,8
		if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
		DO j=1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
				Dblock(1,j+Ncum3) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Dp1
				Dblock(2,j+Ncum3) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Dp2
				Dblock(3,j+Ncum3) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Dp3
				Dblock(4,j+Ncum3) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Dp4
			END DO
			Ncum3 = Ncum3 + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))	
		END DO
		DO neighloop=1,8
		if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
		DO j=1,COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
				Dblock(1,j+Ncum3) = ALL_PARTICLES(IND_recv_danger(j,NEIGHBOUR(neighloop)))%Dp1
				Dblock(2,j+Ncum3) = ALL_PARTICLES(IND_recv_danger(j,NEIGHBOUR(neighloop)))%Dp2
				Dblock(3,j+Ncum3) = ALL_PARTICLES(IND_recv_danger(j,NEIGHBOUR(neighloop)))%Dp3
				Dblock(4,j+Ncum3) = ALL_PARTICLES(IND_recv_danger(j,NEIGHBOUR(neighloop)))%Dp4
			END DO
			Ncum3 = Ncum3 + COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
		END DO
!==============================================================================	
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
		END DO	
		DO neighloop=1,8
		if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
		DO j=1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
				velocity_field(1,j+Ncum4) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Upx
				velocity_field(2,j+Ncum4) = ALL_PARTICLES(IND_recv_overlap(j,NEIGHBOUR(neighloop)))%Upy
			END DO
			Ncum4 = Ncum4 + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
		END DO	
		DO neighloop=1,8
		if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
		DO j=1,COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
				velocity_field(1,j+Ncum4) = ALL_PARTICLES(IND_recv_danger(j,NEIGHBOUR(neighloop)))%Upx
				velocity_field(2,j+Ncum4) = ALL_PARTICLES(IND_recv_danger(j,NEIGHBOUR(neighloop)))%Upy
			END DO
			Ncum4 = Ncum4 + COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
		END DO			
!==============================================================================
	
		call diffusion_field_ltp(Xpart_block,Mblock,Dblock,hx,hy,velocity_field,Npart_block)
				
		call update_d_diffusion(Xpart_block, Dblock, Mblock, Tini, dt, time_scheme, hx, hy, &
     indice_max_norm_Dm1, Norm_inf_Dm1, Dout, Npart_block)
    
!==============================================================================		
		DO i=1,COUNTER_inside + SUM(COUNTER_overlap(:,rank))
			Xpart_block(1,i) = Xpart_block(1,i) - velocity_field(1,i)*dt
			Xpart_block(2,i) = Xpart_block(2,i) - velocity_field(2,i)*dt
		END DO
!==============================================================================
		Ncum5 = 0
		DO i= 1,COUNTER_inside
			ALL_PARTICLES(IND_inside(i))%Xp = Xpart_block(1,i)
			ALL_PARTICLES(IND_inside(i))%Yp = Xpart_block(2,i) 
		END DO
		Ncum5 = Ncum5 + COUNTER_inside

		DO neighloop=1,8
			if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
			DO j=1,COUNTER_overlap(neighloop,rank)					
				ALL_PARTICLES(IND_overlap(j,neighloop))%Xp = Xpart_block(1,j+Ncum5)
				ALL_PARTICLES(IND_overlap(j,neighloop))%Yp = Xpart_block(2,j+Ncum5)
			END DO
			Ncum5 = Ncum5 + COUNTER_overlap(neighloop,rank)
		END DO

!==============================================================================
		Ncum6 = 0
		DO i= 1,COUNTER_inside
			ALL_PARTICLES(IND_inside(i))%Dp1 = Dout(1,i)
			ALL_PARTICLES(IND_inside(i))%Dp2 = Dout(2,i)
			ALL_PARTICLES(IND_inside(i))%Dp3 = Dout(3,i)
			ALL_PARTICLES(IND_inside(i))%Dp4 = Dout(4,i)
		END DO
		Ncum6 = Ncum6 + COUNTER_inside
		DO neighloop=1,8
			if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
			DO j=1,COUNTER_overlap(neighloop,rank)
				ALL_PARTICLES(IND_overlap(j,neighloop))%Dp1 = Dout(1,j+Ncum6)
				ALL_PARTICLES(IND_overlap(j,neighloop))%Dp2 = Dout(2,j+Ncum6)
				ALL_PARTICLES(IND_overlap(j,neighloop))%Dp3 = Dout(3,j+Ncum6)
				ALL_PARTICLES(IND_overlap(j,neighloop))%Dp4 = Dout(4,j+Ncum6)
			END DO
			Ncum6 = Ncum6 + COUNTER_overlap(neighloop,rank)			
		END DO
!==============================================================================
			DEALLOCATE(Dout)
			DEALLOCATE(Xpart_block)
			DEALLOCATE(Mblock)
			DEALLOCATE(Dblock)
			DEALLOCATE(velocity_field)	
		
		END SUBROUTINE block_loop_on_block_global_table	
	
		SUBROUTINE update_displacement_on_just_IND_inside_and_IND_overlap
		integer :: neighloop, Ncum1
		integer :: i,j,ID
		integer,dimension(MPI_STATUS_SIZE) :: status_msg
		integer,dimension(1:nb_proc-1)   ::address
		address(:) = 0
		Ncum1 = 0

		call MPI_PART_TYPE
			DO i=1,COUNTER_inside
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
			
			call MPI_BARRIER(MPI_COMM_WORLD,code)
			
			address(:) = 0
			
			DO neighloop=1,8
				if (NEIGHBOUR(neighloop)<0) then
					cycle
				end if
				DO j=1,COUNTER_overlap(neighloop,rank)
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
			call MPI_BARRIER(MPI_COMM_WORLD,code)
		END SUBROUTINE update_displacement_on_just_IND_inside_and_IND_overlap



		SUBROUTINE update_table_block
			IMPLICIT NONE
			LOGICAL :: overlap_check, totally_inside_check, danger_check, leave_check
			logical :: local_overlapcheck, local_insidecheck,local_dangercheck,local_leave_check
			double precision, dimension(2) :: pointu1,pointu2,pointu3,pointu4
			integer :: j,neighloop
			integer,dimension(MPI_STATUS_SIZE) :: status
			integer, dimension(NB_NEIGHBOURS) :: OPP	
		
			double precision :: axe
			
			COUNTER_recv_overlap = 0
			OPP = (/2,1,4,3,6,5,8,7/)
			
			DO i=1,COUNTER_inside
				call overlap_criterion(ALL_PARTICLES(IND_inside(i)),overlap_check,totally_inside_check,danger_check, & 
				pointu1,pointu2,pointu3,pointu4,axe)
				
					ALL_PARTICLES(IND_inside(i))%pointu1 = pointu1
					ALL_PARTICLES(IND_inside(i))%pointu2 = pointu2
					ALL_PARTICLES(IND_inside(i))%pointu3 = pointu3
					ALL_PARTICLES(IND_inside(i))%pointu4 = pointu4
					ALL_PARTICLES(IND_inside(i))%axe = axe
				if (totally_inside_check) then
					leave_check = .false.
					call particle_screening(ALL_PARTICLES(IND_inside(i)),overlap_check,danger_check,leave_check)
				else
					if(overlap_check) then
						leave_check = .false.
						call particle_screening(ALL_PARTICLES(IND_inside(i)),overlap_check,danger_check,leave_check)
						IND_inside(i:COUNTER_inside-1) = IND_inside(i+1:COUNTER_inside)
						COUNTER_inside = COUNTER_inside -1
					else
						leave_check = .true.
						call particle_screening(ALL_PARTICLES(IND_inside(i)),overlap_check,danger_check,leave_check)
						IND_inside(i:COUNTER_inside-1) = IND_inside(i+1:COUNTER_inside)
						COUNTER_inside = COUNTER_inside -1
					end if
				end if
			END DO		
			
			DO neighloop=1,8
			if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
			DO j=1,COUNTER_overlap(neighloop,rank)					
				call overlap_criterion(ALL_PARTICLES(IND_overlap(j,neighloop)),overlap_check,totally_inside_check,danger_check, & 
				pointu1,pointu2,pointu3,pointu4,axe)
				ALL_PARTICLES(IND_overlap(j,neighloop))%pointu1 = pointu1
				ALL_PARTICLES(IND_overlap(j,neighloop))%pointu2 = pointu2
				ALL_PARTICLES(IND_overlap(j,neighloop))%pointu3 = pointu3
				ALL_PARTICLES(IND_overlap(j,neighloop))%pointu4 = pointu4
				ALL_PARTICLES(IND_overlap(j,neighloop))%axe = axe
				if (totally_inside_check) then
					leave_check = .false.
					call particle_screening(ALL_PARTICLES(IND_overlap(j,neighloop)),overlap_check,danger_check,leave_check)
					COUNTER_inside = COUNTER_inside + 1
					IND_inside(COUNTER_inside) = IND_overlap(j,neighloop)
					IND_overlap(j:COUNTER_overlap(neighloop,rank)-1,neighloop) &
						 = IND_overlap(j+1:COUNTER_overlap(neighloop,rank),neighloop)
						COUNTER_overlap(neighloop,rank) = COUNTER_overlap(neighloop,rank) -1
				else
					if(overlap_check) then
						cycle
					else
						leave_check = .true.
						call particle_screening(ALL_PARTICLES(IND_overlap(j,neighloop)),overlap_check,danger_check,leave_check)
						IND_overlap(j:COUNTER_overlap(neighloop,rank)-1,neighloop) = & 
						IND_overlap(j+1:COUNTER_overlap(neighloop,rank),neighloop)
						COUNTER_overlap(neighloop,rank) = COUNTER_overlap(neighloop,rank) -1
					end if
				end if
			END DO
			END DO
				
			DO neighloop = 1,8
				DO i=1,COUNTER_leave(neighloop,rank)
					if(NEIGHBOUR(neighloop)<0) then
						cycle
					end if
					call MPI_SEND(IND_leave(i,neighloop),1,MPI_INTEGER,NEIGHBOUR(neighloop),8,MPI_COMM_WORLD,code)
				END DO
				
				DO i=1,COUNTER_leave(OPP(neighloop),NEIGHBOUR(neighloop))
					if(NEIGHBOUR(neighloop)<0) then
						cycle
					end if
					call MPI_RECV(IND_recv_leave(i,NEIGHBOUR(neighloop)),1,MPI_INTEGER,NEIGHBOUR(neighloop),8 & 
					,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
				END DO
				if (NEIGHBOUR(neighloop)<0) then
					cycle
				else
					COUNTER_recv_leave = COUNTER_recv_leave + COUNTER_leave(OPP(neighloop),NEIGHBOUR(neighloop))
				end if	
			END DO

			DO neighloop=1,8
				if (NEIGHBOUR(neighloop)<0) then
					cycle
				end if
				DO j=1,COUNTER_leave(OPP(neighloop),NEIGHBOUR(neighloop))
					
					call overlap_criterion(ALL_PARTICLES(IND_recv_leave(j,NEIGHBOUR(neighloop))),local_overlapcheck, & 
					local_insidecheck,local_dangercheck,pointu1,pointu2,pointu3,pointu4,axe)	
					
					if (local_insidecheck) then
						COUNTER_inside = COUNTER_inside + 1
						IND_inside(COUNTER_inside) = i	
						call particle_screening(ALL_PARTICLES(IND_recv_leave(j,NEIGHBOUR(neighloop))), & 
						local_overlapcheck,local_dangercheck,local_leave_check)
					else 
						call particle_screening(ALL_PARTICLES(IND_recv_leave(j,NEIGHBOUR(neighloop))), & 
						local_overlapcheck,local_dangercheck,local_leave_check)
				end if
					
				END DO
		
		END DO
				
		END SUBROUTINE update_table_block
		

		SUBROUTINE update_all_particle_information
		
		if(rank==0) then
			DO i=1,number_of_particles

				Xparticle_read(1,i) = ALL_PARTICLES(i)%Xp
				Xparticle_read(2,i) = ALL_PARTICLES(i)%Yp		

				D_read(1,i) = ALL_PARTICLES(i)%Dp1
				D_read(2,i) = ALL_PARTICLES(i)%Dp2
				D_read(3,i) = ALL_PARTICLES(i)%Dp3 
				D_read(4,i) = ALL_PARTICLES(i)%Dp4
			END DO

			OPEN(unit = 2, FILE = 'trunk/fortran_srcs/coords4fortran.bin',FORM="UNFORMATTED", &
 STATUS="UNKNOWN",POSITION="REWIND", ACTION="READWRITE", ACCESS='STREAM')
			WRITE(2) Xparticle_read
			CLOSE(2)
			
			OPEN(unit= 3, FILE='trunk/fortran_srcs/deformmatrix4fortran.bin',FORM="UNFORMATTED",STATUS="UNKNOWN",POSITION="REWIND" & 
,ACTION="READWRITE",ACCESS='STREAM')			
			
			write(3) D_read
			close(3)
			
		end if
		
		END SUBROUTINE update_all_particle_information



		SUBROUTINE dealloc_all_particle_table
			DEALLOCATE(IND_leave)
			DEALLOCATE(IND_inside)
			DEALLOCATE(IND_overlap)
			DEALLOCATE(ALL_PARTICLES)
			
			DEALLOCATE(COUNTER_leave)
			DEALLOCATE(COUNTER_overlap)
			DEALLOCATE(IND_recv_overlap)
		
			
		END SUBROUTINE dealloc_all_particle_table

END MODULE mpi_2d_structures_modf90
