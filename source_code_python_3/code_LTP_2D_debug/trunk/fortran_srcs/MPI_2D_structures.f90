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
	USE calculsfor_ini_modf90
	USE data_launch
	USE Jacobi_method
	
	IMPLICIT NONE
	
	! DECALARATIONS
	
	! NB_NEIGHBOURS is already declared in pack
	
	INTEGER, DIMENSION(NB_NEIGHBOURS)                     :: NEIGHBOUR
	INTEGER, PARAMETER                                    :: FJ=1,BJ=2,FK=3,BK=4 ! neighbours : Forward_J, Backward_J, ... 
  	INTEGER, PARAMETER                                    :: FJBK=5,BJFK=6,BJBK=7,FJFK=8 ! Common edge neighbour
	
	
	! IND_inside indicates the identity of particles inside the present block
	INTEGER, DIMENSION(:), ALLOCATABLE                    :: IND_inside
	
	! IND_overlap indicates the identity of overlap particles in each direction
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: IND_overlap
	
	! IND_overlap indicates the identity of danger particles in each direction
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: IND_danger
	
	! IND_leave indicates the identity(in ALL_PARTICLES) of particles will leave the present block
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: IND_leave
	
	! IND_recv_danger indicates the identity of danger particles receving from neighbours
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: IND_recv_danger
	
	! IND_recv_danger indicates the identity of leaving particles receving from neighbours
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: IND_recv_leave
	
	! IND_recv_overlap indicates the identity of particles receving from neighbours
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: IND_recv_overlap
	
	
	
	
	! COUNTER_inside indicates the number of particles inside the present block
	INTEGER, DIMENSION(:), ALLOCATABLE                   :: COUNTER_inside 
	
	! COUNTER_recv_danger indicates the number of overlap particles receiving in each direction 

	INTEGER                                               :: COUNTER_recv_danger
	
	! COUNTER_recv_overlap indicates the number of overlap particles receiving in each direction 
	INTEGER                                               :: COUNTER_recv_overlap
	
	! COUNTER_recv_leave indicates the number of overlap particles leaving in each direction
	INTEGER                                               :: COUNTER_recv_leave
	
	! COUNTER_leave indicates the number of particles leaving in each direction
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: COUNTER_leave
	
	! COUNTER_overlap indicates the number of overlap particles in each direction
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: COUNTER_overlap
	
	! COUNTER_overlap indicates the number of danger particles in each direction
	INTEGER, DIMENSION(:,:), ALLOCATABLE                  :: COUNTER_danger
	
	
	
	
	
	
	! ALL_PARTICLES stocks all particle informations and must be broadcasted after being updated all particles displacement.
	TYPE(PARTTYPE), DIMENSION(:), ALLOCATABLE             :: ALL_PARTICLES
	
	! neighbour_limit give neighbour boundaries informations, each block has 8 neighbours, each neighbour must have 4 real values to identify its boundaries 
	DOUBLE PRECISION, DIMENSION(4,8)                      :: neighbour_limit



	

	
	CONTAINS
!================================================================================================	
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

!================================================================================================
		
		SUBROUTINE initiation_table
		
			ALLOCATE(IND_inside(number_of_particles))
			ALLOCATE(IND_leave(number_of_particles,NB_NEIGHBOURS))
			ALLOCATE(IND_overlap(number_of_particles,NB_NEIGHBOURS))
			ALLOCATE(IND_danger(number_of_particles,NB_NEIGHBOURS))	
			
			ALLOCATE(IND_recv_overlap(number_of_particles,NB_NEIGHBOURS))
			ALLOCATE(IND_recv_danger(number_of_particles,NB_NEIGHBOURS))
			ALLOCATE(IND_recv_leave(number_of_particles,NB_NEIGHBOURS))
			
			ALLOCATE(COUNTER_inside(0:nb_proc-1))
			ALLOCATE(COUNTER_danger(NB_NEIGHBOURS,0:nb_proc-1))
			ALLOCATE(COUNTER_overlap(NB_NEIGHBOURS,0:nb_proc-1))
			ALLOCATE(COUNTER_leave(NB_NEIGHBOURS,0:nb_proc-1))
			
			ALLOCATE(ALL_PARTICLES(number_of_particles))
								
		END SUBROUTINE initiation_table
!================================================================================================
		
		SUBROUTINE parttype_convert
			DO i=1,number_of_particles
				ALL_PARTICLES(i)%Xp = Xparticle_read(1,i)
				ALL_PARTICLES(i)%Yp = Xparticle_read(2,i)
				ALL_PARTICLES(i)%Dp1 = D_read(1,i)
				ALL_PARTICLES(i)%Dp2 = D_read(2,i)
				ALL_PARTICLES(i)%Dp3 = D_read(3,i)
				ALL_PARTICLES(i)%Dp4 = D_read(4,i)
				ALL_PARTICLES(i)%ID = i
				ALL_PARTICLES(i)%mass = mass_read(i) 
			END DO
		
		END SUBROUTINE parttype_convert
		
!================================================================================================




		!----------------------------------------!
		! OVERLAP CRITERION FOR JUST ONE PARTICLE!
		!----------------------------------------!
			
		SUBROUTINE overlap_criterion(particle_k,overlap_check,inside_check,danger_check, & 
		pointu1,pointu2,pointu3,pointu4,axe)
		
		! Retrieve the eigenvalues and eigenvectors to determine the direction movement of particles, 	then decide its status
		
		! Input: k-th_particle in "PARTYPE" type 
		! Output: logical state
		
		
			IMPLICIT NONE
			TYPE(PARTTYPE), INTENT(in)                      :: particle_k
			LOGICAL, INTENT(out)                            :: overlap_check
			LOGICAL, INTENT(out)                            :: inside_check
			LOGICAL, INTENT(out)                            :: danger_check
			DOUBLE PRECISION, INTENT(out)                   :: axe
			DOUBLE PRECISION, dimension(2), intent(out)     :: pointu1, pointu2, & 
															   pointu3,pointu4
																	   
			LOGICAL                                         :: partial_inside_check			
			DOUBLE PRECISION                                :: a(2,2), x(2,2)
			DOUBLE PRECISION, PARAMETER                     :: abserr=1e-6
			INTEGER i, j


			double precision, dimension(2)             	    :: eigenvector_1
			double precision, dimension(2)                  :: eigenvector_2
			DOUBLE PRECISION                                :: eigenvalue_1
			DOUBLE PRECISION                                :: eigenvalue_2
			

			DOUBLE PRECISION		                        :: D11,D12,D21,D22, d1,d2,d3,d4
			
			
			LOGICAL                                         :: pointu1inside,pointu2inside, & 
					                                           pointu3inside,pointu4inside
			

				
			D11 = particle_k%Dp1
			D12 = particle_k%Dp2
			D21 = particle_k%Dp3
			D22 = particle_k%Dp4

			! matrix D
  			a(1,1) = D11
  			a(1,2) = D12
  			a(2,1) = D21
  			a(2,2) = D22
		
			

			call Jacobi(a,x,abserr,2)

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
			
			pointu1inside = (pointu1(1)-start_x>=0) &
				.and.(pointu1(1)-END_x<=0)          &
				.and.(pointu1(2)-start_y>=0)        &
				.and.(pointu2(2)-END_y<=0)
			pointu2inside = (pointu2(1)-start_x>=0) &
				.and.(pointu2(1)-END_x<=0)          &
				.and.(pointu2(2)-start_y>=0)        &
				.and.(pointu2(2)-END_y<=0)
			pointu3inside = (pointu3(1)-start_x>=0) &
				.and.(pointu3(1)-END_x<=0)          &
				.and.(pointu3(2)-start_y>=0)        &
				.and.(pointu3(2)-END_y<=0)
			pointu4inside = (pointu4(1)-start_x>=0) &
				.and.(pointu4(1)-END_x<=0)          &
				.and.(pointu4(2)-start_y>=0)        &
				.and.(pointu4(2)-END_y<=0)		
			

		IF (particle_k%Xp>=start_x .and. particle_k%Xp<end_x &
		.and. particle_k%Yp>=start_y .and. particle_k%Yp<end_y) THEN
		
			inside_check = .true.
			
			IF(pointu1inside.and.pointu2inside.and.pointu3inside.and.pointu4inside)then
			
				overlap_check = .false.
				d1 = particle_k%Xp - start_x
				d2 = -particle_k%Xp + end_x
				d3 = particle_k%Yp - start_y
				d4 = -particle_k%Yp + end_y
				
		  		IF (d1<2*axe .or. d2<2*axe .or. d3<2*axe .or. d4 <2*axe) THEN
					danger_check = .true.
				ELSE
					danger_check = .false.
				END IF
			
			ELSE
			
				overlap_check = .true.
				danger_check  = .false.
			
			END IF
			
		
		ELSE 
			overlap_check = .false.
			inside_check  = .false.
			danger_check  = .false.
			
		END IF
		
		END SUBROUTINE overlap_criterion
		
!================================================================================================	
		
		SUBROUTINE particle_distribution
			IMPLICIT NONE
			logical                           :: local_overlapcheck, local_insidecheck, & 
												 local_dangercheck,local_leave_check
			DOUBLE PRECISION, dimension(2)    :: pointu1, pointu2,pointu3,pointu4
			DOUBLE PRECISION                  :: axe
			integer                           :: ID
			
			COUNTER_inside(:)         = 0
			COUNTER_leave(:,rank)     = 0
			COUNTER_overlap(:,rank)   = 0
			COUNTER_danger(:,rank)    = 0

			! AT THE BEGINNING: THERE IS NO PARTICLE LEAVING, WE DISTRIBUTE PARTICLES FOR EVERY BLOCK.
			local_leave_check = .false.

!================================================================================================

			DO i =1,number_of_particles
				CALL overlap_criterion(ALL_PARTICLES(i),local_overlapcheck,local_insidecheck,local_dangercheck & 
				,pointu1,pointu2,pointu3,pointu4,axe)
!				IF(rank==8) THEN
!					
!					print*,ALL_PARTICLES(i)%ID,local_overlapcheck,local_insidecheck,local_dangercheck
!					
!				END IF

				ALL_PARTICLES(i)%pointu1 = pointu1
				ALL_PARTICLES(i)%pointu2 = pointu2
				ALL_PARTICLES(i)%pointu3 = pointu3
				ALL_PARTICLES(i)%pointu4 = pointu4
				ALL_PARTICLES(i)%axe     = axe

				IF (local_insidecheck) then
					COUNTER_inside(rank)         = COUNTER_inside(rank) + 1
					IND_inside(COUNTER_inside(rank)) = ALL_PARTICLES(i)%ID
					

					call particle_screening(ALL_PARTICLES(i),local_overlapcheck,local_dangercheck,local_leave_check)

				END IF
			END DO			
			
!================================================================================================  		
			IF(rank /= 0) then
      			call MPI_SEND(COUNTER_inside(rank),1,MPI_INTEGER,0,1,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1
       				call MPI_RECV(COUNTER_inside(ID),1,MPI_INTEGER,ID,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			END do	
    		END IF
    		call MPI_BCAST(COUNTER_inside,nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)
					
!================================================================================================  		
			IF(rank /= 0) then
      			call MPI_SEND(COUNTER_overlap(1,rank),NB_NEIGHBOURS,MPI_INTEGER,0,1,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1
       				call MPI_RECV(COUNTER_overlap(1,ID),NB_NEIGHBOURS,MPI_INTEGER,ID,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			END do	
    		END IF
    		call MPI_BCAST(COUNTER_overlap,NB_NEIGHBOURS*nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)
!================================================================================================  		
    		IF(rank /= 0) then
      			call MPI_SEND(COUNTER_leave(1,rank),NB_NEIGHBOURS,MPI_INTEGER,0,2,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1
       				call MPI_RECV(COUNTER_leave(1,ID),NB_NEIGHBOURS,MPI_INTEGER,ID,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			END do
    		END IF
    		call MPI_BCAST(COUNTER_leave,NB_NEIGHBOURS*nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)
    		
!================================================================================================    		
    		IF(rank /= 0) then
      			call MPI_SEND(COUNTER_danger(1,rank),NB_NEIGHBOURS,MPI_INTEGER,0,3,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1
       				call MPI_RECV(COUNTER_danger(1,ID),NB_NEIGHBOURS,MPI_INTEGER,ID,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			END do
    		END IF
    		call MPI_BCAST(COUNTER_danger,NB_NEIGHBOURS*nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)
    			
		END SUBROUTINE particle_distribution

!================================================================================================

SUBROUTINE particle_distribution_v2
	! Every block scan its own inside particles in order to obtain 3 tables, "totally_inside_particles" and 
	! "overlap_particels", "danger_particles"
	integer :: step_x,step_y,number_of_rows,number_of_columns,ind,j

	logical                           :: local_overlapcheck, local_insidecheck, & 
												 local_dangercheck,local_leave_check
			DOUBLE PRECISION, dimension(2)    :: pointu1, pointu2,pointu3,pointu4
			DOUBLE PRECISION                  :: axe
			integer                           :: ID
			
			COUNTER_inside(:)         = 0
			COUNTER_leave(:,rank)     = 0
			COUNTER_overlap(:,rank)   = 0
			COUNTER_danger(:,rank)    = 0

			! AT THE BEGINNING: THERE IS NO PARTICLE LEAVING, WE DISTRIBUTE PARTICLES FOR
			! EVERY BLOCK.
			local_leave_check = .false.

!================================================================================================
IF((coords(1)>=1 .and. coords(1)<=dims(1)-2) .and.(coords(2)>=1 .and. coords(2)<=dims(2)-2)) then
	step_x = nxg/(dims(1)-2)
	step_y = number_of_particles/(dims(2)-2)
	number_of_rows = nyg/(dims(2)-2)
	number_of_columns = nxg/(dims(1)-2)


	DO j=1,number_of_rows
		DO i=1,number_of_columns
			ind = step_y*(dims(2)-2-coords(2))+step_x*(coords(1)-1)+i+(j-1)*nxg
			CALL overlap_criterion(ALL_PARTICLES(ind),local_overlapcheck,local_insidecheck,local_dangercheck & 
				,pointu1,pointu2,pointu3,pointu4,axe)
			
				ALL_PARTICLES(ind)%pointu1 = pointu1
				ALL_PARTICLES(ind)%pointu2 = pointu2
				ALL_PARTICLES(ind)%pointu3 = pointu3
				ALL_PARTICLES(ind)%pointu4 = pointu4
				ALL_PARTICLES(ind)%axe     = axe

				IF (local_insidecheck) then
					COUNTER_inside(rank)         = COUNTER_inside(rank) + 1
					IND_inside(COUNTER_inside(rank)) = ind
					
					call particle_screening(ALL_PARTICLES(ind),local_overlapcheck,local_dangercheck,local_leave_check)
					if(rank==8)then
						print*,ind,local_overlapcheck,local_dangercheck,local_leave_check
					end if
				END IF
		END DO
	END DO
END IF

!================================================================================================  		
			IF(rank /= 0) then
      			call MPI_SEND(COUNTER_inside(rank),1,MPI_INTEGER,0,1,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1
       				call MPI_RECV(COUNTER_inside(ID),1,MPI_INTEGER,ID,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			END do	
    		END IF
    		call MPI_BCAST(COUNTER_inside,nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)
					
!================================================================================================  		
			IF(rank /= 0) then
      			call MPI_SEND(COUNTER_overlap(1,rank),NB_NEIGHBOURS,MPI_INTEGER,0,1,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1
       				call MPI_RECV(COUNTER_overlap(1,ID),NB_NEIGHBOURS,MPI_INTEGER,ID,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			END do	
    		END IF
    		call MPI_BCAST(COUNTER_overlap,NB_NEIGHBOURS*nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)
!================================================================================================  		
    		IF(rank /= 0) then
      			call MPI_SEND(COUNTER_leave(1,rank),NB_NEIGHBOURS,MPI_INTEGER,0,2,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1
       				call MPI_RECV(COUNTER_leave(1,ID),NB_NEIGHBOURS,MPI_INTEGER,ID,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			END do
    		END IF
    		call MPI_BCAST(COUNTER_leave,NB_NEIGHBOURS*nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)
    		
!================================================================================================    		
    		IF(rank /= 0) then
      			call MPI_SEND(COUNTER_danger(1,rank),NB_NEIGHBOURS,MPI_INTEGER,0,3,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1
       				call MPI_RECV(COUNTER_danger(1,ID),NB_NEIGHBOURS,MPI_INTEGER,ID,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			END do
    		END IF
    		call MPI_BCAST(COUNTER_danger,NB_NEIGHBOURS*nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)


		END SUBROUTINE particle_distribution_v2

!================================================================================================
	
		SUBROUTINE neighbour_limit_finding

			integer                           :: j,neighloop,ID

			integer, dimension(NB_NEIGHBOURS) :: OPP	
			double precision, dimension(4)    :: limit
			
			OPP = (/2,1,4,3,6,5,8,7/)

			
			limit(1) = start_x
			limit(2) = END_x
			limit(3) = start_y
			limit(4) = END_y
			DO neighloop=1,8
				IF (NEIGHBOUR(neighloop)<0)then
					cycle
				END IF
				call MPI_SEND(limit,4,MPI_DOUBLE_PRECISION,NEIGHBOUR(neighloop),111,MPI_COMM_WORLD,code)
				
			END DO			
			
			DO neighloop=1,8
				IF(NEIGHBOUR(neighloop)<0) then
					cycle
				END IF
				call MPI_RECV(neighbour_limit(1,neighloop),4,MPI_DOUBLE_PRECISION,NEIGHBOUR(neighloop),111,MPI_COMM_WORLD & 
				,MPI_STATUS_IGNORE,code)
			END DO
		END SUBROUTINE neighbour_limit_finding

	
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
			
		
			IF (leave) then
				
				xk = particle_k%Xp
				yk = particle_k%Yp
				
!================================================================================================


				in_up =((yk > neighbour_limit(3,1) .and. yk < neighbour_limit(4,1)) & 
					.and. (xk > neighbour_limit(1,1) .and. xk < neighbour_limit(2,1)))
					
				in_down = ((yk >  neighbour_limit(3,2) .and. yk < neighbour_limit(4,2)) &
					.and. (xk > neighbour_limit(1,2)  .and. xk < neighbour_limit(2,2)))
				
				in_left = ((yk > neighbour_limit(3,4) .and. yk < neighbour_limit(4,4)) &
				    .and. (xk > neighbour_limit(1,4)  .and. xk< neighbour_limit(2,4)))
				
				in_right = ((yk >  neighbour_limit(3,3) .and. yk < neighbour_limit(4,3)) & 
					.and. (xk > neighbour_limit(1,3)  .and. xk < neighbour_limit(2,3))) 
				
				in_up_left = ((yk >  neighbour_limit(3,5) .and. yk <neighbour_limit(4,5)) &
			       .and. (xk > neighbour_limit(1,5)  .and. xk < neighbour_limit(2,5)))
			       
				in_up_right = ((yk >  neighbour_limit(3,8) .and. yk < neighbour_limit(4,8)) &
			       .and. (xk >neighbour_limit(1,8)  .and. xk < neighbour_limit(2,8)))
				
				in_down_left = ((yk >  neighbour_limit(3,7) .and. yk < neighbour_limit(4,7)) &
			       .and. (xk > neighbour_limit(1,7) .and. xk < neighbour_limit(2,7)))
				
				in_down_right = ((yk >  neighbour_limit(3,6) .and. yk < neighbour_limit(4,6)) &
			       .and. (xk > neighbour_limit(1,6)  .and. xk < neighbour_limit(2,6)))

			     
		 		
!================================================================================================
			    
			  	IF (rank_up /= -1) then
			    	IF(in_up) then
			    		COUNTER_leave(FJ,rank) = COUNTER_leave(FJ,rank) + 1
			    		IND_leave(COUNTER_leave(FJ,rank),FJ) = particle_k%ID
			    	END IF
			 	END IF
!================================================================================================
			  	IF (rank_down /= -1) then
			    	IF(in_down) then
			    		COUNTER_leave(BJ,rank) = COUNTER_leave(BJ,rank) + 1
			    		IND_leave(COUNTER_leave(BJ,rank),BJ) = particle_k%ID

			    	END IF
			 	END IF
!================================================================================================
			  	IF (rank_left /= -1) then  
			    	IF(in_left) then
			    		COUNTER_leave(BK,rank) = COUNTER_leave(BK,rank) + 1
			    		IND_leave(COUNTER_leave(BK,rank),BK) = particle_k%ID
			    	END IF
				END IF
!================================================================================================
			   IF (rank_right /= -1) then  
			    	IF(in_right) then
			    		COUNTER_leave(FK,rank) = COUNTER_leave(FK,rank) + 1
			    		IND_leave(COUNTER_leave(FK,rank),FK) = particle_k%ID
			    	END IF
			   END IF

!================================================================================================			    
			    IF (rank_up_left /= -1) then
			    	IF(in_up_left) then
			    		COUNTER_leave(FJBK,rank) =COUNTER_leave(FJBK,rank) + 1
			    		IND_leave(COUNTER_leave(FJBK,rank),FJBK) = particle_k%ID
			    	END IF
			 	END IF
!================================================================================================			    
				IF (rank_up_right /= -1) then
			    	IF(in_up_right) then
			    		COUNTER_leave(FJFK,rank) = COUNTER_leave(FJFK,rank) +1
			    		IND_leave(COUNTER_leave(FJFK,rank),FJFK) = particle_k%ID
			    	END IF
				END IF
!================================================================================================
			    IF (rank_down_left /= -1) then
			    	IF(in_down_left) then
			    		COUNTER_leave(BJBK,rank) = COUNTER_leave(BJBK,rank) + 1
			    		IND_leave(COUNTER_leave(BJBK,rank),BJBK) = particle_k%ID
			    			
			    	END IF 
				END IF
!================================================================================================   
				IF (rank_down_right /= -1) then
					IF(in_down_right) then
						COUNTER_leave(BJFK,rank) = COUNTER_leave(BJFK,rank) + 1
						IND_leave(COUNTER_leave(BJFK,rank),BJFK) = particle_k%ID
					END IF
				END IF
				
			END IF

!================================================================================================
			IF(danger) then
				
				d1 = particle_k%Xp - start_x
				d2 = -particle_k%Xp + END_x
				d3 = particle_k%Yp - start_y
				d4 = -particle_k%Yp + END_y
				axe_k = particle_k%axe
				
				
				IF (rank_up /= -1) then
				IF (d4<2*axe_k ) then 
						COUNTER_danger(FJ,rank) = COUNTER_danger(FJ,rank) + 1
						IND_danger(COUNTER_danger(FJ,rank),FJ) =  particle_k%ID
					END IF
				END IF
!================================================================================================				
				IF (rank_down /= -1) then
				IF (d3<2*axe_k) then
						COUNTER_danger(BJ,rank) = COUNTER_danger(BJ,rank) + 1
						IND_danger(COUNTER_danger(BJ,rank),BJ) =  particle_k%ID
					END IF
				END IF
!================================================================================================				
				IF (rank_left /= -1) then
					IF(d1<2*axe_k) then 
						COUNTER_danger(BK,rank) = COUNTER_danger(BK,rank) + 1
						IND_danger(COUNTER_danger(BK,rank),BK) =  particle_k%ID
					END IF
				END IF
!================================================================================================				
				IF (rank_right /= -1) then	
					IF(d2<2*axe_k) then 
						COUNTER_danger(FK,rank) = COUNTER_danger(FK,rank) + 1
						IND_danger(COUNTER_danger(FK,rank),FK) =  particle_k%ID
					END IF
				END IF
!================================================================================================				
				IF (rank_up_left /= -1) then		
					IF(d4<2*axe_k .and.d1<2*axe_k) then
						COUNTER_danger(FJBK,rank) = COUNTER_danger(FJBK,rank) + 1
						IND_danger(COUNTER_danger(FJBK,rank),FJBK) =  particle_k%ID
					END IF
				END IF
!================================================================================================				
				IF (rank_up_right /= -1) then	    
			        IF (d4<2*axe_k .and. d2<2*axe_k) then
			        	COUNTER_danger(FJFK,rank) = COUNTER_danger(FJFK,rank) + 1
						IND_danger(COUNTER_danger(FJFK,rank),FJFK) =  particle_k%ID
			        END IF
			    END IF
!================================================================================================			    
				IF (rank_down_left /= -1) then	
					IF (d3<2*axe_k .and. d1<2*axe_k) then
						COUNTER_danger(BJBK,rank) = COUNTER_danger(BJBK,rank) + 1
						IND_danger(COUNTER_danger(BJBK,rank),BJBK) =  particle_k%ID
					END IF
				END IF
!================================================================================================					
				IF (rank_down_right /= -1) then
					IF (d3<2*axe_k .and. d2<2*axe_k) then
						COUNTER_danger(BJFK,rank) = COUNTER_danger(BJFK,rank) + 1
						IND_danger(COUNTER_danger(BJFK,rank),BJFK) =  particle_k%ID
					END IF
				END IF
			END IF
!================================================================================================				
			IF(overlap) then
				
				p1 = particle_k%pointu1
				p2 = particle_k%pointu2
				p3 = particle_k%pointu3
				p4 = particle_k%pointu4
				
				in_up =((p1(2) > neighbour_limit(3,1) .and. p1(2) < neighbour_limit(4,1)) & 
					.and. (p1(1) > neighbour_limit(1,1) .and. p1(1) < neighbour_limit(2,1)) .or. &
						(p2(2) > neighbour_limit(3,1) .and. p2(2) < neighbour_limit(4,1)) & 
					.and. (p2(1) > neighbour_limit(1,1) .and. p3(1) < neighbour_limit(2,1)) .or. &
						(p3(2) > neighbour_limit(3,1) .and. p3(2) < neighbour_limit(4,1)) & 
					.and. (p3(1) > neighbour_limit(1,1) .and. p3(1) < neighbour_limit(2,1)) .or. &
						(p4(2) > neighbour_limit(3,1) .and. p4(2) < neighbour_limit(4,1)) & 
					.and. (p4(1) > neighbour_limit(1,1) .and. p4(1) < neighbour_limit(2,1)))
					
				in_down = ((p1(2) >  neighbour_limit(3,2) .and. p1(2) < neighbour_limit(4,2)) &
					.and. (p1(1) > neighbour_limit(1,2)  .and. p1(1) < neighbour_limit(2,2)) .or.&
						   (p2(2) >  neighbour_limit(3,2) .and. p2(2) < neighbour_limit(4,2)) &
					.and. (p2(1) > neighbour_limit(1,2)  .and. p2(1) < neighbour_limit(2,2)) .or.&
						   (p3(2) >  neighbour_limit(3,2) .and. p3(2) < neighbour_limit(4,2)) &
					.and. (p3(1) > neighbour_limit(1,2)  .and. p3(1) < neighbour_limit(2,2)) .or.&
						   (p4(2) >  neighbour_limit(3,2) .and. p4(2) < neighbour_limit(4,2)) &
					.and. (p4(1) > neighbour_limit(1,2)  .and. p4(1) < neighbour_limit(2,2)))
				
				in_left = ((p1(2) > neighbour_limit(3,4) .and. p1(2) < neighbour_limit(4,4)) &
				    .and. (p1(1) > neighbour_limit(1,4)  .and. p1(1)< neighbour_limit(2,4)) .or.&
				    	   (p2(2) > neighbour_limit(3,4) .and. p2(2) < neighbour_limit(4,4)) &
				    .and. (p2(1) > neighbour_limit(1,4)  .and. p2(1)< neighbour_limit(2,4)) .or.&
				    	   (p3(2) > neighbour_limit(3,4) .and. p3(2) < neighbour_limit(4,4)) &
				    .and. (p3(1) > neighbour_limit(1,4)  .and. p3(1)< neighbour_limit(2,4)) .or.&
				    	   (p4(2) > neighbour_limit(3,4) .and. p4(2) < neighbour_limit(4,4)) &
				    .and. (p4(1) > neighbour_limit(1,4)  .and. p4(1)< neighbour_limit(2,4)))
				
				in_right = ((p1(2) >  neighbour_limit(3,3) .and. p1(2) < neighbour_limit(4,3)) & 
					.and. (p1(1) > neighbour_limit(1,3)  .and. p1(1) < neighbour_limit(2,3)) .or. &
							(p2(2) >  neighbour_limit(3,3) .and. p2(2) < neighbour_limit(4,3)) & 
					.and. (p2(1) > neighbour_limit(1,3)  .and. p2(1) < neighbour_limit(2,3)) .or. &
							(p3(2) >  neighbour_limit(3,3) .and. p3(2) < neighbour_limit(4,3)) & 
					.and. (p3(1) > neighbour_limit(1,3)  .and. p3(1) < neighbour_limit(2,3)) .or.&
							(p4(2) >  neighbour_limit(3,3) .and. p4(2) < neighbour_limit(4,3)) & 
					.and. (p4(1) > neighbour_limit(1,3)  .and. p4(1) < neighbour_limit(2,3))) 
				
				in_up_left = ((p1(2) >  neighbour_limit(3,5) .and. p1(2) <neighbour_limit(4,5)) &
			       .and. (p1(1) > neighbour_limit(1,5)  .and. p1(1) < neighbour_limit(2,5)) .or.&
			       			  (p2(2) >  neighbour_limit(3,5) .and. p2(2) <neighbour_limit(4,5)) &
			       .and. (p2(1) > neighbour_limit(1,5)  .and. p2(1) < neighbour_limit(2,5)) .or.&
			       		  	  (p3(2) >  neighbour_limit(3,5) .and. p3(2) <neighbour_limit(4,5)) &
			       .and. (p3(1) > neighbour_limit(1,5)  .and. p3(1) < neighbour_limit(2,5)) .or. &
			       			  (p4(2) >  neighbour_limit(3,5) .and. p4(2) <neighbour_limit(4,5)) &
			       .and. (p4(1) > neighbour_limit(1,5)  .and. p4(1) < neighbour_limit(2,5)))
			       
				in_up_right = ((p1(2) >  neighbour_limit(3,8) .and. p1(2) < neighbour_limit(4,8)) &
			       .and. (p1(1) >neighbour_limit(1,8)  .and. p1(1) < neighbour_limit(2,8)) .or. &
			    				(p2(2) >  neighbour_limit(3,8) .and. p2(2) < neighbour_limit(4,8)) &
			       .and. (p2(1) >neighbour_limit(1,8)  .and. p2(1) < neighbour_limit(2,8)) .or. &
			       				(p3(2) >  neighbour_limit(3,8) .and. p3(2) < neighbour_limit(4,8)) &
			       .and. (p3(1) >neighbour_limit(1,8)  .and. p3(1) < neighbour_limit(2,8)) .or. &
			       				(p4(2) >  neighbour_limit(3,8) .and. p4(2) < neighbour_limit(4,8)) &
			       .and. (p4(1) >neighbour_limit(1,8)  .and. p4(1) < neighbour_limit(2,8)))
				
				in_down_left = ((p1(2) >  neighbour_limit(3,7) .and. p1(2) < neighbour_limit(4,7)) &
			       .and. (p1(1) > neighbour_limit(1,7) .and. p1(1) < neighbour_limit(2,7)) .or. &
			       				(p2(2) >  neighbour_limit(3,7) .and. p2(2) < neighbour_limit(4,7)) &
			       .and. (p2(1) > neighbour_limit(1,7) .and. p2(1) < neighbour_limit(2,7)) .or. &
			       				(p3(2) >  neighbour_limit(3,7) .and. p3(2) < neighbour_limit(4,7)) &
			       .and. (p3(1) > neighbour_limit(1,7) .and. p3(1) < neighbour_limit(2,7)) .or. &
			       				(p4(2) >  neighbour_limit(3,7) .and. p4(2) < neighbour_limit(4,7)) &
			       .and. (p4(1) > neighbour_limit(1,7) .and. p4(1) < neighbour_limit(2,7)))
				
				in_down_right = ((p1(2) >  neighbour_limit(3,6) .and. p1(2) < neighbour_limit(4,6)) &
			       .and. (p1(1) > neighbour_limit(1,6)  .and. p1(1) < neighbour_limit(2,6)) .or. &
			       				(p2(2) >  neighbour_limit(3,6) .and. p2(2) < neighbour_limit(4,6)) &
			       .and. (p2(1) > neighbour_limit(1,6)  .and. p2(1) < neighbour_limit(2,6)) .or. &
			       				(p3(2) >  neighbour_limit(3,6) .and. p3(2) < neighbour_limit(4,6)) &
			       .and. (p3(1) > neighbour_limit(1,6)  .and. p3(1) < neighbour_limit(2,6)) .or. &
			       				(p4(2) >  neighbour_limit(3,6) .and. p4(2) < neighbour_limit(4,6)) &
			       .and. (p4(1) > neighbour_limit(1,6)  .and. p4(1) < neighbour_limit(2,6)))
			       
!================================================================================================			       
			     IF (rank_up /= -1) then  

					IF(in_up ) then
						COUNTER_overlap(FJ,rank) = COUNTER_overlap(FJ,rank) + 1
						IND_overlap(COUNTER_overlap(FJ,rank),FJ) =  particle_k%ID		
					END IF	
				END IF

!================================================================================================			
				IF (rank_down /= -1) then
					IF (in_down  ) then
						COUNTER_overlap(BJ,rank) = COUNTER_overlap(BJ,rank) + 1
						IND_overlap(COUNTER_overlap(BJ,rank),BJ) =  particle_k%ID
					END IF					
				END IF
!================================================================================================
				IF (rank_left /= -1) then
					IF (in_left ) then
						COUNTER_overlap(BK,rank) = COUNTER_overlap(BK,rank) + 1
						IND_overlap(COUNTER_overlap(BK,rank),BK) =  particle_k%ID
					END IF			
				END IF
!================================================================================================
				IF (rank_right/= -1) then
					IF (in_right) then
					 	COUNTER_overlap(FK,rank) = COUNTER_overlap(FK,rank) + 1
						IND_overlap(COUNTER_overlap(FK,rank),FK) =  particle_k%ID
					END IF					
				END IF
!================================================================================================
				IF (rank_up_left /= -1) then
					IF (in_up_left) then
			        	COUNTER_overlap(FJBK,rank) = COUNTER_overlap(FJBK,rank) + 1
						IND_overlap(COUNTER_overlap(FJBK,rank),FJBK) =  particle_k%ID
					END IF
				END IF
!================================================================================================
				IF (rank_up_right /= -1) then
					IF (in_up_right) then
			        	COUNTER_overlap(FJFK,rank) = COUNTER_overlap(FJFK,rank) + 1
						IND_overlap(COUNTER_overlap(FJFK,rank),FJFK) =  particle_k%ID
			        END IF			    
				END IF
!================================================================================================
				IF (rank_down_left /= -1) then
					IF (in_down_left) then
						COUNTER_overlap(BJBK,rank) = COUNTER_overlap(BJBK,rank) + 1
						IND_overlap(COUNTER_overlap(BJBK,rank),BJBK) =  particle_k%ID
					END IF	
				END IF
!================================================================================================
				IF (rank_down_right /= -1) then
					IF (in_down_right) then 
						COUNTER_overlap(BJFK,rank) = COUNTER_overlap(BJFK,rank) + 1
						IND_overlap(COUNTER_overlap(BJFK,rank),BJFK) =  particle_k%ID
					END IF
				END IF

			if(particle_k%ID==3) then
				print*,'ssssssssss',in_up,particle_k%pointu1,particle_k%pointu2,particle_k%pointu3,particle_k%pointu4
			end if			

			END IF 
!================================================================================================
		END SUBROUTINE
		
		
		SUBROUTINE send_overlap_and_danger_particle
		integer :: i,neighloop
		integer, dimension(NB_NEIGHBOURS) :: OPP
		
		OPP = (/2,1,4,3,6,5,8,7/)		

		COUNTER_recv_danger = 0
		
		COUNTER_recv_overlap =0


!================================================================================================
			DO neighloop = 1,8
			IF (NEIGHBOUR(neighloop)<0) THEN
						cycle
					END IF
				DO i=1,COUNTER_overlap(neighloop,rank)
					
					IF (NEIGHBOUR(neighloop)<0) THEN
						cycle
					END IF
					call MPI_SEND(IND_overlap(i,neighloop),1,MPI_INTEGER,NEIGHBOUR(neighloop),4,MPI_COMM_WORLD,code)				
				END DO
			END DO
			
!================================================================================================
			DO neighloop = 1,8
			IF (NEIGHBOUR(neighloop)<0) THEN
						cycle
					END IF
			
				DO i= 1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
					IF (NEIGHBOUR(neighloop)<0) THEN
						CYCLE
					END IF
					call MPI_RECV(IND_recv_overlap(i,neighloop),1,MPI_INTEGER,NEIGHBOUR(neighloop),4, & 
					MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)	
					
				END DO
				
				IF (NEIGHBOUR(neighloop)<0) THEN
						CYCLE
				ELSE
				COUNTER_recv_overlap = COUNTER_recv_overlap + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
				END IF
				
			END DO	
!================================================================================================
			
			DO neighloop = 1,8
				IF (NEIGHBOUR(neighloop)<0) THEN
					cycle
				END IF
				DO i=1,COUNTER_danger(neighloop,rank)
					IF (NEIGHBOUR(neighloop)<0) THEN
						cycle
					END IF
					call MPI_SEND(IND_danger(i,neighloop),1,MPI_INTEGER,NEIGHBOUR(neighloop),5,MPI_COMM_WORLD,code)
					
				END DO
			END DO
			
!================================================================================================			
			DO neighloop = 1,8
				IF (NEIGHBOUR(neighloop)<0) THEN
					cycle
				END IF
				
				DO i= 1,COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
					IF (NEIGHBOUR(neighloop)<0) THEN
						CYCLE
					END IF
					call MPI_RECV(IND_recv_danger(i,neighloop),1,MPI_INTEGER,NEIGHBOUR(neighloop),5, & 
					MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
				END DO
				IF (NEIGHBOUR(neighloop)<0) THEN
						CYCLE
				ELSE
				COUNTER_recv_danger = COUNTER_recv_danger + COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
				END IF
			END DO


!================================================================================================
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

		! We will have "Npart_block" became: COUNTER_inside + sum of COUNTER_overlap in 
		! all direction + COUNTER_recv_overlap + COUNTER_recv_danger
		
		! Where COUNTER_danger is the number of particle inside of neighbour block but can 
		! effect on the particles inside the present block
		
		! Where COUNTER_recv_overlap = sum of COUNTER_overlap from neigbour 

		Npart_block = COUNTER_inside(rank) + COUNTER_recv_danger + COUNTER_recv_overlap  
!		+ SUM(COUNTER_overlap(:,rank))
		

if(rank==12) then
!DO i =1,8
!	
!	print*,'neighbour',i,IND_recv_overlap(:,i)
!	
!END DO*


end if

print*,rank,'Npart',Npart_block
		
			ALLOCATE(Xpart_block(2,Npart_block))
			ALLOCATE(Mblock(Npart_block))
			ALLOCATE(Dblock(4,Npart_block))
			ALLOCATE(velocity_field(2,Npart_block))
			ALLOCATE(Dout(4,Npart_block))
			

!================================================================================================

		
		DO i= 1,COUNTER_inside(rank)
			Xpart_block(1,i) = ALL_PARTICLES(IND_inside(i))%Xp
			Xpart_block(2,i) = ALL_PARTICLES(IND_inside(i))%Yp

		END DO
		Ncum1 = Ncum1 + COUNTER_inside(rank)
		
		DO neighloop=1,8
		if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
		DO j=1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
				Xpart_block(1,j+Ncum1) = ALL_PARTICLES(IND_recv_overlap(j,neighloop))%Xp
				Xpart_block(2,j+Ncum1) = ALL_PARTICLES(IND_recv_overlap(j,neighloop))%Yp	



			END DO
			Ncum1 = Ncum1 + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))

		END DO
		
		DO neighloop=1,8
		IF (NEIGHBOUR(neighloop)<0) then
				cycle
			END IF
		DO j=1,COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
				Xpart_block(1,j+Ncum1) = ALL_PARTICLES(IND_recv_danger(j,neighloop))%Xp
				Xpart_block(2,j+Ncum1) = ALL_PARTICLES(IND_recv_danger(j,neighloop))%Yp

			END DO
			Ncum1 = Ncum1 + COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
		END DO
					




!================================================================================================
		
		DO i= 1,COUNTER_inside(rank)
			Mblock(i) = ALL_PARTICLES(IND_inside(i))%mass

		END DO
		Ncum2 = Ncum2 + COUNTER_inside(rank)
		
		DO neighloop=1,8
		if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
		DO j=1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
				Mblock(j+Ncum2) = ALL_PARTICLES(IND_recv_overlap(j,neighloop))%mass
			END DO
			Ncum2 = Ncum2 + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))	
		END DO
		

		DO neighloop=1,8
		IF (NEIGHBOUR(neighloop)<0) then
				cycle
			END IF
		DO j=1,COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
				Mblock(j+Ncum2) = ALL_PARTICLES(IND_recv_danger(j,neighloop))%mass
			END DO
			Ncum2 = Ncum2 + COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
		END DO

!================================================================================================	


		DO i= 1,COUNTER_inside(rank)

			Dblock(1,i) = ALL_PARTICLES(IND_inside(i))%Dp1
			Dblock(2,i) = ALL_PARTICLES(IND_inside(i))%Dp2
			Dblock(3,i) = ALL_PARTICLES(IND_inside(i))%Dp3
			Dblock(4,i) = ALL_PARTICLES(IND_inside(i))%Dp4
		END DO	

		Ncum3 = Ncum3 + COUNTER_inside(rank)	
		
		DO neighloop=1,8
		if (NEIGHBOUR(neighloop)<0) then
				cycle
			end if
		DO j=1,COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))
				Dblock(1,j+Ncum3) = ALL_PARTICLES(IND_recv_overlap(j,neighloop))%Dp1
				Dblock(2,j+Ncum3) = ALL_PARTICLES(IND_recv_overlap(j,neighloop))%Dp2
				Dblock(3,j+Ncum3) = ALL_PARTICLES(IND_recv_overlap(j,neighloop))%Dp3
				Dblock(4,j+Ncum3) = ALL_PARTICLES(IND_recv_overlap(j,neighloop))%Dp4

			END DO
			Ncum3 = Ncum3 + COUNTER_overlap(OPP(neighloop),NEIGHBOUR(neighloop))	
		END DO


			
		DO neighloop=1,8
		IF (NEIGHBOUR(neighloop)<0) then
				cycle
			END IF
		DO j=1,COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
				Dblock(1,j+Ncum3) = ALL_PARTICLES(IND_recv_danger(j,neighloop))%Dp1
				Dblock(2,j+Ncum3) = ALL_PARTICLES(IND_recv_danger(j,neighloop))%Dp2
				Dblock(3,j+Ncum3) = ALL_PARTICLES(IND_recv_danger(j,neighloop))%Dp3
				Dblock(4,j+Ncum3) = ALL_PARTICLES(IND_recv_danger(j,neighloop))%Dp4

			END DO
			Ncum3 = Ncum3 + COUNTER_danger(OPP(neighloop),NEIGHBOUR(neighloop))
		END DO

!================================================================================================
		velocity_field = 0.
!================================================================================================
!IF(Npic==4.and.rank==4)then
!	print*,'npic:',Npic,Npart_block
!END IF
!IF(Npic==3.and.rank==4)then
!	print*,'npic:',Npic,Npart_block
!END IF	

!if(rank==8) then
!DO i=1,9
!print*,	Xpart_block(:,i)
!END DO
!end if

		call diffusion_field_ltp(Xpart_block,Mblock,Dblock,hx,hy,velocity_field,Npart_block)

if(COUNTER_inside(rank) + SUM(COUNTER_overlap(:,rank))>0) then

		call update_d_diffusion(Xpart_block, Dblock, Mblock, T_start, time_step, time_scheme, hx, hy, &
     indice_max_norm_Dm1, Norm_inf_Dm1, Dout, Npart_block)

end if


!================================================================================================
	
		DO i=1,COUNTER_inside(rank) 
			Xpart_block(1,i) = Xpart_block(1,i) - velocity_field(1,i)*time_step
			Xpart_block(2,i) = Xpart_block(2,i) - velocity_field(2,i)*time_step
		END DO
		

!================================================================================================

		DO i= 1,COUNTER_inside(rank)
			ALL_PARTICLES(IND_inside(i))%Xp = Xpart_block(1,i)
			ALL_PARTICLES(IND_inside(i))%Yp = Xpart_block(2,i) 
		END DO
		
!================================================================================================

		DO i= 1,COUNTER_inside(rank)
			ALL_PARTICLES(IND_inside(i))%Dp1 = Dout(1,i)
			ALL_PARTICLES(IND_inside(i))%Dp2 = Dout(2,i)
			ALL_PARTICLES(IND_inside(i))%Dp3 = Dout(3,i)
			ALL_PARTICLES(IND_inside(i))%Dp4 = Dout(4,i)
		END DO
		
		
		CALL send_all_block_information_to_rank_0_and_broadcast	
			
!================================================================================================
			DEALLOCATE(Dout)
			DEALLOCATE(Xpart_block)
			DEALLOCATE(Mblock)
			DEALLOCATE(Dblock)
			DEALLOCATE(velocity_field)	
		
		END SUBROUTINE block_loop_on_block_global_table	
		
!================================================================================================

		SUBROUTINE update_ALL_particles
			IMPLICIT NONE
			LOGICAL                           :: overlap_check, totally_inside_check, & 
												 danger_check, leave_check
			logical                           :: local_overlapcheck, local_insidecheck, & 
			 									 local_dangercheck,local_leave_check
			double precision, dimension(2)    :: pointu1,pointu2,pointu3,pointu4
			integer                           :: j,neighloop

			integer, dimension(NB_NEIGHBOURS) :: OPP	
		
			double precision                  :: axe
			integer                           :: ID
			
			COUNTER_recv_leave   = 0

			OPP = (/2,1,4,3,6,5,8,7/)

!=================================================================================================
			COUNTER_danger(:,:) = 0
			COUNTER_leave(:,:)  = 0


			
			i=1
			DO WHILE (i<=COUNTER_inside(rank))
			
				call overlap_criterion(ALL_PARTICLES(IND_inside(i)),overlap_check,totally_inside_check,danger_check, & 
				pointu1,pointu2,pointu3,pointu4,axe)

					ALL_PARTICLES(IND_inside(i))%pointu1 = pointu1
					ALL_PARTICLES(IND_inside(i))%pointu2 = pointu2
					ALL_PARTICLES(IND_inside(i))%pointu3 = pointu3
					ALL_PARTICLES(IND_inside(i))%pointu4 = pointu4
					ALL_PARTICLES(IND_inside(i))%axe = axe
				IF (totally_inside_check) then
					leave_check = .false.
					call particle_screening(ALL_PARTICLES(IND_inside(i)),overlap_check,danger_check,leave_check)
					i = i + 1
				else
					IF(overlap_check) then
						leave_check = .false.
						call particle_screening(ALL_PARTICLES(IND_inside(i)),overlap_check,danger_check,leave_check)
						IND_inside(i:COUNTER_inside(rank)-1) = IND_inside(i+1:COUNTER_inside(rank))
						COUNTER_inside(rank) = COUNTER_inside(rank) - 1
					else
						
						leave_check = .true.
						
						call particle_screening(ALL_PARTICLES(IND_inside(i)),overlap_check,danger_check,leave_check)
						
						IND_inside(i:COUNTER_inside(rank)-1) = IND_inside(i+1:COUNTER_inside(rank))
						COUNTER_inside(rank) = COUNTER_inside(rank) - 1
						
						
					END IF
				END IF

			END DO	
		

!================================================================================================

!			DO neighloop=1,8
!			IF (NEIGHBOUR(neighloop)<0) then
!				cycle
!			END IF
!			j=1
!			DO WHILE (j<=COUNTER_overlap(neighloop,rank))

!				call overlap_criterion(ALL_PARTICLES(IND_overlap(j,neighloop)),overlap_check,totally_inside_check,danger_check, & 
!				pointu1,pointu2,pointu3,pointu4,axe)
!				ALL_PARTICLES(IND_overlap(j,neighloop))%pointu1 = pointu1
!				ALL_PARTICLES(IND_overlap(j,neighloop))%pointu2 = pointu2
!				ALL_PARTICLES(IND_overlap(j,neighloop))%pointu3 = pointu3
!				ALL_PARTICLES(IND_overlap(j,neighloop))%pointu4 = pointu4
!				ALL_PARTICLES(IND_overlap(j,neighloop))%axe     = axe
!				
!				IF (totally_inside_check) then
!					leave_check                = .false.
!					call particle_screening(ALL_PARTICLES(IND_overlap(j,neighloop)),overlap_check,danger_check,leave_check)
!					COUNTER_inside(rank)             = COUNTER_inside(rank) + 1
!					IND_inside(COUNTER_inside(rank)) = IND_overlap(j,neighloop)
!					IND_overlap(j:COUNTER_overlap(neighloop,rank)-1,neighloop) &
!						 = IND_overlap(j+1:COUNTER_overlap(neighloop,rank),neighloop)
!					COUNTER_overlap(neighloop,rank) = COUNTER_overlap(neighloop,rank) -1
!				else
!					IF(overlap_check) then
!						j = j + 1
!					else
!						leave_check = .true.
!						call particle_screening(ALL_PARTICLES(IND_overlap(j,neighloop)),overlap_check,danger_check,leave_check)
!						IND_overlap(j:COUNTER_overlap(neighloop,rank)-1,neighloop) = & 
!						IND_overlap(j+1:COUNTER_overlap(neighloop,rank),neighloop)
!						COUNTER_overlap(neighloop,rank) = COUNTER_overlap(neighloop,rank) -1

!					END IF
!				END IF
!			END DO
!			END DO

			IF(rank /= 0) then
      			call MPI_SEND(COUNTER_inside(rank),NB_NEIGHBOURS,MPI_INTEGER,0,7,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1

       				call MPI_RECV(COUNTER_inside(ID),NB_NEIGHBOURS,MPI_INTEGER,ID,7,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			END do
    		END IF
    		call MPI_BCAST(COUNTER_inside,nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)

!================================================================================================    
	
    		IF(rank /= 0) then
      			call MPI_SEND(COUNTER_leave(1,rank),NB_NEIGHBOURS,MPI_INTEGER,0,7,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1

       				call MPI_RECV(COUNTER_leave(1,ID),NB_NEIGHBOURS,MPI_INTEGER,ID,7,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			END do
    		END IF
    		call MPI_BCAST(COUNTER_leave,NB_NEIGHBOURS*nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)



!================================================================================================    		
    		IF(rank /= 0) then
      			call MPI_SEND(COUNTER_danger(1,rank),NB_NEIGHBOURS,MPI_INTEGER,0,8,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1
       				call MPI_RECV(COUNTER_danger(1,ID),NB_NEIGHBOURS,MPI_INTEGER,ID,8,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			END do
    		END IF
    		call MPI_BCAST(COUNTER_danger,NB_NEIGHBOURS*nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)
    			

!================================================================================================		
		
			IF(rank /= 0) then
      			call MPI_SEND(COUNTER_overlap(1,rank),NB_NEIGHBOURS,MPI_INTEGER,0,9,MPI_COMM_WORLD,code)
    		else
      			do ID = 1,nb_proc - 1				
       				call MPI_RECV(COUNTER_overlap(1,ID),NB_NEIGHBOURS,MPI_INTEGER,ID,9,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
      			END do	
    		END IF
    		call MPI_BCAST(COUNTER_overlap,NB_NEIGHBOURS*nb_proc,MPI_INTEGER,0,MPI_COMM_WORLD,code)


!================================================================================================
			DO neighloop = 1,8
			IF(NEIGHBOUR(neighloop)<0) then
						cycle
					END IF
				
				DO i=1,COUNTER_leave(neighloop,rank)
					
					
					call MPI_SEND(IND_leave(i,neighloop),1,MPI_INTEGER,NEIGHBOUR(neighloop),10,MPI_COMM_WORLD,code)
					
				END DO
				
				DO i=1,COUNTER_leave(OPP(neighloop),NEIGHBOUR(neighloop))
					IF(NEIGHBOUR(neighloop)<0) then
						cycle
					END IF
					call MPI_RECV(IND_recv_leave(i,neighloop),1,MPI_INTEGER,NEIGHBOUR(neighloop),10 & 
					,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
				
						
				END DO
				
				IF (NEIGHBOUR(neighloop)<0) then
					cycle
				else
					
					COUNTER_recv_leave = COUNTER_recv_leave + COUNTER_leave(OPP(neighloop),NEIGHBOUR(neighloop))
				END IF	
			END DO

			
!================================================================================================
			
			DO neighloop=1,8
				IF (NEIGHBOUR(neighloop)<0) then
					cycle
				END IF
			
				DO j=1,COUNTER_leave(OPP(neighloop),NEIGHBOUR(neighloop))
				
					call overlap_criterion(ALL_PARTICLES(IND_recv_leave(j,neighloop)),local_overlapcheck, & 
					local_insidecheck,local_dangercheck,pointu1,pointu2,pointu3,pointu4,axe)
					
					IF (local_insidecheck) then
						COUNTER_inside(rank) = COUNTER_inside(rank) + 1
						IND_inside(COUNTER_inside(rank)) = IND_recv_leave(j,neighloop)	
						local_leave_check= .false.
						call particle_screening(ALL_PARTICLES(IND_recv_leave(j,neighloop)), & 
						local_overlapcheck,local_dangercheck,local_leave_check)
					else 
						local_leave_check= .false.
						call particle_screening(ALL_PARTICLES(IND_recv_leave(j,neighloop)), & 
						local_overlapcheck,local_dangercheck,local_leave_check)
				END IF
					
				END DO
		
		END DO
		
		END SUBROUTINE update_ALL_particles
!===============================================================================================
	
		SUBROUTINE sEND_all_block_information_to_rank_0_and_broadcast
		
		integer                          :: neighloop
		integer                          :: i,j,ID
		integer,dimension(1:nb_proc-1)   ::address
		
		
		address(:) = 0

		call MPI_PART_TYPE
		
	
!================================================================================================
		IF (rank /= 0) then
			DO i=1,COUNTER_inside(rank)
				address(rank) = IND_inside(i)
				call MPI_SEND(address(rank),1,MPI_INTEGER,0,11,MPI_COMM_WORLD,code)
				call MPI_SEND(ALL_PARTICLES(address(rank)),1,MPI_PARTICLETYPE,0,12,MPI_COMM_WORLD,code)
			END DO
		else
			DO ID =1,nb_proc-1
				DO i=1,COUNTER_inside(ID)
					call MPI_RECV(address(ID),1,MPI_INTEGER,ID,11,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
					call MPI_RECV(ALL_PARTICLES(address(ID)),1,MPI_PARTICLETYPE,ID,12,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
				END DO
			END DO
		END IF	
!================================================================================================
			
!			address(:) = 0

!		IF (rank /=0) then
!			DO neighloop=1,8
!				IF (NEIGHBOUR(neighloop)<0) then
!							cycle
!				END IF
!				DO j=1,COUNTER_overlap(neighloop,rank)
!					address(rank) = IND_overlap(j,neighloop)
!					call MPI_SEND(address(rank),1,MPI_INTEGER,0,13,MPI_COMM_WORLD,code)
!					call MPI_SEND(ALL_PARTICLES(address(rank)),1,MPI_PARTICLETYPE,0,14,MPI_COMM_WORLD,code)	
!				END DO
!			END DO
!		else 
!			DO ID =1,nb_proc-1
!				DO neighloop=1,8
!					DO j=1,COUNTER_overlap(neighloop,ID)				
!						call MPI_RECV(address(ID),1,MPI_INTEGER,ID,13,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
!						call MPI_RECV(ALL_PARTICLES(address(ID)),1,MPI_PARTICLETYPE,ID,14,MPI_COMM_WORLD,MPI_STATUS_IGNORE,code)
!					END DO
!				END DO
!			END DO
!		END IF


	
	call MPI_BCAST(ALL_PARTICLES,number_of_particles,MPI_PARTICLETYPE,0,MPI_COMM_WORLD,code)



!================================================================================================	
				

		END SUBROUTINE sEND_all_block_information_to_rank_0_and_broadcast
		
!================================================================================================
		SUBROUTINE update_all_particle_information
			
		IF(rank==0) then
			DO i=1,number_of_particles

				Xparticle_read(1,i) = ALL_PARTICLES(i)%Xp
				Xparticle_read(2,i) = ALL_PARTICLES(i)%Yp		

				D_read(1,i) = ALL_PARTICLES(i)%Dp1
				D_read(2,i) = ALL_PARTICLES(i)%Dp2
				D_read(3,i) = ALL_PARTICLES(i)%Dp3 
				D_read(4,i) = ALL_PARTICLES(i)%Dp4
			END DO
	
!================================================================================================
			OPEN(unit = 2, FILE = 'trunk/fortran_srcs/coords4fortran.bin',FORM="UNFORMATTED", &
 STATUS="UNKNOWN",POSITION="REWIND", ACTION="READWRITE", ACCESS='STREAM')
			WRITE(2) Xparticle_read
			CLOSE(2)
!================================================================================================			
			OPEN(unit= 3, FILE='trunk/fortran_srcs/deformmatrix4fortran.bin',FORM="UNFORMATTED",STATUS="UNKNOWN",POSITION="REWIND" & 
,ACTION="READWRITE",ACCESS='STREAM')			
			
			write(3) D_read
			close(3)
!================================================================================================			
		END IF
		
		END SUBROUTINE update_all_particle_information



		SUBROUTINE dealloc_all_particle_table
			DEALLOCATE(IND_leave)
			DEALLOCATE(IND_inside)
			DEALLOCATE(IND_overlap)
			DEALLOCATE(IND_danger)

			DEALLOCATE(ALL_PARTICLES)
			
			DEALLOCATE(COUNTER_danger)		
			DEALLOCATE(COUNTER_leave)
			DEALLOCATE(COUNTER_overlap)

			DEALLOCATE(IND_recv_danger)
			DEALLOCATE(IND_recv_overlap)
		
			
		END SUBROUTINE dealloc_all_particle_table

END MODULE mpi_2d_structures_modf90
