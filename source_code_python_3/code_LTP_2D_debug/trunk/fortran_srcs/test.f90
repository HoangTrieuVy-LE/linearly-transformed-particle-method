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
