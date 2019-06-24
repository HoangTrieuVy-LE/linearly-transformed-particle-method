!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> MPI 2D structures 
!--------------------------------------------------------------------------- 


module mpi_2d_structures_modf90

	! Import modules, TODO decide which module would be necessary
	use calculsfor_rec_modf90
	use calculsfor_var_modf90
	use calculsfor_ini_modf90
	use PARTICLE_2D_modf90

	implicit none

	! DECALARATIONS
	! TODO


	! TODO derived particle dataypes
	!type(PARTTYPE), dimension(:,:), allocatable                      :: PARTICLES

	contains
		subroutine initiation_overlap
		! index of particles which stay on the same processor
		! index of overlap particles 
		! numbers of overlaps in each direction for each processor
		
		! At the beginning, we initialize a VIDE overlap particle table and a PLEIN inside_particlestable
		! array of particle inside a domain
		! array of overlap particle of a domain 
		end subroutine initiation_overlap
	
	
		subroutine particles_counting
		
		end subroutine particles_counting
	
		
		
		
		subroutine overlap_criterion()
		!particle_k,eigenvector_1,eigenvector_2
		! Check a particle if it satifies overlap_criterion
		! Retrieve the eigenvalues and eigenvectors to determine the direction movement of particles, 	then decide it leaves the block or not
		! Input: k-th_particle in "PARTYPE" type 
		! Output: logical
			implicit none
			
			! Initial idea
			integer, parameter                                      :: n=2
			double precision                                        :: a(n,n), x(n,n)
			double precision, parameter                             :: abserr=1.0e-09
			integer i, j

			! type(PARTTYPE), intent(in)   							:: particle_k
			!double precision, dimension(n), intent(out)             :: eigenvector_1
			!double precision, dimension(n), intent(out)             :: eigenvector_2
			double precision                                        :: eigenvalue_1
			double precision                                        :: eigenvalue_2
			logical                                                 :: overlap_check	

			! matrix A
  			! data (a(1,i), i=1,2) /   D11,  D12 /
  			! data (a(2,i), i=1,2) /   D21,  D22 /


			! print a header and the original matrix
			! call Jacobi(a,x,abserr,n)

			! eigenvector_1 = x(1,:) 
			! eigenvector_2 = x(2,:)
			! eigenvalue_1  = a(1,1)
			! eigenvalue_2  = a(2,2)		
			! pointu = eigenvector_1 + eigenvector_2
			! if ((pointu(1)-start_x)*(pointu(2)-start_y)>0 .and. (pointu(1)-end_x)*(pointu(2)-end_y)>0) then 
			! overlap_check = false
			! else 
			! overlap_check =true 
			! end if 
			
			
		end subroutine overlap_criterion


		subroutine particle_numerotation()
			!Xkpart, M, DD, Nk_part,table_block_particle,table_block_overlap
			! Give Identity particle inside block ?
			! Input: Xk(1,:),Xk(2,:), M, DD, Nk_part
			!double precision, dimension(2,Nk_part), intent(in)      :: Xkpart
			!double precision, dimension(Nk_part), intent(in)        :: M
			!double precision, dimension(4, Nk_part), intent(in)     :: DD
			! Output: table_block_particles which has Nk particles of "PARTYPE" type.
		end subroutine


		subroutine particle_screening()
			!table_block_particle
			!type(PARTYPE), dimension(:,:), allocatable, intent(in)  :: table_block_particle
			!type(PARTYPE), dimension(:,:), allocatable, intent(out) :: table_particle_inside
			!type(PARTYPE), dimension(:,:), allocatable, intent(out) :: table_particle_overlap
			! Loop on table_block_particles
			! Check overlap_criterion
			! Yes: append to table_particle_inside
			! No : append to table_particle_overlap
		end subroutine
		
		
		


end module mpi_2d_structures_modf90
