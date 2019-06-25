!---------------------------------------------------------------------------  
!> @author 
!> HoangTrieuVy-LE.
!
! DESCRIPTION: 
!> MPI 2D structures 
!--------------------------------------------------------------------------- 


MODULE mpi_2d_structures_modf90

	IMPLICIT NONE
	
	! Import modules, TODO decide which module would be necessary
	USE mod_particle_2D_modf90
	
	
	! DECALARATIONS
	! TODO


	! TODO derived particle dataypes
	type(PARTTYPE), dimension(:,:), allocatable                      :: PARTICLES

	CONTAINS
		SUBROUTINE initiation_overlap
		! index of particles which stay on the same processor
		! index of overlap particles 
		! numbers of overlaps in each direction for each processor
		
		! At the beginning, we initialize a VIDE overlap particle table and a PLEIN inside_particlestable
		! array of particle inside a domain
		! array of overlap particle of a domain 
		END SUBROUTINE initiation_overlap
	
	
		SUBROUTINE particles_counting
		
		END SUBROUTINE particles_counting
		
		SUBROUTINE block_loop	
		! TODO

		! Input: X_casier_k (k is rank of processor),M_casier_k,D_old_casier_k,hx,hy, N_part_casier_k	
		! Output: U_casier_k, Df2py_casier_k	
	
		END SUBROUTINE block_loop	
		
		
		
		SUBROUTINE overlap_criterion()
		!particle_k,eigenvector_1,eigenvector_2
		! Check a particle if it satifies overlap_criterion
		! Retrieve the eigenvalues and eigenvectors to determine the direction movement of particles, 	then decide it leaves the block or not
		! Input: k-th_particle in "PARTYPE" type 
		! Output: logical
			IMPLICIT NONE
			
			! Initial idea
			INTEGER, PARAMETER                                      :: n=2
			DOUBLE PRECISION                                        :: a(n,n), x(n,n)
			DOUBLE PRECISION, PARAMETER                             :: abserr=1.0e-09
			INTEGER i, j

			! type(PARTTYPE), intent(in)   							:: particle_k
			!double precision, dimension(n), intent(out)             :: eigenvector_1
			!double precision, dimension(n), intent(out)             :: eigenvector_2
			DOUBLE PRECISION                                        :: eigenvalue_1
			DOUBLE PRECISION                                        :: eigenvalue_2
			LOGICAL                                                 :: overlap_check	

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
			
			
		END SUBROUTINE overlap_criterion


		SUBROUTINE particle_numerotation()
			!Xkpart, M, DD, Nk_part,table_block_particle,table_block_overlap
			! Give Identity particle inside block ?
			! Input: Xk(1,:),Xk(2,:), M, DD, Nk_part
			!double precision, dimension(2,Nk_part), intent(in)      :: Xkpart
			!double precision, dimension(Nk_part), intent(in)        :: M
			!double precision, dimension(4, Nk_part), intent(in)     :: DD
			! Output: table_block_particles which has Nk particles of "PARTYPE" type.
		END SUBROUTINE


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
