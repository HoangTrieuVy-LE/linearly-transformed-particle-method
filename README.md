# linearly-transformed-particle-method
Parallelization inside of LTP method

Global review

==================== PARTICLES INITIALISATION ================

Get all particle data from python and transfer into fortran code

======================== SPATIALISATION ====================

1. We assign every block(casier) into a process based on a MPI_CART_CREATE rule.
This step will be call "topology_initialisation"

2. We divide the domain into block, we will call the number of block is nb_proc, which is also the number of process. By using "set_block_grid(mesh(1,1),mesh(1,2),mesh(2,1),mesh(2,2))"
Each block has:
+ start_x, end_x
+ start_y, end_y

3. This step, we will find all possible neighbour blocks of running block. Each block will have 8 maximum neighbours in 2D. By using "neighbour_blocks" subroutine :
+ rank_left 
+ rank_right 
+ rank_up 
+ rank_down 
+ rank_up_left 
+ rank_up_right 
+ rank_down_left 
+ rank_down_right 

=====> We can test by count the number of neighbours of a running block by "neighbour_counter" subroutine.

 
 
 ============== PARTICLES DISTRIBUTION &  OVERLAP PARTICLES LISTS CREATION =============

1. First of all, we would like to renumerize the "ID" of neigbour blocks by "neighbouring" subroutine in order to easier work with MPI_SEND and MPI_RECV.

2. Initiate all necessary tables:
! IND_leave indicates the identity(in ALL_PARTICLES) of particles will leave the present block	
! IND_recv_overlap indicates the identity of particles receving from neighbours	
! IND_inside indicates the identity of particles inside the present block	
! IND_overlap indicates the identity of overlap particles in each direction	

! COUNTER_inside indicates the number of particles inside the present block
! COUNTER_leave indicates the number of particles leaving in each direction
! COUNTER_recv_overlap indicates the number of overlap particles receiving in each direction 

! ALL_PARTICLES stock all particle informations and should be broadcasted after being updated all particles displacement.

3. Convert all particle data information (X,M,D,..) into derived parttype.

4. Using "particle_distribution" to create:
FOR EVERY PROCESS (BLOCK)
+ IND_overlap
+ IND_inside
Them use "send_overlap_particle" to create the last important table:
+ IND_recv_overlap

========= LOOP ON PARTICLES INSIDE BLOCK AND THE OTHERS IN OVERLAP TABLES ============

Take a loop all over the particles inside block and all neighbour block's particles.





