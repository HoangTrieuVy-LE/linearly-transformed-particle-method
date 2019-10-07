!---------------------------------------------------------------------------  
!> @author 
!> Joris BOUTELOUP.
!
! DESCRIPTION: 
!> Decomposition of the domain and determination of the neighbor processors.
!
! REVISION HISTORY:
! 30/05/2013 
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!
!--------------------------------------------------------------------------- 
subroutine MPI_decomposition
	!=================
	!   DECLARATIONS
	!=================
	use mod_dem
	!
	implicit none
	!
	integer                         :: iargc,nb_args,i
	character(len=10), dimension(3) :: arg
	!
	!
	nbdims = 3
	allocate(dims(nbdims),periods(nbdims),coords(nbdims))
	dims(:) = 1
	!
	! lecture eventuelle du decoupage sur la ligne de commande (si rien indique, decoupage suivant x)
	nb_args = iargc()
	if ( nb_args /= 0 ) then
		do i=1,3
			call getarg(i,arg(i))
			read(arg(i),*) dims(i)
		end do
	else
		dims(1) = nb_proc
		dims(2) = 1
		dims(3) = 1
	end if
	!
	! Initialisation
	rang_gauche = MPI_PROC_NULL
	rang_droit = MPI_PROC_NULL
	rang_haut = MPI_PROC_NULL
	rang_bas = MPI_PROC_NULL
	rang_avant = MPI_PROC_NULL
	rang_arriere = MPI_PROC_NULL
	
	rang_h_gauche = MPI_PROC_NULL
	rang_h_droit = MPI_PROC_NULL
	rang_b_gauche = MPI_PROC_NULL
	rang_b_droit = MPI_PROC_NULL
	
	rang_av_h_gauche = MPI_PROC_NULL
	rang_av_h_droit = MPI_PROC_NULL
	rang_av_b_gauche = MPI_PROC_NULL
	rang_av_b_droit = MPI_PROC_NULL
	rang_av_gauche = MPI_PROC_NULL
	rang_av_droit = MPI_PROC_NULL
	rang_av_haut = MPI_PROC_NULL
	rang_av_bas = MPI_PROC_NULL
	
	rang_ar_h_gauche = MPI_PROC_NULL
	rang_ar_h_droit = MPI_PROC_NULL
	rang_ar_b_gauche = MPI_PROC_NULL
	rang_ar_b_droit = MPI_PROC_NULL
	rang_ar_gauche = MPI_PROC_NULL
	rang_ar_droit = MPI_PROC_NULL
	rang_ar_haut = MPI_PROC_NULL
	rang_ar_bas = MPI_PROC_NULL
	!
	!Creation of the cartesian topology
	!
	call MPI_DIMS_CREATE(nb_proc, nbdims, dims, code)
	!
	periods(:) = .false.
	!
	call MPI_CART_CREATE(MPI_COMM_WORLD, nbdims, dims, periods, &
	reorg, comm3d, code)
	!
	!Calculation of coordinates in the topology
	call MPI_CART_COORDS(comm3d, rang, nbdims, coords, code)
	!
	! Recherche de mes voisins Gauche et Droite
	if ( dims(1) > 1) then
		call MPI_CART_SHIFT(comm3d,0,1,rang_gauche,rang_droit,code)
		!
		!> Periodicity
		if (coords(1) == 0 .and. xper == 1) then
			rang_gauche = rang + dims(3)*dims(2)*(dims(1)-1)
		end if
		!
		if (coords(1) == dims(1)-1 .and. xper == 1) then
			rang_droit = rang - dims(3)*dims(2)*(dims(1)-1)
		end if
	end if
	!
	! Recherche de mes voisins Haut et Bas
	if (dims(2) > 1) then
		call MPI_CART_SHIFT(comm3d,1,1,rang_bas,rang_haut,code)
		!
		!> Periodicity
		if (coords(2) == 0 .and. yper == 1) then
			rang_bas = rang + (dims(2) - 1)*dims(3)
		end if
		!
		if (coords(2) == dims(2)-1 .and. yper == 1) then
			rang_haut = rang - (dims(2) - 1)*dims(3)
		end if
	end if
	!
	! Recherche de mes voisins Avant et Arriere
	if (dims(3) > 1) then
		call MPI_CART_SHIFT(comm3d,2,1,rang_avant,rang_arriere,code)
		!
		!> Periodicity
		if (coords(3) == 0 .and. zper == 1) then
			rang_avant = rang + dims(3)-1
		end if
		!
		if (coords(3) == dims(3)-1 .and. zper == 1) then
			rang_arriere = rang - dims(3)+1
		end if
	end if
	!
	if (dims(1) > 1 .and. dims(2) > 1) then
		! Recherche des voisins aux coins
		if (coords(1) .ne. 0 .and. coords(2) .ne. dims(2)-1) then
			rang_h_gauche = rang - (dims(2)-1)*dims(3)
		else
			if (xper == 1 .and. coords(1) == 0 .and. coords(2) .ne. dims(2)-1) then
				rang_h_gauche = rang + dims(3)*dims(2)*(dims(1)-1) + dims(3)
			end if
			if (yper == 1 .and. coords(1) .ne. 0 .and. coords(2) == dims(2)-1) then
				rang_h_gauche = rang - (2*dims(2)-1)*dims(3)
			end if
			if (yper*xper == 1 .and. coords(1) == 0 .and. coords(2) == dims(2)-1) then
				rang_h_gauche = rang + dims(3)*dims(2)*(dims(1)-1) - (dims(2)-1)*dims(3)
			end if
		end if
		!
		if (coords(1) .ne. dims(1)-1 .and. coords(2) .ne. dims(2)-1) then
			rang_h_droit = rang + dims(2)*dims(3) + dims(3)
		else
			if (xper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) .ne. dims(2)-1) then
				rang_h_droit = rang + dims(3) - dims(3)*dims(2)*(dims(1)-1)
			end if
			if (yper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) == dims(2)-1) then
				rang_h_droit = rang + dims(3)
			end if
			if (yper*xper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) == dims(2)-1) then
				rang_h_droit = rang - dims(3)*dims(2)*dims(1) + dims(3)
			end if
		end if
		!
		if (coords(1) .ne. 0 .and. coords(2) .ne. 0) then
			rang_b_gauche = rang - dims(2)*dims(3) - dims(3)
		else
			if (xper == 1 .and. coords(1) == 0 .and. coords(2) .ne. 0) then
				rang_b_gauche = rang - dims(3) + dims(3)*dims(2)*(dims(1)-1)
			end if
			if (yper == 1 .and. coords(1) .ne. 0 .and. coords(2) == 0) then
				rang_b_gauche = rang - dims(3)
			end if
			if (yper*xper == 1 .and. coords(1) == 0 .and. coords(2) == 0) then
				rang_b_gauche = rang + dims(3)*dims(2)*dims(1) - dims(3)
			end if
		end if
		!
		if (coords(1) .ne. dims(1)-1 .and. coords(2) .ne. 0) then
			rang_b_droit = rang + dims(2)*dims(3) - dims(3)
		else
			if (xper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) .ne. 0) then
				rang_b_droit = rang - dims(3) - dims(3)*dims(2)*(dims(1)-1)
			end if
			if (yper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) == 0) then
				rang_b_droit = rang + 2*dims(2)*dims(3) - dims(3)
			end if
			if (yper*xper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) == 0) then
				rang_b_droit = rang - dims(3)*dims(2)*(dims(1)-1) + (dims(2)-1)*dims(3)
			end if
		end if
	end if
	!
	if (dims(3) > 1 .and. dims(1) > 1 .and. dims(2) > 1) then
		!
		if (coords(1) .ne. 0 .and. coords(2) .ne. dims(2)-1 .and. coords(3) .ne. 0) then
			rang_av_h_gauche = rang - (dims(2)-1)*dims(3) - 1
		else
			if (xper == 1 .and. coords(1) == 0 .and. coords(2) .ne. dims(2)-1 .and. coords(3) .ne. 0) then
				rang_av_h_gauche = rang + dims(3)*dims(2)*(dims(1)-1) + dims(3) -1
			end if
			if (yper == 1 .and. coords(1) .ne. 0 .and. coords(2) == dims(2)-1 .and. coords(3) .ne. 0) then
				rang_av_h_gauche = rang -1 - (dims(2)-1)*dims(3) - dims(3)*dims(2)
			end if
			if (zper == 1 .and. coords(1) .ne. 0 .and. coords(2) .ne. dims(2)-1 .and. coords(3) == 0) then
				rang_av_h_gauche = rang - (dims(2)-1)*dims(3) + dims(3) - 1
			end if
			if (xper*yper == 1 .and. coords(1) == 0 .and. coords(2) == dims(2)-1 .and. coords(3) .ne. 0) then
				rang_av_h_gauche = rang + dims(3)*dims(2)*(dims(1)-1) - (dims(2)-1)*dims(3) - 1
			end if
			if (xper*zper == 1 .and. coords(1) == 0 .and. coords(2) .ne. dims(2)-1 .and. coords(3) == 0) then
				rang_av_h_gauche = rang + dims(3)-1 + dims(3)*dims(2)*(dims(1)-1) + dims(3)
			end if
			if (yper*zper == 1 .and. coords(1) .ne. 0 .and. coords(2) == dims(2)-1 .and. coords(3) == 0) then
				rang_av_h_gauche = rang + dims(3)-1 - (dims(2)-1)*dims(3) - dims(2)*dims(3)
			end if
			if (xper*yper*zper == 1 .and. coords(1) == 0 .and. coords(2) == dims(2)-1 .and. coords(3) == 0) then
				rang_av_h_gauche = nb_proc-1 - dims(2)*dims(3) + dims(3)
			end if
		end if
		!
		if (coords(1) .ne. dims(1)-1 .and. coords(2) .ne. dims(2)-1 .and. coords(3) .ne. 0) then
			rang_av_h_droit = rang + dims(2)*dims(3) + dims(3) - 1
		else
			if (xper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) .ne. dims(2)-1 .and. coords(3) .ne. 0) then
				rang_av_h_droit = rang + dims(3) - dims(3)*dims(2)*(dims(1)-1) - 1
			end if
			if (yper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) == dims(2)-1 .and. coords(3) .ne. 0) then
				rang_av_h_droit = rang + dims(3) - 1
			end if
			if (zper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) .ne. dims(2)-1 .and. coords(3) == 0) then
				rang_av_h_droit = rang + dims(2)*dims(3) + dims(3) + dims(3)-1
			end if
			if (xper*yper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) == dims(2)-1 .and. coords(3) .ne. 0) then
				rang_av_h_droit = rang - dims(3)*dims(2)*dims(1) + dims(3) - 1
			end if
			if (xper*zper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) .ne. dims(2)-1 .and. coords(3) == 0) then
				rang_av_h_droit = rang + dims(3) - dims(3)*dims(2)*(dims(1)-1) + dims(3)-1
			end if
			if (yper*zper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) == dims(2)-1 .and. coords(3) == 0) then
				rang_av_h_droit = rang + dims(3) + dims(3)-1
			end if
			if (xper*yper*zper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) == dims(2)-1 .and. coords(3) == 0) then
				rang_av_h_droit = dims(3)-1
			end if
		end if
		!
		if (coords(1) .ne. 0 .and. coords(2) .ne. 0 .and. coords(3) .ne. 0) then
			rang_av_b_gauche = rang - dims(2)*dims(3) - dims(3) - 1
		else
			if (xper == 1 .and. coords(1) == 0 .and. coords(2) .ne. 0 .and. coords(3) .ne. 0) then
				rang_av_b_gauche = rang - dims(3) + dims(3)*dims(2)*(dims(1)-1) - 1
			end if
			if (yper == 1 .and. coords(1) .ne. 0 .and. coords(2) == 0 .and. coords(3) .ne. 0) then
				rang_av_b_gauche = rang - dims(3) - 1
			end if
			if (zper == 1 .and. coords(1) .ne. 0 .and. coords(2) .ne. 0 .and. coords(3) == 0) then
				rang_av_b_gauche = rang - dims(2)*dims(3) - dims(3) + dims(3)-1
			end if
			if (xper*yper == 1 .and. coords(1) == 0 .and. coords(2) == 0 .and. coords(3) .ne. 0) then
				rang_av_b_gauche = rang + dims(3)*dims(2)*dims(1) - dims(3) - 1
			end if
			if (xper*zper == 1 .and. coords(1) == 0 .and. coords(2) .ne. 0 .and. coords(3) == 0) then
				rang_av_b_gauche = rang + dims(3)*dims(2)*(dims(1)-1) - 1
			end if
			if (yper*zper == 1 .and. coords(1) .ne. 0 .and. coords(2) == 0 .and. coords(3) == 0) then
				rang_av_b_gauche = rang - 1
			end if
			if (xper*yper*zper == 1 .and. coords(1) == 0 .and. coords(2) == 0 .and. coords(3) == 0) then
				rang_av_b_gauche = nb_proc-1
			end if
		end if
		!
		if (coords(1) .ne. dims(1)-1 .and. coords(2) .ne. 0 .and. coords(3) .ne. 0) then
			rang_av_b_droit = rang + dims(2)*dims(3) - dims(3) - 1
		else
			if (xper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) .ne. 0 .and. coords(3) .ne. 0) then
				rang_av_b_droit = rang - dims(3) - dims(3)*dims(2)*(dims(1)-1) - 1
			end if
			if (yper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) == 0 .and. coords(3) .ne. 0) then
				rang_av_b_droit = rang + 2*dims(2)*dims(3) - dims(3) - 1
			end if
			if (zper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) .ne. 0 .and. coords(3) == 0) then
				rang_av_b_droit = rang + dims(2)*dims(3) - dims(3) + dims(3)-1
			end if
			if (xper*yper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) == 0 .and. coords(3) .ne. 0) then
				rang_av_b_droit = rang - dims(3)*dims(2)*(dims(1)-1) + (dims(2)-1)*dims(3) - 1
			end if
			if (xper*zper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) .ne. 0 .and. coords(3) == 0) then
				rang_av_b_droit = rang - dims(3)*dims(2)*(dims(1)-1) - 1
			end if
			if (yper*zper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) == 0 .and. coords(3) == 0) then
				rang_av_b_droit = rang + 2*dims(2)*dims(3) - 1
			end if
			if (xper*yper*zper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) == 0 .and. coords(3) == 0) then
				rang_av_b_droit = dims(2)*dims(3)-1
			end if
		end if
		!
		if (coords(1) .ne. 0 .and. coords(2) .ne. dims(2)-1 .and. coords(3) .ne. dims(3)-1) then
			rang_ar_h_gauche = rang - (dims(2)-1)*dims(3) + 1
		else
			if (xper == 1 .and. coords(1) == 0 .and. coords(2) .ne. dims(2)-1 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_h_gauche = rang + dims(3)*dims(2)*(dims(1)-1) + dims(3) + 1
			end if
			if (yper == 1 .and. coords(1) .ne. 0 .and. coords(2) == dims(2)-1 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_h_gauche = rang + 1 - (dims(2)-1)*dims(3) - dims(3)*dims(2)
			end if
			if (zper == 1 .and. coords(1) .ne. 0 .and. coords(2) .ne. dims(2)-1 .and. coords(3) == dims(3)-1) then
				rang_ar_h_gauche = rang - (dims(2)-1)*dims(3) - dims(3)+1
			end if
			if (xper*yper == 1 .and. coords(1) == 0 .and. coords(2) == dims(2)-1 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_h_gauche = rang + dims(3)*dims(2)*(dims(1)-1) - (dims(2)-1)*dims(3) + 1
			end if
			if (xper*zper == 1 .and. coords(1) == 0 .and. coords(2) .ne. dims(2)-1 .and. coords(3) == dims(3)-1) then
				rang_ar_h_gauche = rang - dims(3)+1 + dims(3)*dims(2)*(dims(1)-1) + dims(3)
			end if
			if (yper*zper == 1 .and. coords(1) .ne. 0 .and. coords(2) == dims(2)-1 .and. coords(3) == dims(3)-1) then
				rang_ar_h_gauche = rang - dims(3)+1 - (dims(2)-1)*dims(3) - dims(2)*dims(3)
			end if
			if (xper*yper*zper == 1 .and. coords(1) == 0 .and. coords(2) == dims(2)-1 .and. coords(3) == dims(3)-1) then
				rang_ar_h_gauche = (dims(1)-1)*dims(2)*dims(3)
			end if
		end if
		!
		if (coords(1) .ne. dims(1)-1 .and. coords(2) .ne. dims(2)-1 .and. coords(3) .ne. dims(3)-1) then
			rang_ar_h_droit = rang + dims(2)*dims(3) + dims(3) + 1
		else
			if (xper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) .ne. dims(2)-1 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_h_droit = rang + dims(3) - dims(3)*dims(2)*(dims(1)-1) + 1
			end if
			if (yper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) == dims(2)-1 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_h_droit = rang + dims(3) + 1
			end if
			if (zper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) .ne. dims(2)-1 .and. coords(3) == dims(3)-1) then
				rang_ar_h_droit = rang + dims(2)*dims(3) + dims(3) - dims(3)+1
			end if
			if (xper*yper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) == dims(2)-1 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_h_droit = rang - dims(3)*dims(2)*dims(1) + dims(3) + 1
			end if
			if (xper*zper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) .ne. dims(2)-1 .and. coords(3) == dims(3)-1) then
				rang_ar_h_droit = rang + dims(3) - dims(3)*dims(2)*(dims(1)-1) - dims(3)+1
			end if
			if (yper*zper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) == dims(2)-1 .and. coords(3) == dims(3)-1) then
				rang_ar_h_droit = rang + dims(3) - dims(3)+1
			end if
			if (xper*yper*zper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) == dims(2)-1 .and. coords(3) == dims(3)-1) then
				rang_ar_h_droit = 0
			end if
		end if
		!
		if (coords(1) .ne. 0 .and. coords(2) .ne. 0 .and. coords(3) .ne. dims(3)-1) then
			rang_ar_b_gauche = rang - dims(2)*dims(3) - dims(3) + 1
		else
			if (xper == 1 .and. coords(1) == 0 .and. coords(2) .ne. 0 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_b_gauche = rang - dims(3) + dims(3)*dims(2)*(dims(1)-1) + 1
			end if
			if (yper == 1 .and. coords(1) .ne. 0 .and. coords(2) == 0 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_b_gauche = rang - dims(3) + 1
			end if
			if (zper == 1 .and. coords(1) .ne. 0 .and. coords(2) .ne. 0 .and. coords(3) == dims(3)-1) then
				rang_ar_b_gauche = rang - dims(2)*dims(3) - dims(3) - dims(3)+1
			end if
			if (xper*yper == 1 .and. coords(1) == 0 .and. coords(2) == 0 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_b_gauche = rang + dims(3)*dims(2)*dims(1) - dims(3) + 1
			end if
			if (xper*zper == 1 .and. coords(1) == 0 .and. coords(2) .ne. 0 .and. coords(3) == dims(3)-1) then
				rang_ar_b_gauche = rang - dims(3)+1 - dims(3) + dims(3)*dims(2)*(dims(1)-1)
			end if
			if (yper*zper == 1 .and. coords(1) .ne. 0 .and. coords(2) == 0 .and. coords(3) == dims(3)-1) then
				rang_ar_b_gauche = rang - dims(3)+1 - dims(3)
			end if
			if (xper*yper*zper == 1 .and. coords(1) == 0 .and. coords(2) == 0 .and. coords(3) == dims(3)-1) then
				rang_ar_b_gauche = nb_proc - dims(3)
			end if
		end if
		!
		if (coords(1) .ne. dims(1)-1 .and. coords(2) .ne. 0 .and. coords(3) .ne. dims(3)-1) then
			rang_ar_b_droit = rang + dims(2)*dims(3) - dims(3) + 1
		else
			if (xper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) .ne. 0 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_b_droit = rang - dims(3) - dims(3)*dims(2)*(dims(1)-1) + 1
			end if
			if (yper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) == 0 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_b_droit = rang + 2*dims(2)*dims(3) - dims(3) + 1
			end if
			if (zper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) .ne. 0 .and. coords(3) == dims(3)-1) then
				rang_ar_b_droit = rang + dims(2)*dims(3) - dims(3) - dims(3)+1
			end if
			if (xper*yper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) == 0 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_b_droit = rang - dims(3)*dims(2)*(dims(1)-1) + (dims(2)-1)*dims(3) + 1
			end if
			if (xper*zper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) .ne. 0 .and. coords(3) == dims(3)-1) then
				rang_ar_b_droit = rang - dims(3)+1 - dims(3) - dims(3)*dims(2)*(dims(1)-1)
			end if
			if (yper*zper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(2) == 0 .and. coords(3) == dims(3)-1) then
				rang_ar_b_droit = rang - dims(3)+1 + 2*dims(2)*dims(3)-dims(3)
			end if
			if (xper*yper*zper == 1 .and. coords(1) == dims(1)-1 .and. coords(2) == 0 .and. coords(3) == dims(3)-1) then
				rang_ar_b_droit = dims(3)*(dims(2)-1)
			end if
		end if
	end if
	!
	if (dims(3) > 1 .and. dims(1) > 1) then
		!
		if (coords(1) .ne. 0 .and. coords(3) .ne. 0) then
			rang_av_gauche = rang - dims(3)*dims(2) - 1
		else
			if (xper == 1 .and. coords(1) == 0 .and. coords(3) .ne. 0) then
				rang_av_gauche = rang + dims(3)*dims(2)*(dims(1)-1) - 1
			end if
			if (zper == 1 .and. coords(1) .ne. 0 .and. coords(3) == 0) then
				rang_av_gauche = rang - dims(3)*dims(2) + dims(3)-1
			end if
			if (xper*zper == 1 .and. coords(1) == 0 .and. coords(3) == 0) then
				rang_av_gauche = rang + dims(3)*dims(2)*(dims(1)-1) + dims(3)-1
			end if
		end if
		!
		if (coords(1) .ne. dims(1)-1 .and. coords(3) .ne. 0) then
			rang_av_droit = rang + dims(3)*dims(2) - 1
		else
			if (xper == 1 .and. coords(1) == dims(1)-1 .and. coords(3) .ne. 0) then
				rang_av_droit = rang - dims(3)*dims(2)*(dims(1)-1) - 1
			end if
			if (zper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(3) == 0) then
				rang_av_droit = rang + dims(3)*dims(2) + dims(3)-1
			end if
			if (xper*zper == 1 .and. coords(1) == dims(1)-1 .and. coords(3) == 0) then
				rang_av_droit = rang - dims(3)*dims(2)*(dims(1)-1) + dims(3)-1
			end if
		end if
		!
		if (coords(1) .ne. 0 .and. coords(3) .ne. dims(3)-1) then
			rang_ar_gauche = rang - dims(3)*dims(2) + 1
		else
			if (xper == 1 .and. coords(1) == 0 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_gauche = rang + dims(3)*dims(2)*(dims(1)-1) + 1
			end if
			if (zper == 1 .and. coords(1) .ne. 0 .and. coords(3) == dims(3)-1) then
				rang_ar_gauche = rang - dims(3)*dims(2) - dims(3)+1
			end if
			if (xper*zper == 1 .and. coords(1) == 0 .and. coords(3) == dims(3)-1) then
				rang_ar_gauche = rang + dims(3)*dims(2)*(dims(1)-1) - dims(3)+1
			end if
		end if
		!
		if (coords(1) .ne. dims(1)-1 .and. coords(3) .ne. dims(3)-1) then
			rang_ar_droit = rang + dims(3)*dims(2) + 1
		else
			if (xper == 1 .and. coords(1) == dims(1)-1 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_droit = rang - dims(3)*dims(2)*(dims(1)-1) + 1
			end if
			if (zper == 1 .and. coords(1) .ne. dims(1)-1 .and. coords(3) == dims(3)-1) then
				rang_ar_droit = rang + dims(3)*dims(2) - dims(3)+1
			end if
			if (xper*zper == 1 .and. coords(1) == dims(1)-1 .and. coords(3) == dims(3)-1) then
				rang_ar_droit = rang - dims(3)*dims(2)*(dims(1)-1) - dims(3)+1
			end if
		end if
	end if
	!
	if (dims(3) > 1 .and. dims(2) > 1) then
		!
		if (coords(2) .ne. 0 .and. coords(3) .ne. 0) then
			rang_av_bas = rang - dims(3) - 1
		else
			if (yper == 1 .and. coords(2) == 0 .and. coords(3) .ne. 0) then
				rang_av_bas = rang + (dims(2)-1)*dims(3) - 1
			end if
			if (zper == 1 .and. coords(2) .ne. 0 .and. coords(3) == 0) then
				rang_av_bas = rang - dims(3) + dims(3)-1
			end if
			if (yper*zper == 1 .and. coords(2) == 0 .and. coords(3) == 0) then
				rang_av_bas = rang + (dims(2)-1)*dims(3) + dims(3)-1
			end if
		end if
		!
		if (coords(2) .ne. dims(2)-1 .and. coords(3) .ne. 0) then
			rang_av_haut = rang + dims(3) - 1
		else
			if (yper == 1 .and. coords(2) == dims(2)-1 .and. coords(3) .ne. 0) then
				rang_av_haut = rang - (dims(2)-1)*dims(3) - 1
			end if
			if (zper == 1 .and. coords(2) .ne. dims(2)-1 .and. coords(3) == 0) then
				rang_av_haut = rang + dims(3) + dims(3)-1
			end if
			if (yper*zper == 1 .and. coords(2) == dims(2)-1 .and. coords(3) == 0) then
				rang_av_haut = rang - (dims(2)-1)*dims(3) + dims(3)-1
			end if
		end if
		!
		if (coords(2) .ne. 0 .and. coords(3) .ne. dims(3)-1) then
			rang_ar_bas = rang - dims(3) + 1
		else
			if (yper == 1 .and. coords(2) == 0 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_bas = rang + (dims(2)-1)*dims(3) + 1
			end if
			if (zper == 1 .and. coords(2) .ne. 0 .and. coords(3) == dims(3)-1) then
				rang_ar_bas = rang - dims(3) - dims(3)+1
			end if
			if (yper*zper == 1 .and. coords(2) == 0 .and. coords(3) == dims(3)-1) then
				rang_ar_bas = rang + (dims(2)-1)*dims(3) - dims(3)+1
			end if
		end if
		!
		if (coords(2) .ne. dims(2)-1 .and. coords(3) .ne. dims(3)-1) then
			rang_ar_haut = rang + dims(3) + 1
		else
			if (yper == 1 .and. coords(2) == dims(2)-1 .and. coords(3) .ne. dims(3)-1) then
				rang_ar_haut = rang - (dims(2)-1)*dims(3) + 1
			end if
			if (zper == 1 .and. coords(2) .ne. dims(2)-1 .and. coords(3) == dims(3)-1) then
				rang_ar_haut = rang + dims(3) - dims(3)+1
			end if
			if (yper*zper == 1 .and. coords(2) == dims(2)-1 .and. coords(3) == dims(3)-1) then
				rang_ar_haut = rang - (dims(2)-1)*dims(3) - dims(3)+1
			end if
		end if
	end if
	!
end subroutine MPI_decomposition
