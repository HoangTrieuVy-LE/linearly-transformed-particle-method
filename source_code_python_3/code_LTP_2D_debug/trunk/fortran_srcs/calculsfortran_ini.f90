module calculsfor_ini_modf90

use calculsfor_rec_modf90

implicit none

  double precision, dimension(:,:), allocatable        :: xpart
  double precision, dimension(:), allocatable          :: mpart
  integer                                    :: num_rhoini
  double precision                                     :: r, lx1, lx2, ly1, ly2
  double precision, dimension(2,2)                     :: domaine
  integer                                    :: deg_spline
  integer                                    :: l ! precision pour le remapping : l  pts de grilles entre chaque particules, l>=2
contains

!-------------------------------------------------------
    subroutine  set_num_rhoini(n, rr, llx1, llx2, lly1, lly2)
      integer, intent(in)  :: n
      double precision, intent(in)   :: rr, llx1, llx2, lly1, lly2
      num_rhoini = n
      r = rr
      lx1 = llx1
      lx2 = llx2
      ly1 = lly1
      ly2 = lly2
    end subroutine
!-------------------------------------------------------

    subroutine rhoini(z, x, y)
      double precision, intent(in)   :: x, y
      double precision, intent(out)  :: z
  
      double precision   :: n
   
      if (num_rhoini==1) then 
          n = sqrt(x*x+y*y)
          if (n<=r) then
             z = 1.
          else
             z = 0.
          end if
      else if (num_rhoini==2) then 

      else if (num_rhoini==3) then 
      end if

    end subroutine

!-------------------------------------------------------
!    initialisation sans mettre de particules sur les bords, calcul des position et poids deterministes
    subroutine initialize_spline(doma, Nx_part, Ny_part, deg_spli)
      double precision, dimension(2,2), intent(in)   :: doma
      integer, intent(in)  :: Nx_part, Ny_part, deg_spli  

      double precision  :: hx, hy
      double precision  :: x, y, z, z00, zp0, z0p, zpp, zm0, zmm, z0m, zpm, zmp, t00, t01, t11
      double precision  :: xmin, xmax, ymin, ymax 
      integer :: i, j, k

      xmin = doma(1,1) 
      xmax = doma(1,2)
      ymin = doma(2,1) 
      ymax = doma(2,2)
      hx = (xmax - xmin)/float(Nx_part)  !  distance between particles
      hy = (ymax - ymin)/float(Ny_part)

      k = 1
      do j=1, Ny_part
          if (Ny_part==1) then 
               y = (ymin+ymax)/2.
          else 
               y = ymax - (j + 0.5) * hy
          end if
          do i=1,Nx_part
             if (Nx_part==1) then 
                x = (xmin+xmax)/2.
             else 
                x = xmin + (i + 0.5) * hx
             end if
             xpart(1,k) = x
             xpart(2,k) = y
             if (deg_spli==1) then
                 call rhoini(z, x, y)
                 mpart(k) =  hx * hy * z 
             else 
                call rhoini(z00, x, y)
                call rhoini(zm0, x-hx, y)
                call rhoini(zp0, x+hx, y)
                call rhoini(z0m, x, y-hy)
                call rhoini(z0p, x, y+hy)
                call rhoini(zpp, x+hx, y+hy)
                call rhoini(zmm, x-hx, y-hy)
                call rhoini(zpm, x+hx, y-hy)
                call rhoini(zmp, x-hx, y+hy)
                t00 = z00 * (8./6.)**2
                t01 = -8./36. *(zm0+zp0+z0m+z0p)
                t11 =1./36.* (zmm+zpp+zmp+zpm)
                mpart(k) = hx * hy * max(t00+t01+t11,0.) 
              end if
              k = k+1
            end do
        end do

    end subroutine

!----------------------------------------------------------------
! remapping en utilisant l initialisation par spline
    subroutine set_param_remap(dom, deg_s, ll)
       double precision, dimension(2,2), intent(in)         :: dom
       integer, intent(in)                        :: deg_s, ll
       Domaine = dom
       deg_spline=  deg_s
       l = ll

    end subroutine

    subroutine remap(Xold, Mold, Dold, Nx_part, Ny_part, Nold)
      double precision, dimension(2,Nold), intent(in)   :: Xold
      double precision, dimension(1,Nold), intent(in)   :: Mold
      double precision, dimension(4, Nold), intent(in)  :: Dold
      integer, intent(in)  :: Nx_part, Ny_part,  Nold


      double precision  :: hx, hy, deltax, deltay
      double precision  :: x, y, z, z00, zp0, z0p, zpp, zm0, zmm, z0m, zpm, zmp, t00, t01, t11
      double precision  :: xmin, xmax, ymin, ymax 
      integer :: i, j, k, indice_Xg_x, indice_Yg_y
      double precision, dimension(:), allocatable :: Xg, Yg
      integer :: NgrilleX, NgrilleY

      integer :: N_part
 
      if ((Ny_part==1).or.(Ny_part==1)) then 
           print*, 'prendre plus de particules pour le remapping !!'
           stop
      end if

      xmin = Domaine(1,1) 
      xmax = Domaine(1,2)
      ymin = Domaine(2,1) 
      ymax = Domaine(2,2)
      hx = (xmax - xmin)/float(Nx_part-1)  !  distance between particles
      hy = (ymax - ymin)/float(Ny_part-1)
      N_part = Nx_part * Ny_part

!----- definir la grille sur laquelle on reconstruit -----
! ---  cette grille contient les positions des futurs particules ----
      deltax = hx/(l+1)
      deltay = hy/(l+1)
      NgrilleX = Nx_part +l*(Nx_part-1)+2
      NgrilleY = Ny_part +l*(Ny_part-1)+2
      allocate(Xg(NgrilleX)) 
      allocate(Yg(NgrilleY)) 
      do k=1, NgrilleX
         Xg(k) = xmin-deltax + (k-1)*deltax
      end do
      do k=1, NgrilleY
         Yg(k) = ymin-deltay + (k-1)*deltay
      end do
! ------ reconstruit la densite sur la grille -----------
    call alloc_ro(NgrilleX, NgrilleY)
    call set_degre_spline(deg_spline)
    call set_hx_hy(hx, hy)
    call set_grille(Xg, Yg, NgrilleX, NgrilleY)
    call rec_ltp(Xold, Mold, Dold, Nold)  ! -> met des valeurs dans ro

! ------ calcule les positions et poids des particules

      k = 1
      do j=1, Ny_part
             y =  ymin + hy*(j-1)
             indice_Yg_y = (j-1)*(l+1)+2  ! indice dans le tableau yg
          do i=1,Nx_part
             x = xmin + (i -1)*hx
             indice_Xg_x = (i-1)*(l+1)+2
             xpart(1,k) = x
             xpart(2,k) = y
             if (deg_spline ==1) then
                 z = ro(indice_Xg_x, indice_Yg_y)
                 mpart(k) =  hx * hy * z 
             else 
                z00 = ro(indice_Xg_x, indice_Yg_y)
                zm0 = ro(indice_Xg_x-1, indice_Yg_y)
                zp0 = ro(indice_Xg_x+1, indice_Yg_y)
                z0m = ro(indice_Xg_x, indice_Yg_y-1)
                z0p = ro(indice_Xg_x, indice_Yg_y+1)
                zpp = ro(indice_Xg_x+1, indice_Yg_y+1)
                zmm = ro(indice_Xg_x-1, indice_Yg_y-1)
                zpm = ro(indice_Xg_x+1, indice_Yg_y-1)
                zmp = ro(indice_Xg_x-1, indice_Yg_y+1)
                t00 = z00 * (8./6.)**2
                t01 = -8./36. *(zm0+zp0+z0m+z0p)
                t11 =1./36.* (zmm+zpp+zmp+zpm)
                mpart(k) = hx * hy * max(t00+t01+t11,0.) 
              end if
              k = k+1
            end do
        end do

    call dealloc_ro()
    end subroutine
! -------------------------------------------------------------------
! ================= subroutine d allocation/desallocation
!--------------------------------------------------------------------
    subroutine alloc_vect(N_part) 
      integer, intent(in)                    :: N_part
 
      allocate(xpart(2,N_part))    
      xpart = 0.
      allocate(mpart(N_part))
      mpart=0.
      print*, 'coucou'
    end subroutine

    subroutine dealloc_vect() 

      deallocate(xpart)    
      deallocate(mpart)

    end subroutine


    

    
end module calculsfor_ini_modf90

