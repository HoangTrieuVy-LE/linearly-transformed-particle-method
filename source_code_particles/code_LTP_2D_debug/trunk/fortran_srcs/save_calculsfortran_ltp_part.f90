module calculsfor_rec_modf90

implicit none

  real*8, dimension(:,:), allocatable        :: ro  
  integer                                    :: degre_spline 
  real*8                                     :: radius_phi
  real*8                                     :: hx, hy
  integer                                    :: nxg, nyg
  real*8, dimension (:),  allocatable        :: Xg, Yg
contains

!-------------------------------------------------------
   subroutine set_degre_spline(p)
      integer, intent(in)                    :: p

      degre_spline = p
      if (p==1) then
          radius_phi =1.
      else if (p==3) then
          radius_phi =2.
      end if 

   end subroutine

!-------------------------------------------------------
   subroutine set_hx_hy(val_hx, val_hy)
      real*8,  intent(in)  :: val_hx, val_hy
      hx = val_hx 
      hy = val_hy

   end subroutine

!-------------------------------------------------------
   subroutine set_grille(vect_xg, vect_yg, val_nxg, val_nyg)
      integer  ,  intent(in)  :: val_nxg, val_nyg
      real*8, dimension (val_nxg), intent(in)         :: vect_xg
      real*8, dimension (val_nyg), intent(in)         :: vect_yg

      nxg = val_nxg 
      nyg = val_nyg

      Xg = vect_xg
      Yg = vect_yg
   end subroutine

!-------------------------------------------------------
    subroutine phi(y,x)
       real*8,  intent(in)  :: x
       real*8,  intent(out) :: y
       real*8               :: xx

    xx = abs(x)
    if (degre_spline==1) then
         if (xx>=1) then
             y = 0.
         else
             y = 1 - xx
         end if
    else  
         if (xx>=2) then
             y = 0.
         else if (xx>=1) then
             y =  (2 - xx) * (2 - xx) * (2 - xx) / 6
         else 
             y = (4 - 6 * xx * xx + 3 * xx * xx * xx) / 6
         end if
    end if

    end subroutine



!-------------------------------------------------------
    subroutine rec_ltp(Xpart, M, DD, N_part)
      integer, intent(in)                         :: N_part 
      real*8, dimension (2,N_part), intent(in)    :: Xpart
      real*8, dimension (N_part), intent(in)      :: M
      real*8, dimension (4, N_part), intent(in)   :: DD


      integer(8)   :: i, j, k
      real*8       :: deltax, deltay, a, b, c, d, det_Dk
      real*8       :: rayonk, yj, xi
      integer(8)   :: qkx, qky, ik, jk
      real*8       :: phix, phjy

      deltax = abs(Xg(2) - Xg(1))
      deltay = abs(Yg(2) - Yg(1))
      ro = 0.

      do k=1,N_part
        a = DD(1, k)
        b = DD(2, k)   
        c = DD(3, k)
        d = DD(4, k)
        det_Dk = abs(a*d-b*c)
        if (det_Dk ==0) then
            print*, 'probleme de Dk non inversible'
        end if
        rayonk = max(hx,hy)/det_Dk*max(abs(d)+abs(b),abs(c)+abs(a))
        qkx = int(radius_phi * rayonk / deltax) + 1
        qky = int(radius_phi * rayonk / deltay) + 1
        ik = int((Xpart(1,k) - Xg(1)) / deltax)
        jk = int((Xpart(2,k) - Yg(1) ) / deltay)
        do j=max(1,jk- qky), min(jk+qky+1, nyg)
            yj = Yg(j)  
            do i= max(ik-qkx,1), min(ik+qkx+1,nxg)
                xi = Xg(i)

                call phi(phix, (a*(Xpart(1,k)-xi)+b*(Xpart(2,k) - yj)) / hx)
                call phi(phjy, (c*(Xpart(1,k)-xi)+d*(Xpart(2,k) - yj)) / hy)
                ro(i,j) = ro(i,j) + M(k)*det_Dk/(hx*hy) * phix * phjy                
            end do
        end do
      end do
    end subroutine


!-------------------------------------------------------
    subroutine rec_sansdepot(Xpart, M, DD, N_part)
      integer, intent(in)                         :: N_part 
      real*8, dimension (2,N_part), intent(in)    :: Xpart
      real*8, dimension (N_part), intent(in)      :: M
      real*8, dimension (4, N_part), intent(in)   :: DD


      integer(8)   :: i, j, k
      real*8       :: a, b, c, d, det_Dk
      real*8       :: phix, phjy, xi, yj
      real*8, dimension(2) :: Yk      

      ro = 0.
      print*, 'methode de reconstruction sans depot'
      do j=1,nyg
          yj = Yg(j)  
          do i=1,nxg
              xi = Xg(i)
              do k=1, N_part
                  a = DD(1, k)
                  b = DD(2, k)   
                  c = DD(3, k)
                  d = DD(4, k)
                  det_Dk = abs(a*d-b*c)
                  if (det_Dk ==0) then
                       print*, 'probleme de Dk non inversible'
                  end if
                  Yk(1) = a*(Xpart(1,k)-xi)+b*(Xpart(2,k) - yj)
                  Yk(2) = c*(Xpart(1,k)-xi)+d*(Xpart(2,k) - yj)
                  call phi(phix, Yk(1) / hx)
                  call phi(phjy, Yk(2) / hy)
                  ro(i,j) = ro(i,j) + M(k)*det_Dk/(hx*hy) * phix * phjy                
              end do
          end do
      end do

    end subroutine

!-------------------------------------------------------
    subroutine rec_sp(Xpart, M, epsi, N_part)
    ! ne marche pas, a  debugger
      integer, intent(in)                         :: N_part 
      real*8, dimension (2,N_part), intent(in)    :: Xpart
      real*8, dimension (N_part), intent(in)      :: M
      real*8, intent(in)                          :: epsi
      real*8, dimension(:,:), allocatable           :: depotx, depoty   

      integer(8)   :: i, j, k
      real*8       :: deltax, deltay, a, b, c, d, det_Dk
      real*8       :: rayonk, yj, xi
      integer(8)   :: qx, qy, ik, jk
      real*8       :: phix, phjy

      deltax = abs(Xg(2) - Xg(1))
      deltay = abs(Yg(2) - Yg(1))
      ro = 0.

      allocate(depotx(nxg, 1))
      allocate(depoty(1, nyg))
      !print*, 'epsilon=', epsi
      qx = int(radius_phi * epsi / deltax) + 1
      qy = int(radius_phi * epsi / deltay) + 1

      do k=1,N_part
        ik = int((Xpart(1,k) - Xg(1)) / deltax)
        jk = int((Xpart(2,k) - Yg(1) ) / deltay)
        do j=max(1,jk- qy), min(jk+qy+1, nyg)
            yj = Yg(j)
            call phi(phjy, (Xpart(2,k) - yj)/ epsi) 
            depoty(1,j) = phjy
        end do
        do i= max(ik-qx,1), min(ik+qx+1,nxg)
            xi = Xg(i)
            call phi(phix, (Xpart(1,k)-xi) / epsi)
            depotx(i,1) = phix
        end do
        !print*, 'max des depoty(j)/epsi', max(depoty)/epsi
        !print*, 'max des depotx(i)/epsi', max(depotx)/epsi
        ro = ro + M(k)/(epsi*epsi) * matmul(depotx,depoty)  
             
      end do

      deallocate(depotx)
      deallocate(depoty)

    end subroutine
! -------------------------------------------------------------------
! ================= subroutine d allocation/desallocation
!--------------------------------------------------------------------
    subroutine alloc_diff(N) 
      integer, intent(in)                    :: N      

      allocate(ro(dimx,dimy))    
      ro = 0.
      allocate(Xg(dimx))
      Xg=0.
      allocate(Yg(dimy))
      Yg=0.

    end subroutine

!-------------------------------------------------------------------
    subroutine dealloc_ro() 

      deallocate(ro)
      deallocate(Xg)
      deallocate(Yg)

    end subroutine

end module calculsfor_ltp_modf90



