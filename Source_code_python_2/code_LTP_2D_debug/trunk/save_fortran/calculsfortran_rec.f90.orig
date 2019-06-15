module calculsfor_rec_modf90

implicit none

  double precision, dimension(:,:), allocatable        :: ro 
  integer                                    :: degre_spline 
  integer                                    :: phi_radial
  double precision                                     :: cst_int_phi_exp
  double precision                                     :: radius_phi
  double precision                                     :: hx, hy
  integer                                    :: nxg, nyg
  double precision, dimension (:),  allocatable        :: Xg, Yg
  
contains


!-------------------------------------------------------
! function phi is radial if (phi(x,y) = psi(x) * psi(y)
	subroutine set_phi_radial(bool_radial)
        integer, intent(in)	:: bool_radial
        
        phi_radial = bool_radial		
        
   end subroutine

!-------------------------------------------------------
	subroutine set_degre_spline(p)
		integer, intent(in)                    :: p      		
		degre_spline = p
		if (phi_radial == 1) then ! function phi is radial (phi(x,y) = psi(x) * psi(y)
			if (p==1) then ! order 1 spline shape function phi is radial
				radius_phi =1.
			else if (p==3) then ! order 3 spline shape function phi is radial
				radius_phi =2.
			else 
				print*, "ERROR p value must be either 1 or 3"
			end if			
		else if (phi_radial == 0) then	! exponential with compact support phi(x,y) = exp(-1/(1-x**2-y**2)) / int_{B_{0,1}} {exp(-1/(1-x**2-y**2))}
			radius_phi = 1.		
			cst_int_phi_exp = 0.46651239317 ! value of int_{B_{0,1}} {exp(-1/(1-x**2-y**2))}
		end if 
   end subroutine



!-------------------------------------------------------
   subroutine set_hx_hy(val_hx, val_hy)
      double precision,  intent(in)  :: val_hx, val_hy
      hx = val_hx 
      hy = val_hy

   end subroutine

!-------------------------------------------------------
   subroutine set_grille(vect_xg, vect_yg, val_nxg, val_nyg)
      integer  ,  intent(in)  :: val_nxg, val_nyg
      double precision, dimension (val_nxg), intent(in)         :: vect_xg
      double precision, dimension (val_nyg), intent(in)         :: vect_yg

      nxg = val_nxg 
      nyg = val_nyg

      Xg = vect_xg
      Yg = vect_yg
   end subroutine

!-------------------------------------------------------
	subroutine phi(y,x)
		double precision,  intent(in)  :: x
		double precision,  intent(out) :: y
		double precision               :: xx
		double precision               :: pi        
		pi = 4.*atan(1.)
        xx = abs(x)
		if (phi_radial == 1) then ! phi is a spline
			if (degre_spline==1) then
	            if (xx>=1) then
	                y = 0.
	            else
	                y = 1 - xx
	            end if
	        else if (degre_spline==3) then
	            if (xx>=2) then
	                y = 0.
	            else if (xx>=1) then
	                y =  (2 - xx) * (2 - xx) * (2 - xx) / 6
	            else 
	                y = (4 - 6 * xx * xx + 3 * xx * xx * xx) / 6
	            end if
	        else 
				print *,"ERROR degre_spline value must be either 1 or 3. Current degre_spline = " , degre_spline
	        end if
		else if (phi_radial == 0) then ! (phi is exponential with compact support )
			print *," TO BE DONE"
		else
			print *,"ERROR : phi_radial value must be either 0 or 1 (boolean). Current phi_radial = " , phi_radial 
		end if
    end subroutine
    !1D phi derivative
    ! ADD CASE degre_spline=1
    ! CHANGE AND DO NOT DIVIDE BY xx 
    subroutine derivative_phi(deriv,x)
        double precision, intent(in)  :: x
        double precision, intent(out) :: deriv
        double precision              :: xx
        double precision              :: pi        
        xx = abs(x)        
		pi = 4.*atan(1.)
        
        if (phi_radial == 1) then
	        if (degre_spline == 1) then
	            print *, " ERROR : case degre_spline = " , degre_spline , " not available yet ! "            
	        else if (degre_spline == 3) then 
	            if (xx >= 2) then
	                deriv = 0.
	            else if (xx >= 1) then
	                deriv =  -(x*(xx-2)*(xx-2))/(2*xx)
	            else
	                deriv = 0.5 * x * (3. *xx - 4.)
	            end if
	        else 
				print *,"ERROR degre_spline value must be either 1 or 3. Current degre_spline = " , degre_spline			
	        end if
		else if (phi_radial == 0) then ! (phi is exponential with compact support )
			print *," TO BE DONE"
		else
			print *,"ERROR : phi_radial value must be either 0 or 1 (boolean). Current phi_radial = " , phi_radial 
		end if
    end subroutine
    
    subroutine second_derivative_phi(scd_deriv,x)
        double precision, intent(in)  :: x
        double precision, intent(out) :: scd_deriv
        double precision              :: xx      
        double precision              :: pi        
		pi = 4.*atan(1.)  
        xx = abs(x)
        if (phi_radial == 1) then
	        if (degre_spline == 1) then
	            print *, " ERROR : case degre_spline = " , degre_spline , " not available yet ! "            
	        else if (degre_spline == 3) then 
	            if (xx >= 2) then            
	                scd_deriv = 0.
	            else if (xx >= 1) then            
	                scd_deriv =  2. - xx
	            !scd_deriv = ((3.) * (x**2) / (xx)) - 2.
	            else if (x > 0) then            
	                scd_deriv = 3.*x - 2.
	            else if (x < 0) then            
	                scd_deriv = -3.*x - 2.
	            else
	                scd_deriv=-2
	            end if	        
			else 
				print *,"ERROR degre_spline value must be either 1 or 3. Current degre_spline = " , degre_spline			
			end if
		else if (phi_radial == 0) then ! (phi is exponential with compact support )
			print *," TO BE DONE"
		else
			print *,"ERROR : phi_radial value must be either 0 or 1 (boolean). Current phi_radial = " , phi_radial 
		end if
    end subroutine	        
            


    !phi with 2D variables
    ! identical subroutine if degre_spline = 0 
    ! 
    subroutine phi_2d(phi_out,x,y)
        double precision,  intent(in)  :: x ,y
        double precision,  intent(out) :: phi_out
        double precision               :: phi_x , phi_y
        if (phi_radial == 1) then 
            call phi(phi_x,x)
            call phi(phi_y,y)
            phi_out = phi_x * phi_y
        else if (phi_radial == 0) then 
			!if (x**2 + y**2 >=1) then 
			if (x**2 + y**2 >=radius_phi) then
				phi_out = 0.
			else 				
				phi_out = (1./cst_int_phi_exp) * exp(-1./(radius_phi-(x**2+y**2))) 
			end if
		else
			print *,"ERROR : phi_radial value must be either 0 or 1 (boolean). Current phi_radial = " , phi_radial 
        end if
    end subroutine

    !gradient of phi_2d with 2D variables (2-size array) 
    subroutine grad_phi_2d(grad,x,y)
        double precision,  intent(in)                 :: x, y
        double precision, dimension (2), intent(out)  :: grad
        double precision                              :: phi_x, phi_y, grad_phi_x, grad_phi_y
        if (phi_radial == 1) then 
            call phi(phi_x,x)
            call phi(phi_y,y)
            call derivative_phi(grad_phi_x,x)
            call derivative_phi(grad_phi_y,y)
            grad(1)=grad_phi_x * phi_y
            grad(2)=grad_phi_y * phi_x
        else			
			if (x**2 + y**2 >=radius_phi) then
				grad = 0.
			else 
				grad(1)= (1./cst_int_phi_exp) * (-2*x/(radius_phi-(x**2+y**2))**2) * exp(-1/(radius_phi-x**2-y**2))
				grad(1)= (1./cst_int_phi_exp) * (-2*y/(radius_phi-(x**2+y**2))**2) * exp(-1/(radius_phi-x**2-y**2))
			end if
        end if                        
    end subroutine
    
    !hessien of phi_2d with 2D variables (2x2-size array) 
    subroutine hess_phi_2d(hess,x,y)
        double precision,  intent(in)                     :: x ,y
        double precision, dimension (2,2), intent(out)    :: hess
        double precision                                  :: phi_x, phi_y, deriv_phi_x , deriv_phi_y
        double precision                                  :: second_deriv_phi_x , second_deriv_phi_y
        if (phi_radial == 1) then 
            call phi(phi_x,x)
            call phi(phi_y,y)
            call derivative_phi(deriv_phi_x,x)
            call derivative_phi(deriv_phi_y,y)
            call second_derivative_phi(second_deriv_phi_x,x)
            call second_derivative_phi(second_deriv_phi_y,y)      
    !        print *, "x, y , phi_x, phi_y = " , x, y , phi_x, phi_y
    !        print *, "deriv_phi_x , deriv_phi_y = " , deriv_phi_x , deriv_phi_y
    !        print *, "second_deriv_phi_x, second_deriv_phi_y = " , second_deriv_phi_x, second_deriv_phi_y
            hess(1,1) = second_deriv_phi_x * phi_y
            hess(1,2) = deriv_phi_x * deriv_phi_y
            hess(2,1) = deriv_phi_y * deriv_phi_x
            hess(2,2) = phi_x * second_deriv_phi_y
        else if (phi_radial == 0) then 
			!if (x**2 + y**2 >=1) then
			if (x**2 + y**2 >=radius_phi) then
				hess = 0.
			else 			
	            hess(1,1) = (1./cst_int_phi_exp) * (2 * exp(1/(-radius_phi+x**2+y**2)) &
                * (3 *x**4+2 *x**2*y**2-(-1+y**2)**2))/(-radius_phi+x**2+y**2)**4
	            hess(1,2) = (1./cst_int_phi_exp) * (4*x*y* exp(1/(-radius_phi+x**2+y**2)) &
                * (2*x**2 + 2*y**2 - 1))/(-radius_phi+x**2+y**2)**4
	            hess(2,1) = (1./cst_int_phi_exp) * (4*x*y* exp(1/(-radius_phi+x**2+y**2)) &
                * (2*x**2 + 2*y**2 - 1))/(-radius_phi+x**2+y**2)**4
	            hess(2,2) = (1./cst_int_phi_exp) * (2 * exp(1/(-radius_phi+x**2+y**2)) &
                * (3 *y**4+2 *x**2*y**2-(-1+x**2)**2))/(-radius_phi+x**2+y**2)**4 ! MUST BE CHECKED (!!)
	       end if
        else 
			print *,"ERROR : phi_radial value must be either 0 or 1 (boolean). Current phi_radial = " , phi_radial 			
        end if
            
    end subroutine
        
    !phi_2d_h with 2D variables
    subroutine phi_2d_h(phi_h,x,y, hx_remap ,hy_remap)
        double precision,  intent(in)  :: x ,y , hx_remap , hy_remap
        double precision,  intent(out) :: phi_h
        double precision               :: phi_x_h , phi_y_h
        if (phi_radial == 1) then 
            call phi(phi_x_h,x/hx_remap)
            call phi(phi_y_h,y/hy_remap)
            phi_h = (1./hx_remap)*phi_x_h * (1./hy_remap)*phi_y_h
        else if (phi_radial == 0) then 
            call phi_2d(phi_h,x/hx_remap,y/hy_remap)
            phi_h = phi_h / (hx_remap * hy_remap)
        else 
			print *,"ERROR : phi_radial value must be either 0 or 1 (boolean). Current phi_radial = " , phi_radial
        end if
    end subroutine

    !gradient of phi_2d_h with 2D variables (2-size array) 
    subroutine grad_phi_2d_h(grad_h, x, y, hx_remap, hy_remap)
        double precision,  intent(in)                 :: x ,y , hx_remap , hy_remap
        double precision, dimension (2), intent(out)  :: grad_h    
        if (phi_radial == 1) then 
            call grad_phi_2d(grad_h, x/hx_remap , y/hy_remap)
            grad_h(1) = (1/(hx_remap*hy_remap)) * (1/hx_remap) * grad_h(1)
            grad_h(2) = (1/(hx_remap*hy_remap)) * (1/hy_remap) * grad_h(2)
        else if (phi_radial == 0) then 
            call grad_phi_2d(grad_h, x/hx_remap , y/hy_remap)
            grad_h(1) = grad_h(1) /(hy_remap * hx_remap**3)
            grad_h(2) = grad_h(2) /(hx_remap * hy_remap**3)
        else 
			print *,"ERROR : phi_radial value must be either 0 or 1 (boolean). Current phi_radial = " , phi_radial
        end if
    end subroutine
    
    subroutine hess_phi_2d_h(hess_h, x, y, hx_remap, hy_remap)
        double precision,  intent(in)                     :: x ,y , hx_remap , hy_remap
        double precision, dimension (2,2), intent(out)    :: hess_h        
        if (phi_radial == 1) then 
            call hess_phi_2d(hess_h,x/hx_remap,y/hy_remap)
            hess_h(1,1) = (1/(hx_remap*hy_remap))* (1/(hx_remap*hx_remap)) * hess_h(1,1)
            hess_h(1,2) = (1/(hx_remap*hy_remap)) * (1/(hx_remap*hy_remap)) * hess_h(1,2)
            hess_h(2,1) = (1/(hx_remap*hy_remap)) * (1/(hx_remap*hy_remap)) * hess_h(2,1)
            hess_h(2,2) = (1/(hx_remap*hy_remap))* (1/(hy_remap*hy_remap)) * hess_h(2,2)
        else if (phi_radial == 0) then 
            call hess_phi_2d(hess_h,x/hx_remap,y/hy_remap)
            hess_h(1,1) = 1./(hy_remap * hx_remap**3) * hess_h(1,1)
            hess_h(1,2) = 1./(hx_remap**3 * hy_remap**3)  * hess_h(1,2)
            hess_h(2,1) = 1./(hx_remap**3 * hy_remap**3) * hess_h(2,1)
            hess_h(2,2) = 1./(hx_remap * hy_remap**3) * hess_h(2,2)
        else 
			print *,"ERROR : phi_radial value must be either 0 or 1 (boolean). Current phi_radial = " , phi_radial
        end if
    end subroutine 
    
    
    !__________________________________________________________________!
    !
    ! reconstruct ltp density at a given coordinate
    !__________________________________________________________________!
    subroutine rec_ltp_at_coordinate(rho_hn, Xpart, M, DD, hx_remap, hy_remap,X_coord, N_part)
        integer, intent(in)                         :: N_part
        double precision, dimension (2,N_part), intent(in)    :: Xpart
        double precision, dimension (N_part), intent(in)      :: M
        double precision, dimension (4, N_part), intent(in)   :: DD
        double precision, intent(in)                          :: hx_remap , hy_remap
        double precision, dimension (2), intent(in)           :: X_coord
        double precision, intent(out)                         :: rho_hn
	    
        double precision                  :: phi_hn, det_Dkn
        double precision, dimension (2)   :: vect_var
        double precision, dimension (2,2) :: Dkn
        integer(8)              :: k, nb_proc, chunk
        
        integer, external       :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
        
        rho_hn=0.
        phi_hn=0.
        
!~         !$OMP PARALLEL 
!~         nb_proc = OMP_GET_NUM_THREADS()
!~         chunk = N_part / nb_proc +1
!~         !$OMP END PARALLEL 
        
        
!~ 		!$OMP PARALLEL DO shared(chunk) private(k,Dkn, det_Dkn, vect_var, phi_hn) schedule(dynamic,chunk)
        do k=1,N_part 
            Dkn(1,1) = DD(1,k)
            Dkn(1,2) = DD(2,k)
            Dkn(2,1) = DD(3,k)
            Dkn(2,2) = DD(4,k)
            det_Dkn = Dkn(1,1) * Dkn(2,2) - Dkn(2,1) * Dkn(1,2)
            !vect_var = MATMUL(TRANSPOSE(Dkn),X_coord - Xpart(:,k))
            vect_var = MATMUL(Dkn,X_coord - Xpart(:,k))
            call phi_2d_h( phi_hn , vect_var(1) , vect_var(2) , hx_remap, hy_remap)			
            rho_hn = rho_hn + M(k) * det_Dkn * phi_hn
        end do
!~        !$OMP END PARALLEL DO		
 
    end subroutine
    
    !__________________________________________________________________!
    !
    ! reconstruct ltp density of gradient at a given coordinate
    !__________________________________________________________________!
    subroutine rec_ltp_grad_at_coordinate(grad_rho_hn, Xpart, M, DD, hx_remap, hy_remap,X_coord, N_part, print_log_in)
        integer, intent(in)                         :: N_part
        integer, optional                           :: print_log_in
        double precision, dimension (2,N_part), intent(in)    :: Xpart
        double precision, dimension (N_part), intent(in)      :: M
        double precision, dimension (4, N_part), intent(in)   :: DD
        double precision, intent(in)                          :: hx_remap , hy_remap
        double precision, dimension (2), intent(in)           :: X_coord
        double precision, dimension (2), intent(out)          :: grad_rho_hn
	    
        double precision                  :: det_Dkn
        double precision, dimension (2)   :: vect_var , grad_phi_hn
        double precision, dimension (2,2) :: Dkn
        integer(8)              :: k, print_log, nb_proc, chunk
        
		integer, external       :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

        
        grad_rho_hn=0.
        grad_phi_hn=0.

!~ 		!$OMP PARALLEL 
!~         nb_proc = OMP_GET_NUM_THREADS()
!~         chunk = N_part / nb_proc +1
!~         !$OMP END PARALLEL 
        


		print_log = 0
		if (present(print_log_in)) then
			print_log = print_log_in
		end if		
		if (print_log == 1) then			
!~ 			!$OMP PARALLEL DO shared(chunk) private(k, Dkn, det_Dkn, vect_var, grad_phi_hn) schedule(dynamic,chunk)
			do k=1,N_part
	            Dkn(1,1) = DD(1,k)
	            Dkn(1,2) = DD(2,k)
	            Dkn(2,1) = DD(3,k)
	            Dkn(2,2) = DD(4,k)
	            det_Dkn = Dkn(1,1) * Dkn(2,2) - Dkn(2,1) * Dkn(1,2)	            
	            vect_var = MATMUL(Dkn,X_coord - Xpart(:,k))	            
	            call grad_phi_2d_h(grad_phi_hn , vect_var(1) , vect_var(2) , hx_remap, hy_remap)
	            grad_rho_hn = grad_rho_hn + det_Dkn * M(k) * MATMUL(TRANSPOSE(Dkn), grad_phi_hn)!				
!~ 				print *, " k , grad_rho_hn " , k , grad_rho_hn	            
				if ((abs(grad_phi_hn(1)) > 0) .OR. (abs(grad_phi_hn(2)) > 0)) then
					print *, "k , grad_phi_hn : " , k , grad_phi_hn(1) ,grad_phi_hn(2)
				end if
	        end do
!~ 			!$OMP END PARALLEL DO		
		else
!~ 			!$OMP PARALLEL DO shared(chunk) private(k, Dkn, det_Dkn, vect_var, grad_phi_hn) schedule(dynamic,chunk)
	        do k=1,N_part 
	            Dkn(1,1) = DD(1,k)
	            Dkn(1,2) = DD(2,k)
	            Dkn(2,1) = DD(3,k)
	            Dkn(2,2) = DD(4,k)
	            det_Dkn = Dkn(1,1) * Dkn(2,2) - Dkn(2,1) * Dkn(1,2)	            
	            vect_var = MATMUL(Dkn,X_coord - Xpart(:,k))	        
	            call grad_phi_2d_h(grad_phi_hn , vect_var(1) , vect_var(2) , hx_remap, hy_remap)
	            grad_rho_hn = grad_rho_hn + det_Dkn * M(k) * MATMUL(TRANSPOSE(Dkn), grad_phi_hn)
	            
	            !call rec_ltp_grad_at_coordinate(grad_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,j), N_part)
	            
	        end do        
!~ 	        !$OMP END PARALLEL DO		
        end if
    end subroutine
    
    
    !-------------------------------------------------------    
    subroutine rec_ltp_hess_at_coordinate(hess_rho_hn, Xpart, M, DD, hx_remap, hy_remap,X_coord, N_part)
        integer, intent(in)                         :: N_part
        double precision, dimension (2,N_part), intent(in)    :: Xpart
        double precision, dimension (N_part), intent(in)      :: M
        double precision, dimension (4, N_part), intent(in)   :: DD
        double precision, intent(in)                          :: hx_remap , hy_remap
        double precision, dimension (2), intent(in)           :: X_coord
        double precision, dimension (2,2), intent(out)        :: hess_rho_hn
	    
        double precision                  :: det_Dkn
        double precision, dimension (2)   :: vect_var
        double precision, dimension (2,2) :: Dkn, hess_phi_hn, tmp_mat
        integer(8)              :: k, nb_proc, chunk
       
		integer, external       :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM

        
!~         !$OMP PARALLEL 
!~         nb_proc = OMP_GET_NUM_THREADS()
!~         chunk = N_part / nb_proc +1
!~         !$OMP END PARALLEL 
!~ 		print *, 'nb_proc , chunk' , nb_proc , chunk
                
        hess_rho_hn=0.        
        hess_phi_hn=0.
        
!~ 		!$OMP PARALLEL DO shared(chunk) private(k, Dkn, det_Dkn, vect_var, hess_phi_hn) schedule(dynamic,chunk)
        do k=1,N_part 
			Dkn=0.
			hess_phi_hn=0.
			tmp_mat=0.
			vect_var=0.
            Dkn(1,1) = DD(1,k)
            Dkn(1,2) = DD(2,k)
            Dkn(2,1) = DD(3,k)
            Dkn(2,2) = DD(4,k)
            det_Dkn = Dkn(1,1) * Dkn(2,2) - Dkn(2,1) * Dkn(1,2)            
            vect_var = MATMUL(Dkn,X_coord - Xpart(:,k))            
            call hess_phi_2d_h(hess_phi_hn , vect_var(1) , vect_var(2) , hx_remap, hy_remap)            
            tmp_mat=MATMUL(TRANSPOSE(Dkn), hess_phi_hn)
            hess_rho_hn = hess_rho_hn + det_Dkn * M(k) * MATMUL(tmp_mat,Dkn)                             
        end do
!~ 		!$OMP END PARALLEL Do
    end subroutine

    
    !-------------------------------------------------------
    subroutine rec_ltp(Xpart, M, DD, N_part)
      integer, intent(in)                         :: N_part 
      double precision, dimension (2,N_part), intent(in)    :: Xpart
      double precision, dimension (N_part), intent(in)      :: M
      double precision, dimension (4, N_part), intent(in)   :: DD


      integer(8)   :: i, j, k, nb_proc, chunk
      integer, external :: OMP_GET_NUM_THREADS
      double precision       :: deltax, deltay, a, b, c, d, det_Dk
      double precision       :: rayonk, yj, xi
      integer(8)   :: qkx, qky, ik, jk
      double precision       :: phix, phjy

      deltax = abs(Xg(2) - Xg(1))
      deltay = abs(Yg(2) - Yg(1))
      ro = 0.      

!~       !$OMP PARALLEL 
      nb_proc = OMP_GET_NUM_THREADS()
      chunk = N_part / nb_proc +1        
!~       !$OMP END PARALLEL       
      
      
!~       !$OMP PARALLEL DO shared(chunk) private(i,j,k,a,b,c,d,det_Dk,rayonk,qkx,qky,ik,jk,yj,xi,phix,phjy) schedule(dynamic,chunk)      
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
                !ro(i,j) = ro(i,j) + M(k)/(hx*hy) * phix * phjy                
            end do
        end do
      end do
!~       !$OMP END PARALLEL DO

    end subroutine
    
    !-------------------------------------------------------
    ! routine version without call to phi (1D function) 
    subroutine rec_ltp_v2(Xpart, M, DD, N_part)
      integer, intent(in)                         :: N_part 
      double precision, dimension (2,N_part), intent(in)    :: Xpart
      double precision, dimension (N_part), intent(in)      :: M
      double precision, dimension (4, N_part), intent(in)   :: DD


      integer(8)   :: i, j, k
      double precision       :: deltax, deltay, a, b, c, d, det_Dk
      double precision       :: rayonk, yj, xi
      integer(8)   :: qkx, qky, ik, jk
      double precision       :: phi_hn
      double precision, dimension (2)         :: X_coord, vect_var
      double precision, dimension (2,2)         :: Dkn

      
      
      deltax = abs(Xg(2) - Xg(1))
      deltay = abs(Yg(2) - Yg(1))
      ro = 0.      
    
      do k=1,N_part
        a = DD(1, k)
        b = DD(2, k)   
        c = DD(3, k)
        d = DD(4, k)
        Dkn(1,1) = a
        Dkn(1,2) = b
        Dkn(2,1) = c
        Dkn(2,2) = d
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
                X_coord(1)=xi
                X_coord(2)=yj
                vect_var = MATMUL(Dkn,X_coord - Xpart(:,k))
                call phi_2d_h( phi_hn , vect_var(1) , vect_var(2) , hx, hy)
                ro(i,j) = ro(i,j) + M(k)*det_Dk * phi_hn
                !ro(i,j) = ro(i,j) + M(k)/(hx*hy) * phix * phjy                
            end do
        end do
      end do
    end subroutine

	subroutine rec_ltp_on_convolution_grid(rho_ltp_on_conv_grid, Xpart, M, DD, Xgrid_convol, Ygrid_convol, nq, N_part)
	  integer, intent(in)                         			:: N_part, nq
      double precision, dimension (2,N_part), intent(in)    :: Xpart
      double precision, dimension (N_part), intent(in)      :: M
      double precision, dimension (4, N_part), intent(in)   :: DD
      double precision, dimension (nq), intent(in)   		:: Xgrid_convol, Ygrid_convol
      double precision, dimension (nq,nq), intent(in)   	:: rho_ltp_on_conv_grid


      integer(8)   :: i, j, k
      double precision       :: deltax, deltay, a, b, c, d, det_Dk
      double precision       :: rayonk, yj, xi
      integer(8)   :: qkx, qky, ik, jk
      double precision       :: phi_hn
      double precision, dimension (2)         :: X_coord, vect_var
      double precision, dimension (2,2)         :: Dkn

      
      
      deltax = abs(Xgrid_convol(2) - Xgrid_convol(1))
      deltay = abs(Ygrid_convol(2) - Ygrid_convol(1))
      rho_ltp_on_conv_grid = 0.
    
      do k=1,N_part
        a = DD(1, k)
        b = DD(2, k)
        c = DD(3, k)
        d = DD(4, k)
        Dkn(1,1) = a
        Dkn(1,2) = b
        Dkn(2,1) = c
        Dkn(2,2) = d
        det_Dk = abs(a*d-b*c)
        if (det_Dk ==0) then
            print*, 'probleme de Dk non inversible'
        end if
        rayonk = max(hx,hy)/det_Dk*max(abs(d)+abs(b),abs(c)+abs(a))
        qkx = int(radius_phi * rayonk / deltax) + 1
        qky = int(radius_phi * rayonk / deltay) + 1
        ik = int((Xpart(1,k) - Xgrid_convol(1)) / deltax)
        jk = int((Xpart(2,k) - Ygrid_convol(1) ) / deltay)
        do j=max(1,jk- qky), min(jk+qky+1, nyg)
            yj = Ygrid_convol(j)  
            do i= max(ik-qkx,1), min(ik+qkx+1,nxg)
                xi = Xgrid_convol(i)
                X_coord(1)=xi
                X_coord(2)=yj
                vect_var = MATMUL(Dkn, X_coord - Xpart(:,k))
                call phi_2d_h( phi_hn, vect_var(1) , vect_var(2), hx, hy)
                rho_ltp_on_conv_grid(i,j) = rho_ltp_on_conv_grid(i,j) + M(k)*det_Dk * phi_hn
            end do
        end do
      end do
	end subroutine
    
    !-------------------------------------------------------
    subroutine rec_ltp_phi_support_non_compact(Xpart, M, DD, N_part)
      integer, intent(in)                         :: N_part 
      double precision, dimension (2,N_part), intent(in)    :: Xpart
      double precision, dimension (N_part), intent(in)      :: M
      double precision, dimension (4, N_part), intent(in)   :: DD


      integer(8)   :: i, j, k, nb_proc, chunk
      integer, external :: OMP_GET_NUM_THREADS
      double precision       :: deltax, deltay, a, b, c, d, det_Dk
      double precision       :: rayonk, yj, xi
      integer(8)   :: qkx, qky, ik, jk
      double precision       :: phix, phjy

      deltax = abs(Xg(2) - Xg(1))
      deltay = abs(Yg(2) - Yg(1))
      ro = 0.      
!      !$OMP PARALLEL 
      nb_proc = OMP_GET_NUM_THREADS()
      chunk = N_part / nb_proc +1        
!      !$OMP END PARALLEL       
      
!      !$OMP PARALLEL DO shared(chunk) private(i,j,k,a,b,c,d,det_Dk,rayonk,qkx,qky,ik,jk,yj,xi,phix,phjy) schedule(dynamic,chunk)      
      do k=1,N_part
        a = DD(1, k)
        b = DD(2, k)   
        c = DD(3, k)
        d = DD(4, k)
        det_Dk = abs(a*d-b*c)
        if (det_Dk ==0) then
            print*, 'probleme de Dk non inversible'
        end if
        
        do j=1,nyg
            yj = Yg(j)  
            do i=1,nxg
                xi = Xg(i)
                call phi(phix, (a*(Xpart(1,k)-xi)+b*(Xpart(2,k) - yj)) / hx)
                call phi(phjy, (c*(Xpart(1,k)-xi)+d*(Xpart(2,k) - yj)) / hy)
                ro(i,j) = ro(i,j) + M(k)*det_Dk/(hx*hy) * phix * phjy                
            end do
        end do
      end do
!      !$OMP END PARALLEL DO

    end subroutine

!-------------------------------------------------------
    subroutine rec_sansdepot(Xpart, M, DD, N_part)
      integer, intent(in)                         :: N_part 
      double precision, dimension (2,N_part), intent(in)    :: Xpart
      double precision, dimension (N_part), intent(in)      :: M
      double precision, dimension (4, N_part), intent(in)   :: DD


      integer(8)   :: i, j, k
      double precision       :: a, b, c, d, det_Dk
      double precision       :: phix, phjy, xi, yj
      double precision, dimension(2) :: Yk      

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
      double precision, dimension (2,N_part), intent(in)    :: Xpart
      double precision, dimension (N_part), intent(in)      :: M
      double precision, intent(in)                          :: epsi
      double precision, dimension(:,:), allocatable           :: depotx, depoty   

      integer(8)   :: i, j, k
      double precision       :: deltax, deltay
      double precision       :: yj, xi
      integer(8)   :: qx, qy, ik, jk
      double precision       :: phix, phjy
!
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

!!-------------------------------------------------------    
    subroutine rec_sp_at_coordinate(rho, Xpart, M, epsi, X_coord, N_part)    
        integer, intent(in)                         :: N_part
        double precision, dimension (2,N_part), intent(in)    :: Xpart
        double precision, dimension (2), intent(in)           :: X_coord
        double precision, dimension (N_part), intent(in)      :: M
        double precision, intent(in)                          :: epsi
        double precision, intent(out)                         :: rho

        integer(8)   :: k
        double precision       :: phi
		
        rho = 0.
        do k=1,N_part
            call phi_2d(phi, (X_coord(1)-Xpart(1,k))/epsi , (X_coord(2)-Xpart(2,k))/epsi )
            rho = rho + M(k)/(epsi*epsi) * phi
        end do
    end subroutine

!-------------------------------------------------------    
    subroutine rec_sp_grad_at_coordinate(grad_rho, Xpart, M, epsi, X_coord, N_part)    
        integer, intent(in)                         :: N_part
        double precision, dimension (2,N_part), intent(in)    :: Xpart
        double precision, dimension (2), intent(in)           :: X_coord
        double precision, dimension (N_part), intent(in)      :: M
        double precision, intent(in)                          :: epsi
        double precision, dimension (2), intent(out)          :: grad_rho

        integer(8)            :: k
        double precision, dimension (2) :: grad_phi
		
        grad_rho = 0.
        do k=1,N_part
            call grad_phi_2d(grad_phi, (X_coord(1)-Xpart(1,k))/epsi , (X_coord(2)-Xpart(2,k))/epsi )
            grad_rho = grad_rho + M(k)/(epsi*epsi*epsi) * grad_phi
        end do
    end subroutine

    subroutine rec_sp_hess_at_coordinate(hess_rho, Xpart, M, epsi, X_coord, N_part)    
        integer, intent(in)                         :: N_part
        double precision, dimension (2,N_part), intent(in)    :: Xpart
        double precision, dimension (2), intent(in)           :: X_coord
        double precision, dimension (N_part), intent(in)      :: M
        double precision, intent(in)                          :: epsi
        double precision, dimension (2,2), intent(out)        :: hess_rho

        integer(8)            :: k
        double precision, dimension (2,2) :: hess_phi
		
        hess_phi = 0.
        do k=1,N_part
            call hess_phi_2d(hess_phi, (X_coord(1)-Xpart(1,k))/epsi , (X_coord(2)-Xpart(2,k))/epsi )
            hess_phi = hess_phi + M(k)/(epsi*epsi*epsi*epsi) * hess_phi
        end do
    end subroutine


	!-------------------------------------------------------
	subroutine grad_convolution(epsi, X_coord, Xgrid_convol, Ygrid_convol, fct, quadrature_pts_convol, grad_convol, nq)
		integer, intent(in)                         		  :: nq
        double precision, dimension (2), intent(in)           :: X_coord        
        double precision, intent(in)                          :: epsi
        double precision, dimension (nq,nq), intent(in)       :: quadrature_pts_convol, fct
        double precision, dimension (nq), intent(in)        :: Xgrid_convol, Ygrid_convol        
        
        double precision, dimension (2), intent(out)          :: grad_convol
        
        
        integer(8)            :: i ,j                
        double precision, dimension (2), intent(in)       :: grad_phi
        
        U=0.
        grad_phi = 0.
        do i=1,nq
			do j=1,nq				
				call grad_phi_2d(grad_phi, (X_coord(1)-Xgrid_convol(i)/epsi , (X_coord(2)-Ygrid_convol(j))/epsi )
				grad_phi = (1./(epsi*epsi*epsi)) * grad_phi
				grad_convol += quadrature_pts_convol(i,j) * grad_phi * fct(i,j)
			end do
        end do
        
	end subroutine

  !-------------------------------------------------------
    subroutine diffusion_field_barenblatt_sp(Xpart, M, m_barenblatt,epsi,  diff_field, N_part)
        integer, intent(in)                         :: N_part, m_barenblatt
        double precision, dimension (2,N_part), intent(in)    :: Xpart
        double precision, dimension (N_part), intent(in)      :: M                
        double precision, dimension (2,N_part), intent(out)   :: diff_field
        double precision, intent(in)                          :: epsi
        
        double precision                  :: rho
        double precision, dimension (2)   :: grad_rho
        integer(8)     			          :: j

        integer(8)              :: nb_proc, chunk
        integer, external       :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM


		!$OMP PARALLEL 
        nb_proc = OMP_GET_NUM_THREADS()
        chunk = N_part / nb_proc +1
        !$OMP END PARALLEL 
        
        diff_field = 0.
        
		!$OMP PARALLEL DO shared(chunk) private(j, rho, grad_rho) schedule(dynamic,chunk)
        do j=1,N_part             			
            rho = 0.
            grad_rho = 0.            
            call rec_sp_at_coordinate(rho, Xpart, M, epsi, Xpart(:,j), N_part)
            call rec_sp_grad_at_coordinate(grad_rho, Xpart, M, epsi, Xpart(:,j), N_part)
            diff_field(1,j) = m_barenblatt * rho**(m_barenblatt-2)* grad_rho(1)
            diff_field(2,j) = m_barenblatt * rho**(m_barenblatt-2)* grad_rho(2)            
        end do
		!$OMP END PARALLEL DO    


    end subroutine

!-------------------------------------------------------
    subroutine diffusion_field(Xpart, M, DD, hx_remap, hy_remap, sigma, diff_field, N_part)
        integer, intent(in)                         :: N_part 
        double precision, dimension (2,N_part), intent(in)    :: Xpart
        double precision, dimension (N_part), intent(in)      :: M
        double precision, dimension (4, N_part), intent(in)   :: DD
        double precision, intent(in)                          :: hx_remap, hy_remap, sigma
        double precision, dimension (2,N_part), intent(out)   :: diff_field
	    
        double precision                  :: rho_hn
        double precision, dimension (2)   :: grad_rho_hn 
        integer(8)              :: j
        
        diff_field = 0.
        
!        !$OMP PARALLEL DO
        do j=1,N_part
            !ADD LINKED LISTS (OR NOT ?)
            rho_hn = 0.
            grad_rho_hn = 0.
            call rec_ltp_grad_at_coordinate(grad_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,j), N_part)
            call rec_ltp_at_coordinate(rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,j), N_part)
			!print*," j ,rho_hn , grad_rho_hn = "  ,j ,rho_hn, grad_rho_hn
            diff_field(1,j) = 2*grad_rho_hn(1) !barenblatt with m = 2
            diff_field(2,j) = 2*grad_rho_hn(2)
            !diff_field(1,j) = (grad_rho_hn(1) / rho_hn)
            !diff_field(2,j) = (grad_rho_hn(2) / rho_hn)
            !print *, "j , diff_field(1,j) , diff_field(2,j) = " , j , diff_field(1,j) , diff_field(2,j)
        end do
!        !$OMP END PARALLEL DO
        
        diff_field = sigma * diff_field
               
    end subroutine    
    
    
  !-------------------------------------------------------
    subroutine diffusion_field_barenblatt(Xpart, M, DD, m_barenblatt, hx_remap, hy_remap, diff_field, N_part)
        integer, intent(in)                         :: N_part, m_barenblatt
        double precision, dimension (2,N_part), intent(in)    :: Xpart
        double precision, dimension (N_part), intent(in)      :: M
        double precision, dimension (4, N_part), intent(in)   :: DD
        double precision, intent(in)                          :: hx_remap, hy_remap
        double precision, dimension (2,N_part), intent(out)   :: diff_field
	    
        double precision                  :: rho_hn
        double precision, dimension (2)   :: grad_rho_hn
        integer(8)              :: j        
        
        integer(8)              :: nb_proc, chunk
        integer, external       :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM


		!$OMP PARALLEL 
        nb_proc = OMP_GET_NUM_THREADS()
        chunk = N_part / nb_proc +1
        !$OMP END PARALLEL 
!~         print *, 'nb_proc , chunk' , nb_proc , chunk
        
        diff_field = 0.        

		if (m_barenblatt == 2) then
			!$OMP PARALLEL DO shared(chunk) private(j, grad_rho_hn) schedule(dynamic,chunk)
            do j=1,N_part                             
                grad_rho_hn = 0.                            
                call rec_ltp_grad_at_coordinate(grad_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,j), N_part)
                diff_field(1,j) = m_barenblatt * grad_rho_hn(1)
                diff_field(2,j) = m_barenblatt * grad_rho_hn(2)
            end do
			!$OMP END PARALLEL DO    

        else
            do j=1,N_part             
                rho_hn = 0.
                grad_rho_hn = 0.            
                call rec_ltp_at_coordinate(rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,j), N_part)
                call rec_ltp_grad_at_coordinate(grad_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,j), N_part)            
                diff_field(1,j) = m_barenblatt * rho_hn**(m_barenblatt-2)* grad_rho_hn(1)
                diff_field(2,j) = m_barenblatt * rho_hn**(m_barenblatt-2)* grad_rho_hn(2)                        
            end do
        end if
        
!        !$OMP END PARALLEL DO
    end subroutine

  
    subroutine diffusion_field_jacobian_barenblatt(diff_field_jac, Xpart, M, DD, m_barenblatt, hx_remap, hy_remap, N_part)
        integer, intent(in)                         :: N_part, m_barenblatt
        double precision, dimension (2,N_part), intent(in)    :: Xpart
        double precision, dimension (N_part), intent(in)      :: M
        double precision, dimension (4, N_part), intent(in)   :: DD
        double precision, intent(in)                          :: hx_remap, hy_remap
        double precision, dimension (4,N_part), intent(out)   :: diff_field_jac
	    
        double precision                  :: rho_hn
        double precision, dimension (2)   :: grad_rho_hn
        double precision, dimension (2,2) :: hess_rho_hn, diff_field_tmp, grad_mul
        integer(8)              :: j
        
        diff_field_tmp = 0.
        diff_field_jac = 0.
		!diff_field_jac(1,j) -> coef (1,1) de J_Ah
		!diff_field_jac(2,j) -> coef (1,2) de J_Ah
		!diff_field_jac(3,j) -> coef (2,1) de J_Ah
		!diff_field_jac(4,j) -> coef (2,2) de J_Ah
!        !$OMP PARALLEL DO
        if (m_barenblatt == 2) then
            do j=1,N_part
	            !ADD LINKED LISTS (OR NOT ?)
                rho_hn = 0.
                grad_rho_hn = 0.
                hess_rho_hn = 0.
	            !call rec_ltp_at_coordinate(rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,j), N_part)
	            !call rec_ltp_grad_at_coordinate(grad_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,j), N_part)	            	            
                call rec_ltp_hess_at_coordinate(hess_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,j), N_part)              
                diff_field_tmp = m_barenblatt * hess_rho_hn
                diff_field_jac(1,j) = diff_field_tmp(1,1)
                diff_field_jac(2,j) = diff_field_tmp(1,2)
                diff_field_jac(3,j) = diff_field_tmp(2,1)
                diff_field_jac(4,j) = diff_field_tmp(2,2)
            end do
        else
            do j=1,N_part 
	            !ADD LINKED LISTS (OR NOT ?)
                rho_hn = 0.
                grad_rho_hn = 0.
                hess_rho_hn = 0.
                call rec_ltp_at_coordinate(rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,j), N_part)
                call rec_ltp_grad_at_coordinate(grad_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,j), N_part)
                call rec_ltp_hess_at_coordinate(hess_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,j), N_part)
                print *, " rho_hn" , rho_hn
                print *, " grad_rho_hn" , grad_rho_hn
                print *, " hess_rho_hn" , hess_rho_hn
                grad_mul(1,1) = (grad_rho_hn(1))**2
                grad_mul(1,2) = grad_rho_hn(1)*grad_rho_hn(2)
                grad_mul(2,1) = grad_rho_hn(1)*grad_rho_hn(2)
                grad_mul(2,2) = (grad_rho_hn(2))**2
                diff_field_tmp = (m_barenblatt*(m_barenblatt-2)) * ((rho_hn)**(m_barenblatt-3)) * grad_mul&
                + m_barenblatt * rho_hn**(m_barenblatt-2) * hess_rho_hn
	            print *, " diff_field_tmp first part : " , m_barenblatt*(m_barenblatt-2) * rho_hn**(m_barenblatt-3) * grad_mul
	            print *, " diff_field_tmp grad_mul : " , grad_mul
	            print *, " diff_field_tmp second part : " , m_barenblatt * rho_hn**(m_barenblatt-2) * hess_rho_hn
	            print *, " diff_field_tmp" , diff_field_tmp                 
                diff_field_jac(1,j) = diff_field_tmp(1,1)
                diff_field_jac(2,j) = diff_field_tmp(1,2)
                diff_field_jac(3,j) = diff_field_tmp(2,1)
                diff_field_jac(4,j) = diff_field_tmp(2,2)
            end do
        end if
!        !$OMP END PARALLEL DO
    end subroutine    
    
    
    subroutine diffusion_field_barenblatt_convolution(Xpart, M , DD, m_barenblatt, hx_remap, hy_remap, epsi, &
    Xgrid_convol, Ygrid_convol, quadrature_pts_convol, diff_field, nq, N_part)
		integer, intent(in)                         :: N_part, nq, m_barenblatt
        double precision, dimension (2,N_part), intent(in)    :: Xpart
        double precision, dimension (N_part), intent(in)      :: M
        double precision, dimension (4, N_part), intent(in)   :: DD
        double precision, intent(in)                          :: hx_remap, hy_remap
             
        double precision, dimension (nq), intent(in)        :: Xgrid_convol, Ygrid_convol        
        double precision, dimension (nq,nq), intent(in)       :: quadrature_pts_convol
        
        double precision, dimension (2,N_part), intent(out)   :: diff_field
        
	    
        double precision, dimension (2)   :: result_grad_convol, X_coord
        integer(8)              :: j
        
		if (m_barenblatt == 2) then
			rho_ltp_conv_grid = rec_ltp_on_convolution_grid(Xpart, M, DD, Xgrid_convol, Ygrid_convol, nq, N_part)
            do j=1,N_part	                            
                call grad_convolution(result_grad_convol, epsi, X_coord, Xgrid_convol, Ygrid_convol, &
                rho_ltp_conv_grid, quadrature_pts_convol, nq, N_part)
                diff_field(1,j) = m_barenblatt * result_grad_convol(1)
                diff_field(2,j) = m_barenblatt * result_grad_convol(2)                
            end do
        else
			print *," m != 2 not done yet !! m_barenblatt = " , m_barenblatt
        end if
        
    end subroutine
! -------------------------------------------------------------------
! ================= subroutine d allocation/desallocation
!--------------------------------------------------------------------
    subroutine alloc_ro(dimx,dimy) 
      integer, intent(in)                    :: dimx
      integer, intent(in)                    :: dimy

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

end module calculsfor_rec_modf90



