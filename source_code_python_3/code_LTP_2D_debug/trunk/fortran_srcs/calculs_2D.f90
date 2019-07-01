module calculsfor_var_modf90

use calculsfor_rec_modf90, only : diffusion_field_barenblatt, diffusion_field_jacobian_barenblatt, &
								  rec_ltp_at_coordinate, rec_ltp_grad_at_coordinate, &
								  rec_ltp_hess_at_coordinate, &
                                  rec_ltp_on_convolution_grid, hess_convolution

implicit none
    integer :: flot
    double precision  :: Tmax
    double precision  :: alpha
contains


    !-------------------------------------------------------    
    subroutine set_flow_params(f,T,a)
        integer, intent(in) :: f
        double precision, intent(in)  :: T, a
        flot = f
        Tmax=  T
        alpha = a
    end subroutine

    !-------------------------------------------------------
    function adev(t,x,y) result(U)
        double precision, intent(in)      :: t, x, y
        double precision, dimension(2)    :: U
    
        double precision :: pi
        pi = 4.*atan(1.)
        if (flot == 2) then                
            !U(1) = -(1./pi)*cos(pi*t/Tmax) * (sin(pi*x))**2 * sin(2.*pi*y)
            !U(2) = (1./pi)*cos(pi*t/Tmax) * (sin(pi*y))**2 * sin(2.*pi*x)
            U(1) = - 2 * cos(pi*t/Tmax) * sin(pi*y)*cos(pi*y)*(sin(pi*x))**2
            U(2) = 2 * cos(pi*t/Tmax) * sin(pi*x)*cos(pi*x)*(sin(pi*y))**2
        else if (flot == 5) then
            U(1) = alpha * x
            U(2) = 0
        end if
                            
    end function adev

    !-------------------------------------------------------
    function Dadev(t,x,y) result(M)
        double precision, intent(in)     :: t, x, y
        double precision, dimension(2,2) :: M
    
        double precision :: pi
        pi = 4.*atan(1.)
        if (flot == 2) then
            !M(1,1) = - sin(pi*x) * cos(pi*x) * sin(2*pi*y)
            !M(2,1) = (sin(pi*y))**2 * cos(2*pi*x)
            !M(1,2) = - cos(2*pi*y) * (sin(pi*x))**2
            !M(2,2) =  sin(pi*y) * cos(pi*y) * sin(2*pi*x)
            !M = 2*cos(pi * t / Tmax) * M            
            
            M(1,1) = - 2*sin(pi*x) * cos(pi*x) * sin(pi*y) * cos(pi*y)
            M(2,1) = (sin(pi*y))**2 * ( -(sin(pi*x))**2 + (cos(pi*x))**2 )
            M(1,2) = (sin(pi*x))**2 * ( (sin(pi*y))**2 - (cos(pi*y))**2 )
            M(2,2) =  2*sin(pi*x) * cos(pi*x) * sin(pi*y) * cos(pi*y)
            M = 2*cos(pi * t / Tmax) * M    
        else if (flot == 5) then
            M(1,1) = alpha
            M(2,1) = 0
            M(1,2) = 0
            M(2,2) = 0
        end if
                
    end function Dadev

    !-------------------------------------------------------
    function compute_U_part(Xpart, t, N_part) result(Upart)
        integer, intent(in)                             :: N_part
        double precision, intent(in)                              :: t
        double precision, dimension (2,N_part), intent(in)        :: Xpart        
        double precision, dimension (2,N_part)                    :: Upart        
        
        integer :: k
        !$OMP PARALLEL DO private(k) schedule(dynamic)
        do k=1,N_part
            Upart(:,k) = adev(t, Xpart(1,k), Xpart(2,k))
        end do
        !$OMP END PARALLEL DO
    end function compute_U_part
        
    

    !-------------------------------------------------------
    ! approche exp(A), ou A est une matrice 2*2    
    function compute_expo_mat(A) result(exp_A)
        double precision, dimension(2,2), intent(in) :: A
        double precision, dimension(2,2)             :: exp_A
        
        double precision, dimension(2,2) :: I , B
        double precision                 :: det        
        det=0.
        B=0.
        I(1,1)=1.
        I(2,1)=0.
        I(1,2)=0.
        I(2,2)=1.        
        exp_A = I + A         
        det = exp_A(1,1) * exp_A(2,2) - exp_A(2,1) * exp_A(1,2)
        if (det > 0) then
            exp_A=exp_A            
        else
            B = A*A
            exp_A = exp_A + 0.5*B
            det = exp_A(1,1) * exp_A(2,2) - exp_A(2,1) * exp_A(1,2)
            if (det > 0) then
                exp_A=exp_A
            else
                B=B*A
                exp_A = exp_A + B/6.
            end if
        end if
    end function compute_expo_mat
    
    
    
    !-------------------------------------------------------    
    ! approximation of exp(A) with a prescribed order approximation
    function compute_expo_mat_order(A, N) result(exp_A)
        double precision, dimension(2,2), intent(in) :: A
        integer, intent(in)                :: N                
        double precision, dimension(2,2)             :: exp_A
        
        double precision, dimension(2,2) :: I , B
        double precision                 :: det        
        integer(8)               :: fact_k, k        
        I(1,1)=1.
        I(2,1)=0.
        I(1,2)=0.
        I(2,2)=1.
        B = A        
        exp_A = I + A
        fact_k=1
        do k=1,N
            fact_k = (k+1)*fact_k
            B = MATMUL(B,A)
            exp_A = exp_A + (1./fact_k) * B
!~             exp_A = exp_A + A^2 / 2. + A^3 / 6. + A^4 /4! + A^5
        end do
    end function compute_expo_mat_order


    
    function invMat2D(A) result(inv_A)
        double precision, dimension(2,2), intent(in) :: A
        double precision, dimension(2,2)             :: inv_A
        
        double precision    :: det_A        
        
        det_A = 0.
        inv_A(1,1) = A(2,2)
        inv_A(2,1) = -A(2,1)
        inv_A(1,2) = -A(1,2)
        inv_A(2,2) = A(1,1)
        det_A = A(1,1) * A(2,2) - A(2,1) * A(1,2)
        inv_A = (1./det_A) * inv_A        
        
    end function invMat2D
        
        
    !-------------------------------------------------------    
    subroutine update_d_diffusion(Xpart, DD, M, t, dt, time_scheme, hx_remap, hy_remap, &
     indice_max_norm_Dm1, Norm_inf_Dm1, Dout, N_part)
        integer, intent(in)                                       :: N_part, time_scheme
        double precision, dimension (2,N_part), intent(in)        :: Xpart        
        double precision, dimension (4, N_part), intent(in)       :: DD
        double precision, dimension (N_part), intent(in)          :: M
        double precision, intent(in)                              :: t, dt, hx_remap, hy_remap
        double precision, intent(out)                             :: Norm_inf_Dm1
        integer(8), intent(out)                         :: indice_max_norm_Dm1
        double precision, dimension (4, N_part), intent(out)      :: Dout
            
                
        integer(8)              :: k, invertible        
        double precision                  :: det_Dk, a, b, c, d, norm_inf_Dkm1, det_Jk, det_mat, new_dt
        double precision                  :: rho_hn
        double precision, dimension (2)   :: grad_rho_hn
        double precision, dimension(2,2)  :: Jk, Dk, k1, k2, invJk, Identity, hess_rho_hn, mat_grad_rho
        double precision, dimension (2,2) :: diff_field_jac               
        double precision, dimension(2,2)                  ::  max_hess_rho_hn
        double precision, dimension(2,2)                  ::  indice_max
        
        integer(8)              :: nb_proc, chunk
        integer, external       :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
        
        !$OMP PARALLEL 
        nb_proc = OMP_GET_NUM_THREADS()
        chunk = N_part / nb_proc +1
        !$OMP END PARALLEL 
!        print *, 'nb_proc , chunk' , nb_proc , chunk
        
        Dout=0.
        Norm_inf_Dm1 = 0.
        indice_max_norm_Dm1=0               
        max_hess_rho_hn = -100000.
        
        if (time_scheme == 1) then !Euler explicit (order 1)
            Identity(1,1) = 1.
            Identity(2,2) = 1.
            Identity(2,1) = 0.
            Identity(1,2) = 0.
            !$OMP PARALLEL DO shared(chunk) private(k,invJk,Dk,a,b,c,d,det_Dk,mat_grad_rho, &
            !$OMP& diff_field_jac,hess_rho_hn,grad_rho_hn,rho_hn) schedule(dynamic,chunk)
            do k=1,N_part
                !----- compute Jk ----- 
                hess_rho_hn = 0.
                invJk = 0.
                call rec_ltp_hess_at_coordinate(hess_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,k), N_part)
                call rec_ltp_grad_at_coordinate(grad_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,k), N_part)
                call rec_ltp_at_coordinate(rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,k), N_part)
                mat_grad_rho(1,1) = grad_rho_hn(1) * grad_rho_hn(1) 
                mat_grad_rho(1,2) = grad_rho_hn(1) * grad_rho_hn(2)
                mat_grad_rho(2,1) = grad_rho_hn(2) * grad_rho_hn(1)
                mat_grad_rho(2,2) = grad_rho_hn(2) * grad_rho_hn(2)
                diff_field_jac(1,1) = (rho_hn * hess_rho_hn(1,1) - mat_grad_rho(1,1)) / (rho_hn**2)
                diff_field_jac(1,2) = (rho_hn * hess_rho_hn(1,2) - mat_grad_rho(1,2)) / (rho_hn**2)
                diff_field_jac(2,1) = (rho_hn * hess_rho_hn(2,1) - mat_grad_rho(2,1)) / (rho_hn**2)
                diff_field_jac(2,2) = (rho_hn * hess_rho_hn(2,2) - mat_grad_rho(2,2)) / (rho_hn**2)
                invJk = compute_expo_mat_order(-dt*(- diff_field_jac), 10)
                !----- update D ----- 
                Dk(1,1)=DD(1,k)
                Dk(1,2)=DD(2,k)
                Dk(2,1)=DD(3,k)
                Dk(2,2)=DD(4,k)
                Dk = MATMUL(Dk,invJk)
                !-----  compute max infinity norm of Dk -----              
                a=DD(1,k)
                b=DD(2,k)
                c=DD(3,k)
                d=DD(4,k)
                det_Dk = abs(a*d-b*c)
                norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
                if (norm_inf_Dkm1 > Norm_inf_Dm1) then
                    Norm_inf_Dm1 = norm_inf_Dkm1
                    indice_max_norm_Dm1 = k
                end if
                !-----  update D -----
                Dout(1,k) = Dk(1,1)
                Dout(2,k) = Dk(1,2)
                Dout(3,k) = Dk(2,1)
                Dout(4,k) = Dk(2,2)
                det_Dk =Dk(1,1)*Dk(2,2) - Dk(1,2)*Dk(2,1)
            end do    
            !$OMP END PARALLEL DO         
        end if 
        
    end subroutine     
        
    !-------------------------------------------------------
    subroutine update_D(Xpart, DD, t, dt, indice_max_norm_Dm1, Norm_inf_Dm1, Dout,N_part)
        integer, intent(in)                             :: N_part
        double precision, dimension (2,N_part), intent(in)        :: Xpart        
        double precision, dimension (4, N_part), intent(in)       :: DD
        double precision, intent(in)                              :: t, dt
        double precision, intent(out)                             :: Norm_inf_Dm1
        integer(8), intent(out)                         :: indice_max_norm_Dm1
        double precision, dimension (4, N_part), intent(out)      :: Dout
                            
        integer(8)              :: k, nb_proc, chunk
        integer, external       :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
        double precision                  :: det_Dk, a, b, c, d, norm_inf_Dkm1
        double precision, dimension(2,2)  :: Akn, invJk, Dk
        !$OMP PARALLEL 
        nb_proc = OMP_GET_NUM_THREADS()
        chunk = N_part / nb_proc +1        
        !$OMP END PARALLEL 
        Dout=0.
        Norm_inf_Dm1 = 0.
        indice_max_norm_Dm1=0
        
        !print *, 'nb_proc , chunk' , nb_proc , chunk
        !$OMP PARALLEL DO shared(chunk) private(k,Akn,invJk,Dk,a,b,c,d,det_Dk) schedule(dynamic,chunk)
!        !$OMP PARALLEL DO private(k,Akn,invJk,Dk,a,b,c,d,det_Dk) schedule(dynamic)
        do k=1,N_part 
            !print *, 'OMP_GET_THREAD_NUM() , k : ' , OMP_GET_THREAD_NUM() , k
            !----- compute inverse Matrix of Jkn ----- 
            Akn = Dadev(t,Xpart(1,k),Xpart(2,k))            
            invJk = compute_expo_mat(-dt*Akn)         
            !print *, 'k, invJk : ' , k, invJk
            !----- update D ----- 
            Dk(1,1)=DD(1,k)
            Dk(1,2)=DD(2,k)
            Dk(2,1)=DD(3,k)
            Dk(2,2)=DD(4,k)
            Dk = Dk*invJk
            !-----  calcul de la norme infiny de l'inverse de D -----              
            a=DD(1,k)
            b=DD(2,k)
            c=DD(3,k)
            d=DD(4,k)
            det_Dk = abs(a*d-b*c)
            norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
            !print *,'k , norm_inf_Dkm1 , Norm_inf_Dm1 = ' , k, norm_inf_Dkm1, Norm_inf_Dm1            
            if (norm_inf_Dkm1 > Norm_inf_Dm1) then
                Norm_inf_Dm1 = norm_inf_Dkm1
                indice_max_norm_Dm1 = k
            end if
            !-----  update D -----
            Dout(1,k) = Dk(1,1)
            Dout(2,k) = Dk(1,2)
            Dout(3,k) = Dk(2,1)
            Dout(4,k) = Dk(2,2)
        end do
        !$OMP END PARALLEL DO 
    end subroutine        

    !-------------------------------------------------------
    subroutine update_D_schemes(Xpart, DD, t, dt,time_scheme, indice_max_norm_Dm1, Norm_inf_Dm1, Dout,N_part)
        integer, intent(in)                             :: N_part, time_scheme
        double precision, dimension (2,N_part), intent(in)        :: Xpart        
        double precision, dimension (4, N_part), intent(in)       :: DD
        double precision, intent(in)                              :: t, dt
        double precision, intent(out)                             :: Norm_inf_Dm1
        integer(8), intent(out)                         :: indice_max_norm_Dm1
        double precision, dimension (4, N_part), intent(out)      :: Dout
                            
        integer(8)              :: k, nb_proc, chunk
        integer, external       :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
        double precision                  :: det_Dk, a, b, c, d, norm_inf_Dkm1, det_Jk
        double precision, dimension(2,2)  :: Akn, Jk, Dk, k1, k2, invJk
        !!$OMP PARALLEL 
        nb_proc = OMP_GET_NUM_THREADS()
        chunk = N_part / nb_proc +1        
        !!$OMP END PARALLEL 
        Dout=0.
        Norm_inf_Dm1 = 0.
        indice_max_norm_Dm1=0        
        !print *, 'nb_proc , chunk' , nb_proc , chunk
        !!$OMP PARALLEL DO shared(chunk) private(k,Akn,invJk,Dk,a,b,c,d,det_Dk) schedule(dynamic,chunk)
        if (time_scheme == 1) then!Euler explicit
            do k=1,N_part
                k1(1,1) = 1
                k1(2,2) = 1
                k1(2,1) = 0
                k1(1,2) = 0
                !----- compute Jk ----- 
                Jk = k1 + dt * Dadev(t,Xpart(1,k),Xpart(2,k))
                !----- compute matrix inverse of Jk ----- 
                invJk(1,1) = Jk(2,2)
                invJk(2,1) = -Jk(2,1)
                invJk(1,2) = -Jk(1,2)
                invJk(2,2) = Jk(1,1)
                det_Jk = Jk(1,1) * Jk(2,2) - Jk(2,1) * Jk(1,2)
                invJk = (1/det_Jk) * invJk
                !----- update D ----- 
                Dk(1,1)=DD(1,k)
                Dk(1,2)=DD(2,k)
                Dk(2,1)=DD(3,k)
                Dk(2,2)=DD(4,k)
                Dk = Dk*invJk
                !-----  calcul de la norme infiny de l'inverse de D -----              
                a=DD(1,k)
                b=DD(2,k)
                c=DD(3,k)
                d=DD(4,k)
                det_Dk = abs(a*d-b*c)
                norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
                !print *,'k , norm_inf_Dkm1 , Norm_inf_Dm1 = ' , k, norm_inf_Dkm1, Norm_inf_Dm1            
                if (norm_inf_Dkm1 > Norm_inf_Dm1) then
                    Norm_inf_Dm1 = norm_inf_Dkm1
                    indice_max_norm_Dm1 = k
                end if
                !-----  update D -----
                Dout(1,k) = Dk(1,1)
                Dout(2,k) = Dk(1,2)
                Dout(3,k) = Dk(2,1)
                Dout(4,k) = Dk(2,2)
            end do
        end if

    end subroutine        


    !-------------------------------------------------------
    !!! barenblatt with m = 2 for now
    subroutine update_D_schemes_diffusion_barenblatt(Xpart, DD, M, t, dt, time_scheme, m_barenblatt, hx_remap, hy_remap, &
     indice_max_norm_Dm1, Norm_inf_Dm1, Dout, N_part)
        integer, intent(in)                             :: N_part, time_scheme, m_barenblatt
        double precision, dimension (2,N_part), intent(in)        :: Xpart        
        double precision, dimension (4, N_part), intent(in)       :: DD
        double precision, dimension (N_part), intent(in)          :: M
        double precision, intent(in)                              :: t, dt, hx_remap, hy_remap
        double precision, intent(out)                             :: Norm_inf_Dm1
        integer(8), intent(out)                         :: indice_max_norm_Dm1
        double precision, dimension (4, N_part), intent(out)      :: Dout
                            
        integer(8)              :: k, invertible        
        double precision                  :: det_Dk, a, b, c, d, norm_inf_Dkm1, det_Jk, det_mat, new_dt
        double precision, dimension(2,2)  :: Jk, Dk, k1, k2, invJk, Identity, hess_rho_hn, diff_field_tmp, mat
        double precision, dimension (4, N_part) :: diff_field_jac               
        double precision, dimension(2,2)                  ::  max_hess_rho_hn
        double precision, dimension(2,2)                  ::  indice_max
        
        integer(8)              :: nb_proc, chunk
        integer, external       :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
        
        !$OMP PARALLEL 
        nb_proc = OMP_GET_NUM_THREADS()
        chunk = N_part / nb_proc +1
        !$OMP END PARALLEL 
!~         print *, 'nb_proc , chunk' , nb_proc , chunk
        
        Dout=0.
        Norm_inf_Dm1 = 0.
        indice_max_norm_Dm1=0               
        max_hess_rho_hn = -100000.
        
        if (time_scheme == 1) then!Euler explicit (order 1)
            Identity(1,1) = 1.
            Identity(2,2) = 1.
            Identity(2,1) = 0.
            Identity(1,2) = 0.
            !$OMP PARALLEL DO shared(chunk) private(k,invJk,Dk,a,b,c,d, det_Dk, hess_rho_hn) schedule(dynamic,chunk)
            do k=1,N_part
                !----- compute Jk ----- 
                hess_rho_hn = 0.
                invJk = 0.
                call rec_ltp_hess_at_coordinate(hess_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,k), N_part)
                invJk = compute_expo_mat_order(-dt*(-m_barenblatt*hess_rho_hn), 10)
                !----- update D ----- 
                Dk(1,1)=DD(1,k)
                Dk(1,2)=DD(2,k)
                Dk(2,1)=DD(3,k)
                Dk(2,2)=DD(4,k)
                Dk = MATMUL(Dk,invJk)
                !-----  calcul de la norme infiny de l'inverse de D -----              
                a=DD(1,k)
                b=DD(2,k)
                c=DD(3,k)
                d=DD(4,k)
                det_Dk = abs(a*d-b*c)
                norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
                if (norm_inf_Dkm1 > Norm_inf_Dm1) then
                    Norm_inf_Dm1 = norm_inf_Dkm1
                    indice_max_norm_Dm1 = k
                end if
                !-----  update D -----
                Dout(1,k) = Dk(1,1)
                Dout(2,k) = Dk(1,2)
                Dout(3,k) = Dk(2,1)
                Dout(4,k) = Dk(2,2)
                det_Dk =Dk(1,1)*Dk(2,2) - Dk(1,2)*Dk(2,1)
!~                 if (abs(det_Dk) <= 0.0001)  then
!~                     print *, " k , det_Dk , Dk " , k , det_Dk , Dk(1,1) , Dk(1,2) , Dk(2,1) , Dk(2,2)
!~                 end if
            end do    
            !$OMP END PARALLEL DO    
        else if (time_scheme == 2) then!midpoint scheme 
            Identity(1,1) = 1.
            Identity(2,2) = 1.
            Identity(2,1) = 0.
            Identity(1,2) = 0.
        else if (time_scheme == 3) then!Euler implicit
            Identity(1,1) = 1.
            Identity(2,2) = 1.
            Identity(2,1) = 0.
            Identity(1,2) = 0.
            do k=1,N_part
                !----- compute Jk ----- 
                call rec_ltp_hess_at_coordinate(hess_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,k), N_part)
                invertible = 0
                new_dt=dt
                do while  (invertible == 0)
                    mat = Identity + new_dt * m_barenblatt * hess_rho_hn
                    det_mat = mat(1,1) * mat(2,2) - mat(2,1) * mat(1,2)
                    !if ((abs(det_mat) - 0.00000001) <= 0) then
                    if (det_mat <= 0) then
                        new_dt = new_dt/2.
                        print *, "k, new dt = " , k, new_dt
                    else
                        invertible = 1
                    end if                    
                end do
                
                invJk = mat
                !print *," diff_field_jac_square : " , diff_field_jac_square
                !print *," Xpart(:,k) : " , Xpart(1,k) , Xpart(2,k)
                !----- compute matrix inverse of Jk ----- 
!                invJk(1,1) = Jk(2,2)
!                invJk(2,1) = -Jk(2,1)
!                invJk(1,2) = -Jk(1,2)
!                invJk(2,2) = Jk(1,1)
!                det_Jk = Jk(1,1) * Jk(2,2) - Jk(2,1) * Jk(1,2)
!                invJk = (1./det_Jk) * invJk
                !invJk = invMat2D(Jk)
                !----- update D ----- 
                Dk(1,1)=DD(1,k)
                Dk(1,2)=DD(2,k)
                Dk(2,1)=DD(3,k)
                Dk(2,2)=DD(4,k)
                Dk = MATMUL(Dk,invJk)
                !-----  calcul de la norme infiny de l'inverse de D -----              
                a=DD(1,k)
                b=DD(2,k)
                c=DD(3,k)
                d=DD(4,k)
                det_Dk = abs(a*d-b*c)
                norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
                !print *,'k , norm_inf_Dkm1 , Norm_inf_Dm1 = ' , k, norm_inf_Dkm1, Norm_inf_Dm1            
                if (norm_inf_Dkm1 > Norm_inf_Dm1) then
                    Norm_inf_Dm1 = norm_inf_Dkm1
                    indice_max_norm_Dm1 = k
                end if
                !-----  update D -----
                Dout(1,k) = Dk(1,1)
                Dout(2,k) = Dk(1,2)
                Dout(3,k) = Dk(2,1)
                Dout(4,k) = Dk(2,2)
            end do            
        end if 
        
!~         print *,"max_hess_rho_hn(1,1), indice_max(1,1) = " , max_hess_rho_hn(1,1), indice_max(1,1)
!~         print *,"max_hess_rho_hn(2,1), indice_max(2,1) = " , max_hess_rho_hn(2,1), indice_max(2,1)
!~         print *,"max_hess_rho_hn(1,2), indice_max(1,2) = " , max_hess_rho_hn(1,2), indice_max(1,2)
!~         print *,"max_hess_rho_hn(2,2), indice_max(2,2) = " , max_hess_rho_hn(2,2), indice_max(2,2)
    end subroutine       
    
    subroutine update_d_hybrid_diffusion_barenblatt(Xpart, DD, M, t, dt, m_barenblatt, &
    epsilon_SP, hx_remap, hy_remap, indice_max_norm_Dm1, Norm_inf_Dm1, Dout, N_part)
        integer, intent(in)                             :: N_part, m_barenblatt
        double precision, dimension (2,N_part), intent(in)        :: Xpart        
        double precision, dimension (4, N_part), intent(in)       :: DD
        double precision, dimension (N_part), intent(in)          :: M
        double precision, intent(in)                              :: t, dt, hx_remap, hy_remap, epsilon_SP
!~         double precision, intent(out)                             :: Norm_inf_Dm1
        double precision, intent(out)                             :: Norm_inf_Dm1
        integer(8), intent(out)                         :: indice_max_norm_Dm1
        double precision, dimension (4, N_part), intent(out)      :: Dout
                            
        integer(8)              :: k, nb_proc, chunk, invertible        
        double precision                  :: det_Dk, a, b, c, d, norm_inf_Dkm1, det_Jk, det_mat, new_dt
        double precision, dimension(2,2)  :: Akn, Jk, Dk, k1, k2, invJk, Identity, hess_rho_hn, diff_field_tmp, mat
        double precision, dimension (4, N_part) :: diff_field_jac
        Identity(1,1) = 1.
        Identity(2,2) = 1.
        Identity(2,1) = 0.
        Identity(1,2) = 0.
        do k=1,N_part
            !----- compute Jk ----- 
            call rec_ltp_hess_at_coordinate(hess_rho_hn, Xpart, M, DD, epsilon_SP, epsilon_SP, Xpart(:,k), N_part)
!~                 invJk = compute_expo_mat(-dt*(-m_barenblatt*hess_rho_hn))
            invJk = compute_expo_mat_order(-dt*(-m_barenblatt*hess_rho_hn), 10)
            !----- update D ----- 
            Dk(1,1)=DD(1,k)
            Dk(1,2)=DD(2,k)
            Dk(2,1)=DD(3,k)
            Dk(2,2)=DD(4,k)
            Dk = MATMUL(Dk,invJk)
            !-----  calcul de la norme infiny de l'inverse de D -----              
            a=DD(1,k)
            b=DD(2,k)
            c=DD(3,k)
            d=DD(4,k)
            det_Dk = abs(a*d-b*c)
            norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
            !print *,'k , norm_inf_Dkm1 , Norm_inf_Dm1 = ' , k, norm_inf_Dkm1, Norm_inf_Dm1            
            if (norm_inf_Dkm1 > Norm_inf_Dm1) then
                Norm_inf_Dm1 = norm_inf_Dkm1
                indice_max_norm_Dm1 = k
            end if
            !-----  update D -----
            Dout(1,k) = Dk(1,1)
            Dout(2,k) = Dk(1,2)
            Dout(3,k) = Dk(2,1)
            Dout(4,k) = Dk(2,2)
            det_Dk =Dk(1,1)*Dk(2,2) - Dk(1,2)*Dk(2,1)
            if (abs(det_Dk) <= 0.0001)  then
                print *, " k , det_Dk , Dk " , k , det_Dk , Dk(1,1) , Dk(1,2) , Dk(2,1) , Dk(2,2)
            end if
        end do                    
    end subroutine 

    !-------------------------------------------------------
    subroutine update_d_convolution_method_diffusion_barenblatt(Xpart, DD, M, t, dt, m_barenblatt, &
               epsilon_SP, hx_remap, hy_remap, Xgrid_convol, Ygrid_convol, quadrature_pts_convol, &
               indice_max_norm_Dm1, Norm_inf_Dm1, Dout, nqx, nqy, N_part)
        integer, intent(in)                                       :: N_part, m_barenblatt, nqx, nqy
        double precision, dimension (2,N_part), intent(in)        :: Xpart        
        double precision, dimension (4, N_part), intent(in)       :: DD
        double precision, dimension (N_part), intent(in)          :: M
        double precision, intent(in)                              :: t, dt, hx_remap, hy_remap, epsilon_SP        
        double precision, dimension (nqx), intent(in)             :: Xgrid_convol        
        double precision, dimension (nqy), intent(in)             :: Ygrid_convol
        double precision, dimension (nqx,nqy), intent(in)         :: quadrature_pts_convol
        
        integer(8), intent(out)                                   :: indice_max_norm_Dm1
        double precision, intent(out)                             :: Norm_inf_Dm1        
        double precision, dimension (4, N_part), intent(out)      :: Dout
                            
        integer(8)                              :: k, nb_proc, chunk, invertible        
        double precision                        :: det_Dk, a, b, c, d, norm_inf_Dkm1, &
                                                   det_Jk, det_mat, new_dt
        double precision, dimension(2,2)        :: Akn, Jk, Dk, k1, k2, invJk, Identity, &
                                                   hess_convol, diff_field_tmp, mat
        double precision, dimension (4, N_part) :: diff_field_jac
        double precision, dimension (nqx,nqy)   :: rho_ltp_conv_grid
        
        Identity(1,1) = 1.
        Identity(2,2) = 1.
        Identity(2,1) = 0.
        Identity(1,2) = 0.
        rho_ltp_conv_grid=0.
        hess_convol = 0.
        call rec_ltp_on_convolution_grid(rho_ltp_conv_grid, Xpart, M, DD, Xgrid_convol, Ygrid_convol, &
                 hx_remap, hy_remap, nqx, nqy, N_part)
        do k=1,N_part
            !----- compute Jk ----- 
!~             call rec_ltp_hess_at_coordinate(hess_rho_hn, Xpart, M, DD, epsilon_SP, epsilon_SP, Xpart(:,k), N_part)
            
            call hess_convolution(hess_convol, epsilon_SP, Xpart(:,k), Xgrid_convol, Ygrid_convol, rho_ltp_conv_grid, &
            quadrature_pts_convol, nqx, nqy)
!~                 invJk = compute_expo_mat(-dt*(-m_barenblatt*hess_rho_hn))
            invJk = compute_expo_mat_order(-dt*(-m_barenblatt*hess_convol), 10)
            !----- update D ----- 
            Dk(1,1)=DD(1,k)
            Dk(1,2)=DD(2,k)
            Dk(2,1)=DD(3,k)
            Dk(2,2)=DD(4,k)
            Dk = MATMUL(Dk,invJk)
            !-----  calcul de la norme infiny de l'inverse de D -----              
            a=DD(1,k)
            b=DD(2,k)
            c=DD(3,k)
            d=DD(4,k)
            det_Dk = abs(a*d-b*c)
            norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
            !print *,'k , norm_inf_Dkm1 , Norm_inf_Dm1 = ' , k, norm_inf_Dkm1, Norm_inf_Dm1            
            if (norm_inf_Dkm1 > Norm_inf_Dm1) then
                Norm_inf_Dm1 = norm_inf_Dkm1
                indice_max_norm_Dm1 = k
            end if
            !-----  update D -----
            Dout(1,k) = Dk(1,1)
            Dout(2,k) = Dk(1,2)
            Dout(3,k) = Dk(2,1)
            Dout(4,k) = Dk(2,2)
            det_Dk =Dk(1,1)*Dk(2,2) - Dk(1,2)*Dk(2,1)
            if (abs(det_Dk) <= 0.0001)  then
                print *, " k , det_Dk , Dk " , k , det_Dk , Dk(1,1) , Dk(1,2) , Dk(2,1) , Dk(2,2)
            end if
        end do                    
    end subroutine 

!~     subroutine update_d_barenblatt_mid_scheme_v2(Xpart, Xpart_middle, DD, M, t, dt, time_scheme, m_barenblatt, hx_remap, &
!~     hy_remap, D_middle, indice_max_norm_Dm1, Norm_inf_Dm1, Dout, N_part)
!~         integer, intent(in)                             :: N_part, time_scheme, m_barenblatt
!~         double precision, dimension (2,N_part), intent(in)        :: Xpart, Xpart_middle        
!~         double precision, dimension (4, N_part), intent(in)       :: DD
!~         double precision, dimension (N_part), intent(in)          :: M
!~         double precision, intent(in)                              :: t, dt, hx_remap, hy_remap
!~         double precision, intent(out)                             :: Norm_inf_Dm1
!~         integer(8), intent(out)                         :: indice_max_norm_Dm1
!~         double precision, dimension (4, N_part), intent(out)      :: D_middle, Dout
        
                            
!~         integer(8)              :: k, invertible        
!~         double precision                  :: det_Dk, a, b, c, d, norm_inf_Dkm1, det_Jk, det_mat, new_dt
!~         double precision, dimension(2,2)  :: Akn, Jk, Dk, k1, k2, invJk, Identity, hess_rho_hn
!~         double precision, dimension(2,2)  :: hess_rho_hn_middle, diff_field_tmp, mat, Jk_mid
!~         double precision, dimension (4, N_part) :: diff_field_jac, J_mid
        
!~         integer(8)              :: nb_proc, chunk
!~         integer, external       :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
        
!~         !$OMP PARALLEL 
!~         nb_proc = OMP_GET_NUM_THREADS()
!~         chunk = N_part / nb_proc +1
!~         !$OMP END PARALLEL 
        
!~         Dout=0.
!~         Norm_inf_Dm1 = 0.
!~         indice_max_norm_Dm1=0
        
!~         Identity(1,1) = 1.
!~         Identity(2,2) = 1.
!~         Identity(2,1) = 0.
!~         Identity(1,2) = 0.   
                
!~         ! Compute Dk at intermediate time step n + 1/2
!~         !$OMP PARALLEL DO shared(chunk) private(k,invJk,Dk,a,b,c,d, det_Dk, hess_rho_hn) schedule(dynamic,chunk)
!~         do k=1,N_part
!~             !----- compute Jk ----- 
!~             call rec_ltp_hess_at_coordinate(hess_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,k), N_part)                
!~             invJk = compute_expo_mat(-(dt/2.)*(-m_barenblatt*hess_rho_hn))         
!~             Jk_mid = compute_expo_mat_order((dt/2.)*(-m_barenblatt*hess_rho_hn), 10)
!~             !----- update D ----- 
!~             J_mid(0,k) = Jk_mid(0,0)
!~             J_mid(1,k) = Jk_mid(0,1)    
!~             J_mid(2,k) = Jk_mid(1,0)
!~             J_mid(3,k) = Jk_mid(1,1)
!~             Dk(1,1)=DD(1,k)
!~             Dk(1,2)=DD(2,k)
!~             Dk(2,1)=DD(3,k)
!~             Dk(2,2)=DD(4,k)
!~             Dk = MATMUL(Dk,invJk)
!~             !-----  update D -----
!~             D_middle(1,k) = Dk(1,1)
!~             D_middle(2,k) = Dk(1,2)
!~             D_middle(3,k) = Dk(2,1)
!~             D_middle(4,k) = Dk(2,2)
!~             det_Dk =Dk(1,1)*Dk(2,2) - Dk(1,2)*Dk(2,1) 
!~             if (abs(det_Dk) <= 0.0001)  then
!~                 print *, " k , det_Dk , Dk " , k , det_Dk , Dk(1,1) , Dk(1,2) , Dk(2,1) , Dk(2,2)                
!~             end if
!~         end do 
!~         !$OMP END PARALLEL DO

!~         ! compute Dkn+1               
!~         !$OMP PARALLEL DO shared(chunk) private(k,invJk,Dk,a,b,c,d, det_Dk, hess_rho_hn_middle) schedule(dynamic,chunk)
!~         do k=1,N_part
!~             !----- compute Jk ----- 
!~             Jk_mid(0,0) = J_mid(0,k)
!~             Jk_mid(0,1) = J_mid(1,k)
!~             Jk_mid(1,0) = J_mid(2,k)
!~             Jk_mid(1,1) = J_mid(3,k)
!~             call rec_ltp_hess_at_coordinate(hess_rho_hn_middle, Xpart_middle, M, D_middle, hx_remap, hy_remap, &
!~             Xpart_middle(:,k), N_part)            
!~             invJk = compute_expo_mat_order(-dt*(-m_barenblatt*MATMUL(hess_rho_hn_middle*Jk_mid) ),10)

!~             !----- update D ----- 
!~             Dk(1,1)=DD(1,k)
!~             Dk(1,2)=DD(2,k)
!~             Dk(2,1)=DD(3,k)
!~             Dk(2,2)=DD(4,k)
!~             Dk = MATMUL(Dk,invJk)
!~             !-----  calcul de la norme infiny de l'inverse de D -----              
!~             a=DD(1,k)
!~             b=DD(2,k)
!~             c=DD(3,k)
!~             d=DD(4,k)
!~             det_Dk = abs(a*d-b*c)
!~             norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
!~             !print *,'k , norm_inf_Dkm1 , Norm_inf_Dm1 = ' , k, norm_inf_Dkm1, Norm_inf_Dm1            
!~             if (norm_inf_Dkm1 > Norm_inf_Dm1) then
!~                 Norm_inf_Dm1 = norm_inf_Dkm1
!~                 indice_max_norm_Dm1 = k
!~             end if
!~             !-----  update D -----
!~             Dout(1,k) = Dk(1,1)
!~             Dout(2,k) = Dk(1,2)
!~             Dout(3,k) = Dk(2,1)
!~             Dout(4,k) = Dk(2,2)
!~             det_Dk =Dk(1,1)*Dk(2,2) - Dk(1,2)*Dk(2,1) 
!~             if (abs(det_Dk) <= 0.0001)  then
!~                 print *, " k , det_Dk , Dk " , k , det_Dk , Dk(1,1) , Dk(1,2) , Dk(2,1) , Dk(2,2)                
!~             end if
!~         end do                
!~         !$OMP END PARALLEL DO  
!~     end subroutine   
    
    !~         print *, "fortran : Norm_inf_Dm1, indice_max_norm_Dm1 = " , Norm_inf_Dm1, indice_max_norm_Dm1
!~         do k=1,N_part
!~             print *, k, Dout(1,k) , Dout(2,k), Dout(3,k), Dout(4,k)
!~             print *, k, D_middle(1,k) , D_middle(2,k), D_middle(3,k), D_middle(4,k)
!~         end do
    
!~     subroutine update_positions_and_deformations_rk4(Xpart, DD, M, t, dt, m_barenblatt, hx_remap, &
!~     hy_remap, X_new, indice_max_norm_Dm1, Norm_inf_Dm1, D_new, N_part)
!~         integer, intent(in)                                       :: N_part, m_barenblatt
!~         double precision, dimension (2,N_part), intent(in)        :: Xpart
!~         double precision, dimension (4, N_part), intent(in)       :: DD
!~         double precision, dimension (N_part), intent(in)          :: M
!~         double precision, intent(in)                              :: t, dt, hx_remap, hy_remap
!~         double precision, intent(out)                             :: Norm_inf_Dm1
!~         integer(8), intent(out)                                   :: indice_max_norm_Dm1
!~         double precision, dimension (4, N_part), intent(out)      :: D_new
!~         double precision, dimension (2, N_part), intent(out)      :: X_new
!~                                     
!~         integer(8)              :: k, invertible        
!~         double precision                  :: det_Dk, a, b, c, d, norm_inf_Dkm1, det_Jk, det_mat, new_dt
!~         double precision, dimension(2,2)  :: Akn, Jk, Dk, k1, k2, invJk, Identity, hess_rho_hn
!~         double precision, dimension(2,2)  :: hess_rho_hn_middle, diff_field_tmp, mat
!~         double precision, dimension (4, N_part) :: diff_field_jac
!~         
!~         Identity(1,1) = 1.
!~         Identity(2,2) = 1.
!~         Identity(2,1) = 0.
!~         Identity(1,2) = 0.   
!~         
!~         ! COMPUTE X1_mid and D1_mid
!~         do k=1,N_part            
!~             call diffusion_field_barenblatt(U, Xpart, M, DD, m_barenblatt, hx_remap, hy_remap, N_part)
!~             U = -U
!~             X1_new(:,k) = Xpart(:,k) + dt * U(:,k)
!~             X1_mid(:,k) = (1./2.) * (Xpart(:,k) + X1_new(:,k))
!~             call rec_ltp_hess_at_coordinate(hess_rho_hn, Xpart, M, DD, hx_remap, hy_remap, Xpart(:,k), N_part)
!~             invJk = compute_expo_mat_order(-dt*(-m_barenblatt*hess_rho_hn),10)
!~             DU(1,k) = -m_barenblatt*hess_rho_hn(1,1)
!~             DU(2,k) = -m_barenblatt*hess_rho_hn(1,2)
!~             DU(3,k) = -m_barenblatt*hess_rho_hn(2,1)
!~             DU(4,k) = -m_barenblatt*hess_rho_hn(2,2)
!~             Dk(1,1)=DD(1,k)
!~             Dk(1,2)=DD(2,k)
!~             Dk(2,1)=DD(3,k)
!~             Dk(2,2)=DD(4,k)
!~             Dk = MATMUL(Dk,invJk)
!~             D1_new(1,k) =  Dk(1,1)
!~             D1_new(2,k) =  Dk(1,2)
!~             D1_new(3,k) =  Dk(2,1)
!~             D1_new(4,k) =  Dk(2,2)
!~             D1_mid(:,k) = (1./2.) * (DD(:,k) + D1_new(:,k))
!~         end do
!~            
!~         ! COMPUTE X2_mid and D2_mid from U1_mid
!~         do k=1,N_part            
!~             call diffusion_field_barenblatt(U1_mid, X1_mid, M, D1_mid, m_barenblatt, hx_remap, hy_remap, N_part)
!~             U1_mid = -U1_mid
!~             X2_new(:,k) = Xpart(:,k) + dt * U1_mid(:,k)
!~             X2_mid(:,k) = (1./2.) * (Xpart(:,k) + X2_new(:,k))
!~             call rec_ltp_hess_at_coordinate(hess_rho_hn, X1_mid, M, D1_mid, hx_remap, hy_remap, X1_mid(:,k), N_part)
!~             invJk = compute_expo_mat_order(-dt*(-m_barenblatt*hess_rho_hn),10)
!~             DU_mid(1,k) = -m_barenblatt*hess_rho_hn(1,1)
!~             DU(2,k) = -m_barenblatt*hess_rho_hn(1,2)
!~             DU(3,k) = -m_barenblatt*hess_rho_hn(2,1)
!~             DU(4,k) = -m_barenblatt*hess_rho_hn(2,2)
!~             Dk(1,1)=DD(1,k)
!~             Dk(1,2)=DD(2,k)
!~             Dk(2,1)=DD(3,k)
!~             Dk(2,2)=DD(4,k)
!~             Dk = MATMUL(Dk,invJk)
!~             D2_new(1,k) =  Dk(1,1)
!~             D2_new(2,k) =  Dk(1,2)
!~             D2_new(3,k) =  Dk(2,1)
!~             D2_new(4,k) =  Dk(2,2)
!~             D2_mid(:,k) = (1./2.) * (DD(:,k) + D2_new(:,k))
!~         end do
!~ 
!~         ! COMPUTE X3_mid and D3_mid from U2_mid
!~         do k=1,N_part
!~             call diffusion_field_barenblatt(U2_mid, X2_mid, M, D2_mid, m_barenblatt, hx_remap, hy_remap, N_part)
!~             U2_mid = -U2_mid
!~             X3_new(:,k) = Xpart(:,k) + dt * U2_mid(:,k)
!~             X3_mid(:,k) = (1./2.) * (Xpart(:,k) + X3_new(:,k))
!~             call rec_ltp_hess_at_coordinate(hess_rho_hn, X2_mid, M, D2_mid, hx_remap, hy_remap, X2_mid(:,k), N_part)
!~             invJk = compute_expo_mat_order(-dt*(-m_barenblatt*hess_rho_hn),10)
!~             DU2_mid(1,k) = invJk(1,1)
!~             DU2_mid(2,k) = invJk(1,2)
!~             DU2_mid(3,k) = invJk(2,1)
!~             DU2_mid(4,k) = invJk(2,2)
!~             Dk(1,1)=DD(1,k)
!~             Dk(1,2)=DD(2,k)
!~             Dk(2,1)=DD(3,k)
!~             Dk(2,2)=DD(4,k)
!~             Dk = MATMUL(Dk,invJk)
!~             D3_new(1,k) =  Dk(1,1)
!~             D3_new(2,k) =  Dk(1,2)
!~             D3_new(3,k) =  Dk(2,1)
!~             D3_new(4,k) =  Dk(2,2)
!~             D3_mid(:,k) = (1./2.) * (DD(:,k) + D3_new(:,k))
!~         end do
!~         
!~         ! COMPUTE  k4 = U3_new from X3_new  and D3_new
!~         do k=1,N_part
!~             k1_X(:,k) = U(:,k)
!~             k2_X(:,k) = U1_mid(:,k)
!~             k3_X(:,k) = U2_mid(:,k)
!~             
!~             call diffusion_field_barenblatt(U3_new, X3_new, M, D3_new, m_barenblatt, hx_remap, hy_remap, N_part)
!~             U3_new = -U3_new
!~             k4_X(:,k) = U3_new(:,k)
!~             
!~             k1_D(:,k) = DU(:,k)
!~             k2_D(:,k) = DU1_mid(:,k)
!~             k3_D(:,k) = DU2_mid(:,k)
!~             
!~             X_new(:,k) = Xpart(:,k) + (dt./6.) * (k1_X(:,k) + 2*k2_X(:,k) + 2*k3_X(:,k) + k4_X(:,k))
!~             D_new(:,k) = DD(:,k) + (dt./6.) * (k1_D(:,k) + 2*k2_D(:,k) + 2*k3_D(:,k) + k4_D(:,k))
!~             
!~         end do
!~         
!~         
!~         
!~         
!~         
!~         
!~             a=DD(1,k)
!~             b=DD(2,k)
!~             c=DD(3,k)
!~             d=DD(4,k)
!~             det_Dk = abs(a*d-b*c)
!~             norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk           
!~             if (norm_inf_Dkm1 > Norm_inf_Dm1) then
!~                 Norm_inf_Dm1 = norm_inf_Dkm1
!~                 indice_max_norm_Dm1 = k
!~             end if
!~ 
!~     end subroutine   
    
    !-------------------------------------------------------
    subroutine move_particles(Xold, t, dt, time_scheme, Xnew, N_part)
        integer, intent(in)                             :: N_part
        double precision, dimension (2,N_part), intent(in)        :: Xold        
        double precision, dimension (2, N_part), intent(out)      :: Xnew
        double precision, intent(in)                              :: t, dt        
        integer, intent(in)                             :: time_scheme
        
        double precision, dimension(2,N_part)         :: k1, k2, k3, k4
       
        Xnew=0.
        if (time_scheme == 1) then
            k1=compute_U_part(Xold, t, N_part)
            Xnew(1, :) = Xold(1, :) + dt*k1(1,:)
            Xnew(2, :) = Xold(2, :) + dt*k1(2,:)
        else if (time_scheme == 2) then !Runge Kutta of order 2 
            k1=compute_U_part(Xold, t, N_part)
            k2=compute_U_part(Xold +(dt/2.)*k1, t+(dt/2.), N_part)
            Xnew(1, :) = Xold(1, :) + dt*k2(1,:)
            Xnew(2, :) = Xold(2, :) + dt*k2(2,:)
        else if (time_scheme == 3) then !Runge Kutta of order 4
            k1=compute_U_part(Xold, t, N_part)
            k2=compute_U_part(Xold +(dt/2.)*k1, t+(dt/2.), N_part)
            k3=compute_U_part(Xold +(dt/2.)*k2, t+(dt/2.), N_part)
            k4=compute_U_part(Xold +dt*k3, t+dt, N_part )
            Xnew(1, :) = Xold(1, :) + (dt/6.) * (k1(1, :) + 2.*k2(1, :) + 2.*k3(1, :) + k4(1, :))
            Xnew(2, :) = Xold(2, :) + (dt/6.) * (k1(2, :) + 2.*k2(2, :) + 2.*k3(2, :) + k4(2, :))
    
        else 
            print *, 'no correct time_scheme given in data_code_LTP_2D.py'
        end if        
    end subroutine
    
    !!!segmentation fault for N=1000000 . dont know why.
    subroutine move_particles_test(Xold, t, dt, time_scheme, Xnew, N_part)
        integer, intent(in)                             :: N_part
        double precision, dimension (2,N_part), intent(in)        :: Xold        
        double precision, dimension (2, N_part), intent(out)      :: Xnew
        double precision, intent(in)                              :: t, dt        
        integer, intent(in)                             :: time_scheme
        
        double precision, dimension(2,N_part)         :: k1, k2, k3, k4
        integer(8)              :: k, nb_proc, chunk
        integer, external       :: OMP_GET_NUM_THREADS
       !$OMP PARALLEL 
        nb_proc = OMP_GET_NUM_THREADS()
        chunk = N_part / nb_proc +1        
        !$OMP END PARALLEL 
        
        Xnew=0.
        if (time_scheme == 2) then  
            k1=compute_U_part(Xold, t, N_part)
            k2=compute_U_part(Xold +(dt/2.)*k1, t+(dt/2.), N_part)
            Xnew(1, :) = Xold(1, :) + dt*k2(1,:)
            Xnew(2, :) = Xold(2, :) + dt*k2(2,:)            
        else if (time_scheme == 3) then                           
            !$OMP PARALLEL DO shared(chunk) private(k,k1,k2,k3,k4) schedule(dynamic,chunk)
            do k=1,N_part 
                k1(:,k)=adev(t, Xold(1,k), Xold(2,k))
                k2(:,k)=adev(t+(dt/2.), Xold(1,k)+(dt/2.)*k1(1,k), Xold(2,k)+(dt/2.)*k1(2,k))
                k3(:,k)=adev(t+(dt/2.), Xold(1,k)+(dt/2.)*k2(1,k), Xold(2,k)+(dt/2.)*k2(2,k))
                k4(:,k)=adev(t+dt, Xold(1,k)+dt*k3(1,k), Xold(2,k)+dt*k3(2,k))
                Xnew(1, k) = Xold(1, k) + (dt/6.) * (k1(1, k) + 2.*k2(1, k) + 2.*k3(1, k) + k4(1, k))
                Xnew(2, k) = Xold(2, k) + (dt/6.) * (k1(2, k) + 2.*k2(2, k) + 2.*k3(2, k) + k4(2, k))
            end do 
            !$OMP END PARALLEL DO
        else 
            print *, 'no correct time_scheme given in data_code_LTP_2D.py'
        end if        
    end subroutine

    !-------------------------------------------------------    
    !! gives the number of the particle at position (X0,X1) 
    !! if there is no particle the number is -1
    !! if there is a particle the number is k -1 (for python use must add + 1 for fortran use)
    function number_from_pos(X,X0,X1,N) result (part_pos)
        integer, intent(in)                 :: N
        double precision, dimension(2,N), intent(in)  :: X
        double precision, intent(in)                  :: X0, X1
        integer                             :: part_pos
        
        integer :: k
        double precision  :: eps
        !NOT PRETTY !!
        eps=0.00000001
        part_pos = -1
        do k=1,N
            if ((abs(X(1,k) - X0) <= eps) .AND. (abs(X(2,k) - X1) <=eps)) then
                ! there is a -1 to get indices in python where it starts from 0
                part_pos=k-1
            end if
        end do
    end function number_from_pos

    !-------------------------------------------------------    
    !! creates a matrix 'tab' of 4 rows of length N
    !! it contains the cardinal particles at North (tab(k,1)), South (tab(k,2)), 
    !! East (tab(k,3))and West (tab(k,4))
    !! if there is no particle, at north of particle k for example tab(k,1) = -1
    !! if there is a particle at north of particle k then tab(k,1) = indic_particle_north - 1
    !! the -1 is because number_from_pos translates indices in python arrays (which start from 0) 
    subroutine tab_cardinal_particles(X,hx_remap,hy_remap,tab,N)
        integer, intent(in)                    :: N
        double precision, dimension(2,N), intent(in)     :: X
        double precision, intent(in)                     :: hx_remap, hy_remap
        integer*8, dimension(N,4), intent(out) :: tab
        
        integer :: k 
        
        do k=1,N
            tab(k,1) = number_from_pos(X,X(1,k),X(2,k)+hy_remap,N)
            tab(k,2) = number_from_pos(X,X(1,k),X(2,k)-hy_remap,N)
            tab(k,3) = number_from_pos(X,X(1,k)+hx_remap,X(2,k),N)
            tab(k,4) = number_from_pos(X,X(1,k)-hx_remap,X(2,k),N)
        end do
    end subroutine
    
    !-------------------------------------------------------    
    !! update matrix D (so all matrix Dkn) by finite differences on particles positions
    !! there are various schemes used (centered, forward, backward) considering the particle 
    !! doesn't have any neighboors
    subroutine update_D_finite_differences(Xpart, tab, dt ,hx , hy, indice_max_norm_Dm1, Norm_inf_Dm1, Dout, N_part)
        integer, intent(in)                          :: N_part
        double precision, dimension(2,N_part), intent(in)      :: Xpart
        integer*8, dimension(N_part,4), intent(in)   :: tab
        double precision, intent(in)                           :: hx, hy, dt                
        double precision, dimension(4,N_part), intent(out)     :: Dout
        double precision, intent(out)                          :: Norm_inf_Dm1
        integer(8), intent(out)                      :: indice_max_norm_Dm1

        integer*8               :: k, ind_N, ind_S, ind_E, ind_W
        double precision                  :: det_Jk, det_Dk, a, b, c, d, norm_inf_Dkm1
        double precision, dimension(2,2)  :: Jk, Dk
                

        !print *," tab : " , tab
        do k=1,N_part
            !----- calcul de Jkn -----
            Jk=0.        
            Dk=0.
            ind_N=tab(k,1) !+ 1
            ind_S=tab(k,2) !+ 1
            ind_E=tab(k,3) !+ 1
            ind_W=tab(k,4) !+ 1
            !print *,"k ,  ind_N , ind_S , ind_E , ind_W : " , k  , ind_N , ind_S , ind_E , ind_W

            !!! full West (one particle on East)
            if (((ind_W) == -1)  .AND. ((ind_N) == -1) .AND. ((ind_S) == -1) .AND. ((ind_E) /= -1)) then
                Jk(1,1)=(Xpart(1,ind_E) - Xpart(1,k))/hx
                Jk(2,1)=(Xpart(2,ind_E) - Xpart(2,k))/hx
                !!! Jk(1,2)=(Xpart(1,ind_N) - Xpart(1,ind_S))/hy
                !!! Jk(2,2)=(Xpart(2,ind_N) - Xpart(2,ind_S))/hy    
            !!! full North (one particle on South)
            else if (((ind_W) == -1)  .AND. ((ind_N) == -1) .AND. ((ind_E) == -1) .AND. ((ind_S) /= -1)) then
                !!! Jk(1,1)=(Xpart(1,ind_E) - Xpart(1,k))/hx
                !!! Jk(2,1)=(Xpart(2,ind_E) - Xpart(2,k))/hx
                Jk(1,2)=(Xpart(1,k) - Xpart(1,ind_S))/hy
                Jk(2,2)=(Xpart(2,k) - Xpart(2,ind_S))/hy
            !!! full East (one particle on West)
            else if (((ind_E) == -1)  .AND. ((ind_N) == -1) .AND. ((ind_S) == -1) .AND. ((ind_W) /= -1)) then
                Jk(1,1)=(Xpart(1,k) - Xpart(1,ind_W))/hx
                Jk(2,1)=(Xpart(2,k) - Xpart(2,ind_W))/hx
                !!! Jk(1,2)=(Xpart(1,ind_N) - Xpart(1,ind_S))/hy
                !!! Jk(2,2)=(Xpart(2,ind_N) - Xpart(2,ind_S))/hy
            !!! full South (one particle on North)
            else if (((ind_W) == -1)  .AND. ((ind_E) == -1) .AND. ((ind_S) == -1) .AND. ((ind_N) /= -1)) then
                !!! Jk(1,1)=(Xpart(1,ind_E) - Xpart(1,k))/hx
                !!! Jk(2,1)=(Xpart(2,ind_E) - Xpart(2,k))/hx
                Jk(1,2)=(Xpart(1,ind_N) - Xpart(1,k))/hy
                Jk(2,2)=(Xpart(2,ind_N) - Xpart(2,k))/hy
            !!! between North and South        
            else if (((ind_W) == -1) .AND. ((ind_E) == -1) .AND. ((ind_N) /= -1) .AND. ((ind_S) /= -1)) then
                !!! derivative in x ??
                !!! Jk(1,1)=(Xpart(1,ind_E) - Xpart(1,k))/hx
                !!! Jk(2,1)=(Xpart(2,ind_E) - Xpart(2,k))/hx
                Jk(1,2)=(Xpart(1,ind_N) - Xpart(1,ind_S))/(2*hy)
                Jk(2,2)=(Xpart(2,ind_N) - Xpart(2,ind_S))/(2*hy)
            !!! between West and East
            else if (((ind_N) == -1) .AND. ((ind_S) == -1) .AND. ((ind_W) /= -1) .AND. ((ind_E) /= -1)) then
                !!! derivative in y ??
                Jk(1,1)=(Xpart(1,ind_E) - Xpart(1,ind_W))/(2*hx)
                Jk(2,1)=(Xpart(2,ind_E) - Xpart(2,ind_W))/(2*hx)
                !!! Jk(1,2)=(Xpart(1,ind_N) - Xpart(1,ind_S))/hy
                !!! Jk(1,1)=(Xpart(1,ind_N) - Xpart(1,ind_S))/hy
            !!! corner NW
            else if ((ind_N == -1) .AND. (ind_W == -1) .AND. (ind_E /= -1) .AND. (ind_S /= -1)) then
                Jk(1,1)=(Xpart(1,ind_E+1) - Xpart(1,k))/hx
                Jk(2,1)=(Xpart(2,ind_E+1) - Xpart(2,k))/hx
                Jk(1,2)=(Xpart(1,k) - Xpart(1,ind_S+1))/hy
                Jk(2,2)=(Xpart(2,k) - Xpart(2,ind_S+1))/hy
            !!! corner NE
            else if ((ind_N == -1) .AND. (ind_E == -1) .AND. (ind_W /= -1) .AND. (ind_S /= -1)) then
                Jk(1,1)=(Xpart(1,k) - Xpart(1,ind_W+1))/hx
                Jk(2,1)=(Xpart(2,k) - Xpart(2,ind_W+1))/hx
                Jk(1,2)=(Xpart(1,k) - Xpart(1,ind_S+1))/hy
                Jk(2,2)=(Xpart(2,k) - Xpart(2,ind_S+1))/hy
            !!! corner SE
            else if ((ind_S == -1) .AND. (ind_E == -1) .AND. (ind_W /= -1) .AND. (ind_N /= -1)) then
                Jk(1,1)=(Xpart(1,k) - Xpart(1,ind_W+1))/hx
                Jk(2,1)=(Xpart(2,k) - Xpart(2,ind_W+1))/hx
                Jk(1,2)=(Xpart(1,ind_N+1) - Xpart(1,k))/hy
                Jk(2,2)=(Xpart(2,ind_N+1) - Xpart(2,k))/hy
            !!! corner SW
            else if ((ind_S == -1) .AND. (ind_W == -1) .AND. (ind_N /= -1) .AND. (ind_E /= -1)) then
                Jk(1,1)=(Xpart(1,ind_E+1) - Xpart(1,k))/hx
                Jk(2,1)=(Xpart(2,ind_E+1) - Xpart(2,k))/hx
                Jk(1,2)=(Xpart(1,ind_N+1) - Xpart(1,k))/hy
                Jk(2,2)=(Xpart(2,ind_N+1) - Xpart(2,k))/hy
            !!! north but not in corner
            else if ((ind_N == -1) .AND. (ind_W /= -1) .AND. (ind_E /= -1) .AND. (ind_S /= -1)) then
                Jk(1,1)=(Xpart(1,ind_E+1) - Xpart(1,ind_W+1))/(2*hx)
                Jk(2,1)=(Xpart(2,ind_E+1) - Xpart(2,ind_W+1))/(2*hx)
                Jk(1,2)=(Xpart(1,k) - Xpart(1,ind_S+1))/hy
                Jk(2,2)=(Xpart(2,k) - Xpart(2,ind_S+1))/hy
            !!! east but not in corner
            else if ((ind_E == -1) .AND. (ind_N /= -1) .AND. (ind_S /= -1) .AND. (ind_W /= -1)) then
                Jk(1,1)=(Xpart(1,k) - Xpart(1,ind_W+1))/hx
                Jk(2,1)=(Xpart(2,k) - Xpart(2,ind_W+1))/hx
                Jk(1,2)=(Xpart(1,ind_N+1) - Xpart(1,ind_S+1))/(2*hy)
                Jk(2,2)=(Xpart(2,ind_N+1) - Xpart(2,ind_S+1))/(2*hy)
            !!! south but not in corner
            else if ((ind_S == -1) .AND. (ind_W /= -1) .AND. (ind_E /= -1) .AND. (ind_N /= -1)) then
                Jk(1,1)=(Xpart(1,ind_E+1) - Xpart(1,ind_W+1))/(2*hx)
                Jk(2,1)=(Xpart(2,ind_E+1) - Xpart(2,ind_W+1))/(2*hx)
                Jk(1,2)=(Xpart(1,ind_N+1) - Xpart(1,k))/hy
                Jk(2,2)=(Xpart(2,ind_N+1) - Xpart(2,k))/hy
            !!! west but not in corner
            else if ((ind_W == -1) .AND. (ind_N /= -1) .AND. (ind_S /= -1) .AND. (ind_E /= -1)) then
                Jk(1,1)=(Xpart(1,ind_E+1) - Xpart(1,k))/hx
                Jk(2,1)=(Xpart(2,ind_E+1) - Xpart(2,k))/hx
                Jk(1,2)=(Xpart(1,ind_N+1) - Xpart(1,ind_S+1))/(2*hy)
                Jk(2,2)=(Xpart(2,ind_N+1) - Xpart(2,ind_S+1))/(2*hy)
            !!! not in "boundary"
            else 
                Jk(1,1)=(Xpart(1,ind_E+1) - Xpart(1,ind_W+1))/(2*hx)
                Jk(2,1)=(Xpart(2,ind_E+1) - Xpart(2,ind_W+1))/(2*hx)
                Jk(1,2)=(Xpart(1,ind_N+1) - Xpart(1,ind_S+1))/(2*hy)
                Jk(2,2)=(Xpart(2,ind_N+1) - Xpart(2,ind_S+1))/(2*hy)
            end if
            !----- inverse of Jk -----            
            Dk(1,1) = Jk(2,2)
            Dk(2,1) = -Jk(2,1)
            Dk(1,2) = -Jk(1,2)
            Dk(2,2) = Jk(1,1)
            det_Jk = Jk(1,1) * Jk(2,2) - Jk(2,1) * Jk(1,2)
            Dk = (1/det_Jk) * Dk            
            a = Dk(1,1)
            b = Dk(1,2)
            c = Dk(2,1)
            d = Dk(2,2)
            det_Dk = abs(a*d-b*c)
            norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
            !print *,'k , norm_inf_Dkm1 , Norm_inf_Dm1 = ' , k, norm_inf_Dkm1, Norm_inf_Dm1            
            if (norm_inf_Dkm1 > Norm_inf_Dm1) then
                Norm_inf_Dm1 = norm_inf_Dkm1
                indice_max_norm_Dm1 = k
            end if
            !-----  update D -----
            Dout(1,k) = Dk(1,1)
            Dout(2,k) = Dk(1,2)
            Dout(3,k) = Dk(2,1)
            Dout(4,k) = Dk(2,2)            
        end do
    end subroutine

    !-------------------------------------------------------    
    !! compute trapezoidal quadrature points for integral computation
    !! on convolution product
    !! We assume that the Xgrid_convol and Ygrid_convol are well ordered grid 
    !! in the sense that : Xgrid_convol(1) == MINVAL(Xgrid_convol)
    !!                     Xgrid_convol(nqx) == MAXVAL(Xgrid_convol)
    !!                     Ygrid_convol(1) == MINVAL(Ygrid_convol)
    !!                     Ygrid_convol(nqy) == MINVAL(Ygrid_convol)    
    subroutine quadrature_pts_trapezoidal_method(quadrature_pts_convol, Xgrid_convol, Ygrid_convol, nqx, nqy)
        integer, intent(in)                                  :: nqx, nqy
        double precision, dimension (nqx), intent(in)        :: Xgrid_convol
        double precision, dimension (nqy), intent(in)        :: Ygrid_convol
        double precision, dimension (nqx,nqy), intent(out)   :: quadrature_pts_convol
        
        integer(8)  :: i, j
                
        do i=1,nqx
            do j=1,nqy
                !on borders set quadture points to 2 (corner will be changed right after)
                if ((i == 1) .OR. (j==1) .OR. (i==nqx) .OR. (j==nqy)) then
                    quadrature_pts_convol(i,j) = 2.
                ! everywhere else set quadrature points to 4
                else                
                    quadrature_pts_convol(i,j) = 4.
                end if
            end do
        end do
        
        ! on corner set quadrature points to 1
        quadrature_pts_convol(1,1) = 1.
        quadrature_pts_convol(1,nqy) = 1.
        quadrature_pts_convol(nqx,1) = 1.
        quadrature_pts_convol(nqx,nqy) = 1.
    end subroutine
    
end module calculsfor_var_modf90
