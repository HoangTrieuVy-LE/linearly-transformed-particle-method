#  !/usr/bin/env python
# - * - codin: utf-8 -*-
# derniere modif le 09-10-15 par Frederique
from __future__ import division
import numpy  #as np
from scipy import *  #as sp

import config

#~ from shapefunction2D import *
#~ import shapefunction2D

import barenblatt_fcts
import diffusion_fcts

from _calculsfor_f90 import  calculsfor_ini_modf90
from _calculsfor_f90 import  calculsfor_rec_modf90
from _calculsfor_f90 import  calculsfor_var_modf90

def compute_solution(Xg,Yg, t) :
    mx = size(Xg)
    my = size(Yg)
    sol = numpy.zeros((mx,my))       
    for i in range(mx) :        
        for j in range(my) :
            if (config.name_solution == 'barenblatt') :
                sol[i,j] = barenblatt_fcts.barenblatt(Xg[i],Yg[j],t)            
            elif (config.name_solution == 'diffusion') :
                #~ sol[i,j] = diffusion_fcts.diffusion_sol(Xg[i],Yg[j],t)
                sol[i,j] = diffusion_fcts.diffusion_sol_analytic(Xg[i],Yg[j],t)
    return sol

def compute_solution_heat_kernel(Xg,Yg, t) :
    mx = size(Xg)
    my = size(Yg)
    sol = numpy.zeros((mx,my))       
    for i in range(mx) :        
        print (i)
        for j in range(my) :
            if (config.name_solution == 'barenblatt') :
                sol[i,j] = barenblatt_fcts.barenblatt(Xg[i],Yg[j],t)            
            elif (config.name_solution == 'diffusion') :
                sol[i,j] = diffusion_fcts.diffusion_sol(Xg[i],Yg[j],t)
                #~ sol[i,j] = diffusion_fcts.diffusion_sol_analytic(Xg[i],Yg[j],t)
    return sol

def diffusion_field_barenblatt(X, N , M , D,m_barenblatt,hx_remap,hy_remap):
    diffusion_barenblatt_f2py  = calculsfor_rec_modf90.diffusion_field_barenblatt(X, M , D, m_barenblatt, hx_remap, hy_remap, N)
    return diffusion_barenblatt_f2py

def error(analytic_sol, approx_sol, error_type='L_infinity'): 
    size=analytic_sol.shape[0] * analytic_sol.shape[1]
    vect_analytic_sol = analytic_sol.reshape(size,1)
    vect_approx_sol = approx_sol.reshape(size,1)
    if (error_type == 'L_infinity') :
        error_sol = max(abs(vect_analytic_sol - vect_approx_sol)) / max(abs(vect_analytic_sol))
    elif (error_type == 'absolute_error') :
        error_sol = max(abs(vect_analytic_sol - vect_approx_sol))			    
    else :
        raise ValueError("incorrect type for error_type = "+str(error_type))
    return error_sol

def error_Lp(Xg, Yg, analytic_sol, approx_sol, norm_Lp_order): 
    Z=abs(approx_sol - analytic_sol)**norm_Lp_order
    N1 = analytic_sol.shape[0] 
    N2 = analytic_sol.shape[1]

    a1 = min(Xg) 
    b1 = max(Xg) 
    a2 = min(Yg)
    b2 = max(Yg)
    hgx = (b1-a1)/N1;
    hgy = (b2-a2)/N2;
    
    ## this is the 2D trapezoidal integration rule based on the step of grid reconstruction
    Lp_Z_int = hgx*hgy*(sum(Z[:,:]) - 0.5*(sum(Z[1,:]+Z[N1-1,:]) + sum(Z[:,1]+Z[:,N2-1])) + 0.25*(Z[1,1]+Z[1,N2-1]+Z[N1-1,1]+Z[N1-1,N2-1]))
    Lp_Z_int = (Lp_Z_int)**(1./norm_Lp_order)
    
    Lp_analytic_sol = hgx*hgy*(sum(analytic_sol[:,:]) - 0.5*(sum(analytic_sol[1,:]+analytic_sol[N1-1,:]) + sum(analytic_sol[:,1]+analytic_sol[:,N2-1])) + 0.25*(analytic_sol[1,1]+analytic_sol[1,N2-1]+analytic_sol[N1-1,1]+analytic_sol[N1-1,N2-1]))
    Lp_analytic_sol = Lp_analytic_sol**(1./norm_Lp_order)
    error = Lp_Z_int / Lp_analytic_sol
    return error
    
#~ f = @(x,y) sqrt((exp(sin(x.*pi - y.*pi)))+(exp(cos(x.*y.*pi))));
#~ a1 = 0; b1 = 3; N1 = 20;
#~ h1 = (b1-a1)/N1;
#~ a2 = -1; b2 = 4; N2 = 30;
#~ h2 = (b2-a2)/N2;
#~ [X,Y] = meshgrid(a1:h1:b1, a2:h2:b2);
#~ Z = f(X,Y); 
#~ T = h1*h2*(sum(Z(:)) - 0.5*(sum(Z(1,:)+Z(end,:)) + sum(Z(:,1)+Z(:,end))) + 0.25*(Z(1,1)+Z(1,end)+Z(end,1)+Z(end,end)))
#~ 

def MaJ_D_diffusion(D, X, M, t, dt,hx_remap, hy_remap):
	N = len(X[0,:])
	if N==0:
		print( 'Plus de particules !!!')
    # parameter 1 in 'update_d_schemes' is the number for euler explicit time scheme	
    #diff_field_jac = calculsfor_rec_modf90.diffusion_field_jacobian_barenblatt(X, M, D, m_barenblatt, hx_remap, hy_remap, N)
	update_d_scheme = 0	
	if config.D_method_scheme  == 'euler_explicit' :
		update_d_scheme = 1
		indice_max_norm_Dm1, Norm_inf_Dm1, Df2py  = calculsfor_var_modf90.update_d_diffusion(X, D, M, t, dt, update_d_scheme, hx_remap, hy_remap, N)    	
		indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)
	else :
		raise ValueError('wrong D_method_scheme value given : '+str(D_method_scheme))

	return Df2py, Norm_inf_Dm1, indice_max_norm_Dm1    



def MaJ_D_barenblatt(D, X, M, t, dt,hx_remap, hy_remap):
	N = len(X[0,:])
	if N==0:
		print ('plus de particules !!!')
    # parameter 1 in 'update_d_schemes' is the number for euler explicit time scheme	
    #diff_field_jac = calculsfor_rec_modf90.diffusion_field_jacobian_barenblatt(X, M, D, m_barenblatt, hx_remap, hy_remap, N)
	update_d_scheme = 0	
	if config.D_method_scheme == 'euler_implicit' :
		update_d_scheme = 3
		indice_max_norm_Dm1, Norm_inf_Dm1, Df2py  = calculsfor_var_modf90.update_d_schemes_diffusion_barenblatt(X, D, M, t, dt, update_d_scheme, config.m_barenblatt, hx_remap, hy_remap, N)    	
		indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)
		print (indice_max_norm_Dm1, Norm_inf_Dm1, Df2py)
	elif config.D_method_scheme  == 'euler_explicit' :
		update_d_scheme = 1
		indice_max_norm_Dm1, Norm_inf_Dm1, Df2py  = calculsfor_var_modf90.update_d_schemes_diffusion_barenblatt(X, D, M, t, dt, update_d_scheme, config.m_barenblatt, hx_remap, hy_remap, N)    	
		indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)
	elif config.D_method_scheme  == 'RK2' :
		update_d_scheme = 2
		indice_max_norm_Dm1, Norm_inf_Dm1, Df2py  = calculsfor_var_modf90.update_d_schemes_diffusion_barenblatt(X, D, M, t, dt, update_d_scheme, config.m_barenblatt, hx_remap, hy_remap, N)    	
		indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)	
	else :
		raise ValueError('wrong D_method_scheme value given : '+str(D_method_scheme))

	return Df2py, Norm_inf_Dm1, indice_max_norm_Dm1    

def MaJ_D_barenblatt_mid_scheme(D, X, X_middle, M, t, dt, hx_remap, hy_remap):
	N = len(X[0,:])
	if N==0:
		print( 'plus de particules !!!'    )
	update_d_scheme = 0
	if config.time_scheme == 'middle_point' :		
		D_middle, indice_max_norm_Dm1, Norm_inf_Dm1, D_new = calculsfor_var_modf90.update_d_barenblatt_mid_scheme(X, X_middle, D, M, t, config.dt,1, config.m_barenblatt, hx_remap, hy_remap, N)
		indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)
	else :
		raise ValueError('wrong time_scheme value given : '+str(time_scheme))	
		
	return D_middle, D_new, Norm_inf_Dm1, indice_max_norm_Dm1    


def MaJ_D_convolution_barenblatt(D, X, M, t, dt,hx_remap, hy_remap, epsilon, Xgrid_convol, Ygrid_convol, quadrature_pts_convol):
    N = len(X[0,:])
    if N==0:
        print ('plus de particules !!!')
    if (config.D_method == 'D_convol') : # an other check        
        indice_max_norm_Dm1, Norm_inf_Dm1, D_new  = calculsfor_var_modf90.update_d_convolution_method_diffusion_barenblatt(X, D, M, t, dt, config.m_barenblatt, epsilon, hx_remap, hy_remap, Xgrid_convol, Ygrid_convol, quadrature_pts_convol, len(Xgrid_convol), len(Ygrid_convol), N)
        indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)
    else :
        raise ValueError('wrong D_method value given : '+str(D_method))

    return D_new , Norm_inf_Dm1, indice_max_norm_Dm1



def MaJ_D_hybrid_barenblatt(D, X, M, t, dt,hx_remap, hy_remap):
    N = len(X[0,:])
    if N==0:
        print ('plus de particules !!!')
    update_d_scheme = 0
    if config.D_method == "hybrid" :
        print ("launch :  update_d_hybrid_diffusion_barenblatt")
        print ("epsilon = " , epsilon)
        indice_max_norm_Dm1, Norm_inf_Dm1, Df2py  = calculsfor_var_modf90.update_d_hybrid_diffusion_barenblatt(X, D, M, t, dt, m_barenblatt, epsilon, hx_remap, hy_remap, N)
    else :
        raise ValueError('wrong D_method_scheme value given : '+str(config.D_method_scheme))

    return Df2py, Norm_inf_Dm1, indice_max_norm_Dm1   



def MaJ_D_FD_scheme_avec_f2py(X, tab, dt ,hx , hy) :
	Norm_inf_Dm1 = 0.
	N = len(X[0,:])	
	if N==0:
		print ('plus de particules !!!')
	indice_max_norm_Dm1, Norm_inf_Dm1, Df2py = calculsfor_var_modf90.update_d_finite_differences(X, tab, dt ,hx , hy, N)
	indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)            
	return Df2py, Norm_inf_Dm1, indice_max_norm_Dm1		

def MaJ_D_FD_scheme(X, D , tab, t, dt ,hx , hy):
	Dfor, Norm_inf_Dm1for, indice_max_norm_Dm1for = MaJ_D_FD_scheme_avec_f2py(X, tab, dt ,hx , hy)
	return Dfor, Norm_inf_Dm1for, indice_max_norm_Dm1for


def MAJ_grille(X,Ix,Iy):
    modif_grille = 0
    xmin = X[0,:].min() 
    if xmin<Ix[0]:
        Ix[0] = xmin
        modif_grille = 1
    xmax = X[0,:].max()
    if xmax>Ix[1]:
        Ix[1] = xmax
        modif_grille = 1        
    ymin = X[1,:].min()
    if ymin<Iy[0]:
        Iy[0] = ymin
        modif_grille = 1
    ymax = X[1,:].max()
    if ymax>Iy[1]:
        Iy[1] = ymax
        modif_grille = 1 
    return  Ix, Iy, modif_grille



# ===== Grille =====
def make_grid_unif(Ix, Iy, m):    
    xmin = Ix[0]
    xmax = Ix[1] #min(X[0,:]) + 3 * epsilon
    ymin = Iy[0]  #min(X[0,:]) - 3 * epsilon
    ymax = Iy[1] #min(X[0,:]) + 3 * epsilon
    dx = abs(xmax - xmin) / (m - 1)
    dy = abs(ymax - ymin) / (m - 1)
    delta = min(dx, dy)
    mx = abs(xmax - xmin) / delta + 1
    my = abs(ymax - ymin) / delta + 1
    Xgrid = numpy.linspace(xmin, xmax, mx)  # numpy.array([xmin + i * delta for i in range(mx)])
    Ygrid = numpy.linspace(ymin, ymax, my)  # numpy.array([ymin + j * delta for j in range(my)])
    return [Xgrid, Ygrid]
    

def move_particles_diffusion_barenblatt(X, N , M , D,m_barenblatt,hx_remap,hy_remap) :    
	if time_scheme == 'euler_explicit' :
		k1 = - diffusion_field_barenblatt(X, N , M , D,m_barenblatt,hx_remap,hy_remap)  
		X[0, :] += dt*k1[0,:]
		X[1, :] += dt*k1[1,:]	
	#~ elif time_scheme == 'RK2' : 
		#~ k1=-diffusion_field_barenblatt(X, N , M , D,m_barenblatt,hx_remap,hy_remap)  		
		#~ k2=-diffusion_field_barenblatt(X+(dt/2.)*k1, N , M , D,m_barenblatt,hx_remap,hy_remap)
		#~ X[0, :] += dt*k2[0,:]
		#~ X[1, :] += dt*k2[1,:]        	
	#~ elif time_scheme == 'RK4' : 
		#~ k1 = - diffusion_field_barenblatt(X, N , M , D,m_barenblatt,hx_remap,hy_remap)  		
		#~ k2 = - diffusion_field_barenblatt(X + (dt/2.)*k1, N , M , D,m_barenblatt,hx_remap,hy_remap)
		#~ k3 = - diffusion_field_barenblatt(X + (dt/2.)*k2, N , M , D,m_barenblatt,hx_remap,hy_remap)
		#~ k4 = - diffusion_field_barenblatt(X + dt*k3, N , M , D,m_barenblatt,hx_remap,hy_remap)
		#~ X[0, :] += (dt/6.) * (k1[0,:] + 2*k2[0,:] + 2*k3[0,:] + k4[0,:])
		#~ X[1, :] += (dt/6.) * (k1[1,:] + 2*k2[1,:] + 2*k3[1,:] + k4[1,:])
	else : 
		raise ValueError('uncorrect time scheme given')
	return X
	
def update_d_hessian_analytic(X, N, D, m_barenblatt, t, dt, hx_remap, hy_remap) :
    hess_rho_analytic = numpy.zeros([2,2])
    Dk = numpy.zeros([2,2])
    Dnew=numpy.zeros([4,N])
    Norm_inf_Dm1 = 0.
    indice_max_norm_Dm1=0
    for k in range(0,N) :
        hess_rho_analytic = barenblatt_fcts.hess_barenblatt(X[0,k],X[1,k],t)
        invJk = calculsfor_var_modf90.compute_expo_mat_order(-dt*(-m_barenblatt*hess_rho_analytic), 30)
        Dk[0,0]=D[0,k]
        Dk[0,1]=D[1,k]
        Dk[1,0]=D[2,k]
        Dk[1,1]=D[3,k]
        Dk = numpy.dot(Dk,invJk)
        a=D[0,k]
        b=D[1,k]
        c=D[2,k]
        d=D[3,k]
        det_Dk = abs(a*d-b*c)
        norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
        if (norm_inf_Dkm1 > Norm_inf_Dm1) :
            Norm_inf_Dm1 = norm_inf_Dkm1
            indice_max_norm_Dm1 = k
        Dnew[0,k] = Dk[0,0]
        Dnew[1,k] = Dk[0,1]
        Dnew[2,k] = Dk[1,0]
        Dnew[3,k] = Dk[1,1]
    return Dnew, Norm_inf_Dm1, indice_max_norm_Dm1
    
def update_d_hessian_analytic_middle_point(X, X_mid, N, D, m_barenblatt, t, dt, hx_remap, hy_remap) :
    hess_rho_analytic = numpy.zeros([2,2])
    Dk = numpy.zeros([2,2])
    Dk_mid = numpy.zeros([2,2])
    Jk_mid = numpy.zeros([2,2])
    
    D_new = numpy.zeros([4,N])
    D_mid = numpy.zeros([4,N])
    J_mid = numpy.zeros([4,N])
    
    for k in range(0,N) :
        hess_rho_analytic = barenblatt_fcts.hess_barenblatt(X[0,k],X[1,k],t)
        invJk = calculsfor_var_modf90.compute_expo_mat_order(-dt/2.*(-m_barenblatt * hess_rho_analytic ), 10)
        Jk_mid = calculsfor_var_modf90.compute_expo_mat_order(dt/2.*(-m_barenblatt * hess_rho_analytic ), 10)
        J_mid[0,k] = Jk_mid[0,0]
        J_mid[1,k] = Jk_mid[0,1]
        J_mid[2,k] = Jk_mid[1,0]
        J_mid[3,k] = Jk_mid[1,1]
        Dk[0,0]=D[0,k]
        Dk[0,1]=D[1,k]
        Dk[1,0]=D[2,k]
        Dk[1,1]=D[3,k]
        Dk = numpy.dot(Dk,invJk)
        D_mid[0,k] = Dk[0,0]
        D_mid[1,k] = Dk[0,1]
        D_mid[2,k] = Dk[1,0]
        D_mid[3,k] = Dk[1,1]
        
    Norm_inf_Dm1 = 0.
    indice_max_norm_Dm1=0    
    for k in range(0,N) :       
        Jk_mid[0,0] = J_mid[0,k]
        Jk_mid[0,1] = J_mid[1,k]
        Jk_mid[1,0] = J_mid[2,k]
        Jk_mid[1,1] = J_mid[3,k]
        hess_rho_analytic = barenblatt_fcts.hess_barenblatt(X_mid[0,k],X_mid[1,k],t+dt/2.)
        #~ invJk = calculsfor_var_modf90.compute_expo_mat_order(-dt*(-m_barenblatt*numpy.dot(hess_rho_analytic, Jk_mid) ), 10)
        invJk = calculsfor_var_modf90.compute_expo_mat_order(-dt/2.*(-m_barenblatt*hess_rho_analytic), 10)
        Dk_mid[0,0]=D_mid[0,k]
        Dk_mid[0,1]=D_mid[1,k]
        Dk_mid[1,0]=D_mid[2,k]
        Dk_mid[1,1]=D_mid[3,k]    
        Dk = numpy.dot(Dk_mid,invJk)
        a=D[0,k]
        b=D[1,k]
        c=D[2,k]
        d=D[3,k]
        det_Dk = abs(a*d-b*c)
        norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
        if (norm_inf_Dkm1 > Norm_inf_Dm1) :
            Norm_inf_Dm1 = norm_inf_Dkm1
            indice_max_norm_Dm1 = k            
        D_new[0,k] = Dk[0,0]
        D_new[1,k] = Dk[0,1]
        D_new[2,k] = Dk[1,0]
        D_new[3,k] = Dk[1,1]
        
    return D_mid, D_new, Norm_inf_Dm1, indice_max_norm_Dm1
    
def update_d_and_x_middle_point_gradient_and_hessian_analytic(X, M, N, D, m_barenblatt, t, dt, hx_remap, hy_remap) :
    X_old = numpy.zeros([2,N])
    X_mid = numpy.zeros([2,N])
    X_new = numpy.zeros([2,N])
    U     = numpy.zeros([2,N])
    U_mid = numpy.zeros([2,N])
    U_new = numpy.zeros([2,N])
    D_old = numpy.zeros([4,N])
    D_mid = numpy.zeros([4,N])
    D_new = numpy.zeros([4,N])
    
    X_old = X.copy()
    D_old = D.copy()
    for k in range(0,N) :
        U[:,k] = -config.m_barenblatt * barenblatt_fcts.grad_barenblatt(X_old[0,k], X_old[1,k], t)
    X_mid[0,:] = X_old[0,:] + (dt/2.) * U[0,:]
    X_mid[1,:] = X_old[1,:] + (dt/2.) * U[1,:]
    
    D_mid, Norm_inf_Dm1_mid, indice_max_norm_Dm1_mid = update_d_hessian_analytic(X_old, N, D_old, m_barenblatt, t, dt/2., hx_remap, hy_remap)
    
    for k in range(0,N) :
        U_mid[:,k] = -config.m_barenblatt * barenblatt_fcts.grad_barenblatt(X_mid[0,k], X_mid[1,k], t+dt/2.)
    X_new[0,:] = X_old[0,:] + dt * U_mid[0,:]
    X_new[1,:] = X_old[1,:] + dt * U_mid[1,:]
    
    D_new, Norm_inf_Dm1, indice_max_norm_Dm1 = update_d_hessian_analytic(X_mid, N, D_old, m_barenblatt, t+dt/2., dt, hx_remap, hy_remap)
    
    return X_new, D_new, Norm_inf_Dm1, indice_max_norm_Dm1

    
def update_d_and_x_middle_point_gradient_and_hessian_analytic_v2(X, M, N, D, m_barenblatt, t, dt, hx_remap, hy_remap) :
    X_old = numpy.zeros([2,N])
    X_mid = numpy.zeros([2,N])
    X_new = numpy.zeros([2,N])
    U     = numpy.zeros([2,N])
    U_mid = numpy.zeros([2,N])
    U_new = numpy.zeros([2,N])
    D_old = numpy.zeros([4,N])
    D_mid = numpy.zeros([4,N])
    D_new = numpy.zeros([4,N])
    
    X_old = X.copy()
    D_old = D.copy()
    for k in range(0,N) :
        U[:,k] = -config.m_barenblatt * barenblatt_fcts.grad_barenblatt(X_old[0,k], X_old[1,k], t)
    X_mid[0,:] = X_old[0,:] + (dt/2.) * U[0,:]
    X_mid[1,:] = X_old[1,:] + (dt/2.) * U[1,:]
    
    D_mid, D_new, Norm_inf_Dm1, indice_max_norm_Dm1 = update_d_hessian_analytic_middle_point(X_old, X_mid, N, D_old, m_barenblatt, t, dt, hx_remap, hy_remap)
    
    for k in range(0,N) :
        U_mid[:,k] = -config.m_barenblatt * barenblatt_fcts.grad_barenblatt(X_mid[0,k], X_mid[1,k], t+dt/2.)
    X_new[0,:] = X_old[0,:] + dt * U_mid[0,:]
    X_new[1,:] = X_old[1,:] + dt * U_mid[1,:]        
    
    return X_new, D_new, Norm_inf_Dm1, indice_max_norm_Dm1
    
    
def update_d_and_x_middle_point_hessian_analytic(X, M, N, D, m_barenblatt, t, dt, hx_remap, hy_remap) :
    X_old = numpy.zeros([2,N])
    X_mid = numpy.zeros([2,N])
    X_new = numpy.zeros([2,N])
    U     = numpy.zeros([2,N])
    U_mid = numpy.zeros([2,N])
    U_new = numpy.zeros([2,N])
    D_old = numpy.zeros([4,N])
    D_mid = numpy.zeros([4,N])
    D_new = numpy.zeros([4,N])
    
    X_old = X.copy()
    D_old = D.copy()
    
    U = - calculsfor_rec_modf90.diffusion_field_barenblatt(X_old, M , D_old, m_barenblatt, hx_remap, hy_remap, N)
    X_mid[0,:] = X_old[0,:] + (dt/2.) * U[0,:]
    X_mid[1,:] = X_old[1,:] + (dt/2.) * U[1,:]
    
    D_mid, Norm_inf_Dm1_mid, indice_max_norm_Dm1_mid = update_d_hessian_analytic(X_old, N, D_old, m_barenblatt, t, dt/2., hx_remap, hy_remap)
    
    U_mid = - calculsfor_rec_modf90.diffusion_field_barenblatt(X_mid, M , D_mid, m_barenblatt, hx_remap, hy_remap, N)    
    X_new[0,:] = X_old[0,:] + dt * U_mid[0,:]
    X_new[1,:] = X_old[1,:] + dt * U_mid[1,:]
    
    D_new, Norm_inf_Dm1, indice_max_norm_Dm1 = update_d_hessian_analytic(X_mid, N, D_old, m_barenblatt, t+dt/2., dt, hx_remap, hy_remap)
    
    return X_new, D_new, Norm_inf_Dm1, indice_max_norm_Dm1

# in : numpy.array([2,N])
# should also work with a numpy.array([4,N])
def norm_max_all_particles(X) :
	N = len(X[0,:])
	list_max_coordinates = []
	for k in range(0,N) :
		list_max_coordinates.append(max (abs(X[:,k])))	
	return max(list_max_coordinates)
    
# in : numpy.array([2,N])
# should also work with a numpy.array([4,N])
def norm_max_all_matrices(D) :
    N = len(D[0,:])
    list_max = []
    Dk = numpy.zeros([2,2])
    for k in range(0,N) :
        Dk[0,0] = D[0,k]
        Dk[0,1] = D[1,k]
        Dk[1,0] = D[2,k]
        Dk[1,1] = D[3,k]
        
        list_max.append(max(sum(abs(Dk[0,:])) , sum(abs(Dk[1,:])) ))        
    
    return max(list_max)

def relative_error_norm_max_all_particles(X_exact, X):
	return norm_max_all_particles(X-X_exact) / norm_max_all_particles(X_exact)

def absolute_error_norm_max_all_particles(X_exact, X):
	return norm_max_all_particles(X-X_exact)
    
def relative_error_norm_max_all_matrices(D_exact, D):
	return norm_max_all_matrices(D-D_exact) / norm_max_all_matrices(D_exact)

def absolute_error_norm_max_all_matrices(D_exact, D):
	return norm_max_all_matrices(D-D_exact)

def get_quadrature_pts_convol() :
    xg_convol_min = config.Sx[0]-abs(config.Sx[1]-config.Sx[0])/2.
    xg_convol_max = config.Sx[1]+abs(config.Sx[1]-config.Sx[0])/2.
    yg_convol_min = config.Sy[0]-abs(config.Sy[1]-config.Sy[0])/2.
    yg_convol_max = config.Sy[1]+abs(config.Sy[1]-config.Sy[0])/2.
    Ix_convol = numpy.array([xg_convol_min, xg_convol_max])
    Iy_convol = numpy.array([yg_convol_min, yg_convol_max])
    delta_x = config.hx_remap / 2.
    delta_y = config.hy_remap / 2.        
    nqx = abs(xg_convol_max - xg_convol_min)/delta_x
    nqy = abs(yg_convol_max - yg_convol_min)/delta_y
    nqx , nqy = 50, 50
    Nmesh_convol = int(round(max(nqx,nqy),0))    
    [Xgrid_convol, Ygrid_convol] = make_grid_unif(Ix_convol, Iy_convol, Nmesh_convol)
    quadrature_pts_convol = calculsfor_var_modf90.quadrature_pts_trapezoidal_method(Xgrid_convol, Ygrid_convol, len(Xgrid_convol), len(Ygrid_convol))        
    return quadrature_pts_convol , Xgrid_convol, Ygrid_convol
	
