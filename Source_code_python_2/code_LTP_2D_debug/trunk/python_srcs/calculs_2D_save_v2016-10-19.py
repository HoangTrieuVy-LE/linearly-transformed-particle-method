#  !/usr/bin/env python
# - * - codin: utf-8 -*-
# derniere modif le 09-10-15 par Frederique
from __future__ import division
import numpy  #as np
from scipy import *  #as sp

from shapefunction2D import *
from data_code_LTP_2D  import *
from graphes2D import draw_analytical_solution
from graphes2D import fait_series_dessins
from reconstruction_densites_2D import rec_densite_grille_unif_LTP

from _calculsfor_f90 import  calculsfor_ini_modf90
from _calculsfor_f90 import  calculsfor_rec_modf90
from _calculsfor_f90 import  calculsfor_var_modf90


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

#===== Appels des calculs de U =====
def calcule_U_PIC(X, N, t):    
    start1 = time.clock()
    R = rec_densite_grille_unif_sp(X, M, epsilon, Xgrille, Ygrille)     
    Ugrille = calcul_U_sur_grille(R, gradw, Xgrille, Ygrille)
    U = depose_sur_part_sp(Ugrille, X, Xgrille, Ygrille)
    print 'temps de calcul de U', time.clock() - start1  

def calcule_U_part(X, N, t):
    U = zeros((2,N))
    for k in range(N):
        U[:,k] = advection_field(t,X[0,k],X[1,k])
    return U

#===== Calculs sur une matrice 2*2 =====
def calculdet(A):
    return A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]

def inverse(A):
    det = A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
    if det==0:
        print 'Matrice non inversible'
        exit()
    else:
        B = numpy.zeros((2, 2))
        B[0, 0] = A[1, 1]
        B[1, 1] = A[0, 0]
        B[1, 0] = -A[1, 0]
        B[0, 1] = -A[0, 1]
        B = B/det
        return B

def calcul_expo_mat(A):
# approche exp(A), ou A est une matrice 2*2
    I= numpy.zeros((2,2))      
    I[0,0] = 1.
    I[1,1] = 1.    
    E = I + A    
    d = calculdet(E)
    if d>0:
        return E        
    else:
        B = numpy.dot(A,A)
        E += B/2.
        d = calculdet(E)
        if d>0:                    
            return E
        else:
            B = numpy.dot(B,A)
            E += B/6.            
            return E

def MaJ_D(D, X, t, dt) :    
    Dfor, Norm_inf_Dm1, indice_max_norm_Dm1 = MaJ_D_avec_f2py(D,X, t, dt)    
    #Dpy, Norm_inf_Dm1, indice_max_norm_Dm1 = MaJ_D_python(D, X, t, dt)    
    return Dfor, Norm_inf_Dm1, indice_max_norm_Dm1

def MaJ_D_avec_f2py(D, X, t, dt):
    N = len(X[0,:])
    if N==0:
        print 'plus de particules !!!'
    # parameter 1 in 'update_d_schemes' is the number for euler explicit time scheme	
    indice_max_norm_Dm1, Norm_inf_Dm1, Df2py,  = calculsfor_var_modf90.update_d_schemes(X,D, t, dt,1, N)
    indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)            
    return Df2py, Norm_inf_Dm1, indice_max_norm_Dm1
    
def MaJ_D_barenblatt(D, X, M, t, dt,hx_remap, hy_remap):
	N = len(X[0,:])
	if N==0:
		print 'plus de particules !!!'
    # parameter 1 in 'update_d_schemes' is the number for euler explicit time scheme	
    #diff_field_jac = calculsfor_rec_modf90.diffusion_field_jacobian_barenblatt(X, M, D, m_barenblatt, hx_remap, hy_remap, N)
	update_d_scheme = 0
	if D_method_scheme == 'euler_implicit' :
		update_d_scheme = 3
		indice_max_norm_Dm1, Norm_inf_Dm1, Df2py  = calculsfor_var_modf90.update_d_schemes_diffusion_barenblatt(X, D, M, t, dt, update_d_scheme, m_barenblatt, hx_remap, hy_remap, N)    	
		indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)
	elif D_method_scheme  == 'euler_explicit' :
		update_d_scheme = 1
		indice_max_norm_Dm1, Norm_inf_Dm1, Df2py  = calculsfor_var_modf90.update_d_schemes_diffusion_barenblatt(X, D, M, t, dt, update_d_scheme, m_barenblatt, hx_remap, hy_remap, N)    	
		indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)
	elif D_method_scheme  == 'RK2' :
		update_d_scheme = 2
		indice_max_norm_Dm1, Norm_inf_Dm1, Df2py  = calculsfor_var_modf90.update_d_schemes_diffusion_barenblatt(X, D, M, t, dt, update_d_scheme, m_barenblatt, hx_remap, hy_remap, N)    	
		indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)
	else :
		raise ValueError('wrong D_method_scheme value given : '+str(D_method_scheme))
	#for k in range(0, len(X[0,:])) :
	#	print "k , Df2py[:,k] = " , k , Df2py[:,k]
	return Df2py, Norm_inf_Dm1, indice_max_norm_Dm1    

def MaJ_D_barenblatt_mid_scheme(D, X, X_middle, M, t, dt, hx_remap, hy_remap):
	N = len(X[0,:])
	if N==0:
		print 'plus de particules !!!'    
	update_d_scheme = 0
	if time_scheme == 'middle_point' :
		int_time_scheme = 1
		D_middle, indice_max_norm_Dm1, Norm_inf_Dm1, D_new = calculsfor_var_modf90.update_d_barenblatt_mid_scheme(X, X_middle, D, M, t, dt,1, m_barenblatt, hx_remap, hy_remap, N)				
		#~ D_middle, indice_max_norm_Dm1, Norm_inf_Dm1, Dout
		indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)
	else :
		raise ValueError('wrong time_scheme value given : '+str(time_scheme))	
		
	return D_middle, D_new, Norm_inf_Dm1, indice_max_norm_Dm1    
	
	
#______________________________________________________________________#
#
# Compute matrix Jkn with exact reconstruction of density rho  to test 
# against LTP reconstruction of Jkn
#______________________________________________________________________#
def MaJ_D_barenblatt_analytic(D, X, M, t, dt):
	N = len(X[0,:])
	if N==0:
		print 'plus de particules !!!'
	Id = np.zeros([2,2])
	hess_sol_barenblatt = np.zeros([2,2])
	Jk = np.zeros([2,2])
	invJk = np.zeros([2,2])
	Id[0,0] = 1.
	Id[1,1] = 1.
	Norm_inf_Dm1 = 0.
	for k in range(N):        
        #----- calcul de l inverse de Jkn ----- 
		hess_sol_barenblatt = hess_barenblatt(X[0,k],X[1,k],t)
		Jk = Id - 2 * dt * hess_sol_barenblatt
		invJk = inverse(Jk)
        
        #print 'k, invJk : ' , k, invJk
        #----- mise a jour D ----- 
		Dk = D[:, k].reshape((2, 2))
		Dk = numpy.dot(Dk, invJk)
        #-----  calcul de la norme infiny de l'inverse de D ----- 
		a = D[0 , k]
		b = D[1 , k]   
		c = D[2 , k]
		d = D[3 , k]
		det_Dk = abs(a*d-b*c)      
		norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
        #print 'k , norm_inf_Dkm1 , Norm_inf_Dm1 = ' , k, norm_inf_Dkm1, Norm_inf_Dm1            
		if norm_inf_Dkm1 > Norm_inf_Dm1:
			Norm_inf_Dm1 = norm_inf_Dkm1
			indice_max_norm_Dm1 = k
            
		D[:, k] = numpy.reshape(Dk, 4)
        
    #print 'D[Part_suivies[0]]', D[:,Part_suivies[0]]
        
	return D, Norm_inf_Dm1, indice_max_norm_Dm1
		
		

def error_Maj_D_exponential_versus_time_scheme(D, X, t, dt):
    N = len(X[0,:])
    if N==0:
        print 'plus de particules !!!'
    calculsfor_var_modf90.set_flow_params(int(flot), Tmax, alpha)        
    indice_max_norm_Dm1, Norm_inf_Dm1, Df2py,  = calculsfor_var_modf90.update_d(X,D, t, dt, N)
    indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)
    #print "Df2py : " , Df2py
    #print "Norm_inf_Dm1 : " , Norm_inf_Dm1
    #print "indice_max_norm_Dm1 : " , indice_max_norm_Dm1 - 1
    
	# parameter 1 in 'update_d_schemes' is the number for euler explicit time scheme
	# maybe implement 
    indice_max_norm_Dm1_2, Norm_inf_Dm1_2, Df2py_2,  = calculsfor_var_modf90.update_d_schemes(X,D, t, dt,1, N)
    indice_max_norm_Dm1_2 = int(indice_max_norm_Dm1_2 - 1) #so this array makes sense in pythonic world (not sure if useful)
    #print "\nDf2py_2 : " , Df2py_2
    #print "Norm_inf_Dm1_2 : " , Norm_inf_Dm1_2
    #print "indice_max_norm_Dm1_2 : " , indice_max_norm_Dm1_2 - 1
    
    print " \nerror = " , error(Df2py, Df2py_2)
    
    return Df2py, Norm_inf_Dm1, indice_max_norm_Dm1
    
#===== Calcul de la matrice de deformation =====
def MaJ_D_python(D, X, t, dt):
    # on met a jour Dkn par la formule
    # Dk^{n+1}=D_k^n * invJkn
    Norm_inf_Dm1 = 0.
    N = len(X[0,:])
    if N==0:
        print 'plus de particules !!!'
    for k in range(N):        
        #----- calcul de l inverse de Jkn ----- 
        Akn = Dadev(t,X[0,k],X[1,k])        
        invJk = calcul_expo_mat(-Akn*dt)
        #print 'k, invJk : ' , k, invJk
        #----- mise a jour D ----- 
        Dk = D[:, k].reshape((2, 2))
        Dk = numpy.dot(Dk, invJk)
        #-----  calcul de la norme infiny de l'inverse de D ----- 
        a = D[0 , k]
        b = D[1 , k]   
        c = D[2 , k]
        d = D[3 , k]
        det_Dk = abs(a*d-b*c)      
        norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
        #print 'k , norm_inf_Dkm1 , Norm_inf_Dm1 = ' , k, norm_inf_Dkm1, Norm_inf_Dm1            
        if norm_inf_Dkm1 > Norm_inf_Dm1:
            Norm_inf_Dm1 = norm_inf_Dkm1
            indice_max_norm_Dm1 = k
            
        D[:, k] = numpy.reshape(Dk, 4)
        
    #print 'D[Part_suivies[0]]', D[:,Part_suivies[0]]
        
    return D, Norm_inf_Dm1, indice_max_norm_Dm1
    
#needs to add remapping    


# subroutine update_D_finite_differences(Xpart, tab, dt ,hx , hy, indice_max_norm_Dm1, Norm_inf_Dm1, Dout, N_part)

def MaJ_D_FD_scheme(X, D , tab, t, dt ,hx , hy):
	#Dpy, Norm_inf_Dm1py, indice_max_norm_Dm1py = MaJ_D_FD_scheme_avec_python(X, D , tab, t, dt ,hx , hy)
	#print "Norm_inf_Dm1py, indice_max_norm_Dm1py : " , Norm_inf_Dm1py, indice_max_norm_Dm1py , "\n\n"
	#for k in range(0, len(X[0,:])) :
		#print "Dpy[:,k] = " , Dpy[:,k]
	Dfor, Norm_inf_Dm1for, indice_max_norm_Dm1for = MaJ_D_FD_scheme_avec_f2py(X, tab, dt ,hx , hy)
	#print " Norm_inf_Dm1for, indice_max_norm_Dm1for : " , Norm_inf_Dm1for, indice_max_norm_Dm1for , "\n\n"
	#for k in range(0,len(X[0,:])) :
		#print k , " , " , Dpy[:,k] , " , ", Dfor[:,k]
	#	print k , " , " , Dfor[:,k]
	#print "\nerror : " , error(Dpy, Dfor,error_type='absolute_error')
	#return Dpy, Norm_inf_Dm1py, indice_max_norm_Dm1py
	return Dfor, Norm_inf_Dm1for, indice_max_norm_Dm1for
		
def MaJ_D_FD_scheme_avec_f2py(X, tab, dt ,hx , hy) :
	Norm_inf_Dm1 = 0.
	N = len(X[0,:])	
	if N==0:
		print 'plus de particules !!!'	
	
	indice_max_norm_Dm1, Norm_inf_Dm1, Df2py = calculsfor_var_modf90.update_d_finite_differences(X, tab, dt ,hx , hy, N)
	#for k in range(0,len(X[0,:])) :
		#print k , " , " , Df2py[:,k]
	indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)            
	return Df2py, Norm_inf_Dm1, indice_max_norm_Dm1		
	
def MaJ_D_FD_scheme_avec_python(X, D , tab, t, dt ,hx , hy):
	#hx = (Sx[1] - Sx[0]) / float(Nx)  #  distance between particles
	#hy = (Sy[1] - Sy[0]) / float(Ny)
    # on met a jour Dkn par la formule
    # Dk^{n}=invJkn+1 ^{-1}
    #tab[k][0] = indice of particle at north of particle k at t=0
    #tab[k][1] = indice of particle at south of particle k at t=0
    #tab[k][2] = indice of particle at east of particle k at t=0
    #tab[k][3] = indice of particle at west of particle k at t=0
	Norm_inf_Dm1 = 0.
	N = len(X[0,:])
	if N==0:
		print 'plus de particules !!!'
	for k in range(N):        		
		
        #----- calcul de Jkn ----- 
		Jkn = numpy.zeros([2,2])
		ind_N=tab[k][0]
		ind_S=tab[k][1]
		ind_E=tab[k][2]
		ind_W=tab[k][3]
		#print "k ,  ind_N , ind_S , ind_E , ind_W : " , k  , ind_N , ind_S , ind_E , ind_W
		#full West (one particle on East)
		if (int(ind_W) == -1)  & (int(ind_N) == -1) & (int(ind_S) == -1) & (int(ind_E) != -1) :
			Jkn[0,0]=(X[0,ind_E] - X[0,k])/hx
			Jkn[1,0]=(X[1,ind_E] - X[1,k])/hx
			#Jkn[0,1]=(X[0,ind_N] - X[0,ind_S])/hy
			#Jkn[1,1]=(X[1,ind_N] - X[1,ind_S])/hy	
		#full North (one particle on South)
		if (int(ind_W) == -1)  & (int(ind_N) == -1) & (int(ind_E) == -1) & (int(ind_S) != -1) :
			#Jkn[0,0]=(X[0,ind_E] - X[0,k])/hx
			#Jkn[1,0]=(X[1,ind_E] - X[1,k])/hx
			Jkn[0,1]=(X[0,k] - X[0,ind_S])/hy
			Jkn[1,1]=(X[1,k] - X[1,ind_S])/hy
		#full East (one particle on West)
		if (int(ind_E) == -1)  & (int(ind_N) == -1) & (int(ind_S) == -1) & (int(ind_W) != -1) :
			Jkn[0,0]=(X[0,k] - X[0,ind_W])/hx
			Jkn[1,0]=(X[1,k] - X[1,ind_W])/hx
			#Jkn[0,1]=(X[0,ind_N] - X[0,ind_S])/hy
			#Jkn[1,1]=(X[1,ind_N] - X[1,ind_S])/hy
		#full South (one particle on North)
		if (int(ind_W) == -1)  & (int(ind_E) == -1) & (int(ind_S) == -1) & (int(ind_N) != -1) :
			#Jkn[0,0]=(X[0,ind_E] - X[0,k])/hx
			#Jkn[1,0]=(X[1,ind_E] - X[1,k])/hx
			Jkn[0,1]=(X[0,ind_N] - X[0,k])/hy
			Jkn[1,1]=(X[1,ind_N] - X[1,k])/hy
		#between North and South		
		if (int(ind_W) == -1) & (int(ind_E) == -1) & (int(ind_N) != -1) & (int(ind_S) != -1):
			#derivative in x ??
			#Jkn[0,0]=(X[0,ind_E] - X[0,k])/hx
			#Jkn[1,0]=(X[1,ind_E] - X[1,k])/hx
			Jkn[0,1]=(X[0,ind_N] - X[0,ind_S])/(2*hy)
			Jkn[1,1]=(X[1,ind_N] - X[1,ind_S])/(2*hy)
		#between West and East
		if (int(ind_N) == -1) & (int(ind_S) == -1) & (int(ind_W) != -1) & (int(ind_E) != -1):
			#derivative in y ??
			Jkn[0,0]=(X[0,ind_E] - X[0,ind_W])/(2*hx)
			Jkn[1,0]=(X[1,ind_E] - X[1,ind_W])/(2*hx)
			#Jkn[0,1]=(X[0,ind_N] - X[0,ind_S])/hy
			#Jkn[1,1]=(X[1,ind_N] - X[1,ind_S])/hy
        #corner NW
		if (int(ind_N) == -1) & (int(ind_W) == -1) & (int(ind_E) != -1) & (int(ind_S) != -1):
			Jkn[0,0]=(X[0,ind_E] - X[0,k])/hx
			Jkn[1,0]=(X[1,ind_E] - X[1,k])/hx
			Jkn[0,1]=(X[0,k] - X[0,ind_S])/hy
			Jkn[1,1]=(X[1,k] - X[1,ind_S])/hy
		#corner NE
		elif (int(ind_N) == -1) & (int(ind_E) == -1) & (int(ind_W) != -1) & (int(ind_S) != -1):			
			Jkn[0,0]=(X[0,k] - X[0,ind_W])/hx
			Jkn[1,0]=(X[1,k] - X[1,ind_W])/hx
			Jkn[0,1]=(X[0,k] - X[0,ind_S])/hy
			Jkn[1,1]=(X[1,k] - X[1,ind_S])/hy
		#corner SE
		elif (int(ind_S) == -1) & (int(ind_E) == -1) & (int(ind_W) != -1) & (int(ind_N) != -1):
			Jkn[0,0]=(X[0,k] - X[0,ind_W])/hx
			Jkn[1,0]=(X[1,k] - X[1,ind_W])/hx
			Jkn[0,1]=(X[0,ind_N] - X[0,k])/hy
			Jkn[1,1]=(X[1,ind_N] - X[1,k])/hy
		#corner SW
		elif (int(ind_S) == -1) & (int(ind_W) == -1) & (int(ind_N) != -1) & (int(ind_E) != -1):
			Jkn[0,0]=(X[0,ind_E] - X[0,k])/hx
			Jkn[1,0]=(X[1,ind_E] - X[1,k])/hx
			Jkn[0,1]=(X[0,ind_N] - X[0,k])/hy
			Jkn[1,1]=(X[1,ind_N] - X[1,k])/hy
		#north but not in corner
		elif (int(ind_N) == -1) & (int(ind_W) != -1) & (int(ind_E) != -1) & (int(ind_S) != -1):
			Jkn[0,0]=(X[0,ind_E] - X[0,ind_W])/(2*hx)
			Jkn[1,0]=(X[1,ind_E] - X[1,ind_W])/(2*hx)
			Jkn[0,1]=(X[0,k] - X[0,ind_S])/hy
			Jkn[1,1]=(X[1,k] - X[1,ind_S])/hy
		#east but not in corner
		elif (int(ind_E) == -1) & (int(ind_N) != -1) & (int(ind_S) != -1) & (int(ind_W) != -1):
			Jkn[0,0]=(X[0,k] - X[0,ind_W])/hx
			Jkn[1,0]=(X[1,k] - X[1,ind_W])/hx
			Jkn[0,1]=(X[0,ind_N] - X[0,ind_S])/(2*hy)
			Jkn[1,1]=(X[1,ind_N] - X[1,ind_S])/(2*hy)
		#south but not in corner
		elif (int(ind_S) == -1) & (int(ind_W) != -1) & (int(ind_E) != -1) & (int(ind_N) != -1):
			Jkn[0,0]=(X[0,ind_E] - X[0,ind_W])/(2*hx)
			Jkn[1,0]=(X[1,ind_E] - X[1,ind_W])/(2*hx)
			Jkn[0,1]=(X[0,ind_N] - X[0,k])/hy
			Jkn[1,1]=(X[1,ind_N] - X[1,k])/hy
		#west but not in corner
		elif (int(ind_W) == -1) & (int(ind_N) != -1) & (int(ind_S) != -1) & (int(ind_E) != -1):
			Jkn[0,0]=(X[0,ind_E] - X[0,k])/hx
			Jkn[1,0]=(X[1,ind_E] - X[1,k])/hx
			Jkn[0,1]=(X[0,ind_N] - X[0,ind_S])/(2*hy)
			Jkn[1,1]=(X[1,ind_N] - X[1,ind_S])/(2*hy)
		# not in "boundary"
		else :
			Jkn[0,0]=(X[0,ind_E] - X[0,ind_W])/(2*hx)
			Jkn[1,0]=(X[1,ind_E] - X[1,ind_W])/(2*hx)
			Jkn[0,1]=(X[0,ind_N] - X[0,ind_S])/(2*hy)
			Jkn[1,1]=(X[1,ind_N] - X[1,ind_S])/(2*hy)
		
		#print "\nJkn : " , Jkn 
		det_Jkn=Jkn[0,0]* Jkn[1,1] - Jkn[1,0] * Jkn[0,1]
		if abs(det_Jkn) <= 0.000001 :
			print "k , det_Jkn , ind_N , ind_S , ind_E , ind_W : " , k , det_Jkn , ind_N , ind_S , ind_E , ind_W
		#print "det : Jkn " , Jkn[0,0]* Jkn[1,1] - Jkn[1,0] * Jkn[0,1]
        #print 'invJk', invJk 
        #----- mise a jour D ----- 
		Dk = numpy.linalg.inv(Jkn)
        #-----  calcul de la norme infiny de l'inverse de D ----- 
		a = Dk[0,0]
		b = Dk[0,1]
		c = Dk[1,0]
		d = Dk[1,1]
		det_Dk = abs(a*d-b*c)      
		norm_inf_Dkm1 = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk
        
		if norm_inf_Dkm1 > Norm_inf_Dm1:
			Norm_inf_Dm1 = norm_inf_Dkm1
			indice_max_norm_Dm1 = k
            
		D[:, k] = numpy.reshape(Dk, 4)
        
    #print 'D[Part_suivies[0]]', D[:,Part_suivies[0]]
	return D, Norm_inf_Dm1, indice_max_norm_Dm1      
  
#===== TRUCS INUTILES =====          
#def MAJ_I(X,T):
#    xmin=min(X[0,:])
#    xmax=max(X[0,:])
#    ymin=min(X[1,:])
#    ymax=max(X[1,:])    
#    return 0.
#def det_vect(M):
    # pour calcul le determinant d une matrice sous forme array
    #  sans avoir a la redimensionner

#def calcul_nv_poids(M, D):
#    N = len(M)
# #    for k in range(N):
#        det_Dk = det_vect(D[:,k])
#        M[k] = M[k] * det_Dk
#    return M
    
#===== Localisation -> ne sert pas (et non teste pour verifie si pas de bug) =====
def localise(z, Xgrille, Ygrille):
    deltax = abs(Xgrille[1] - Xgrille[0])
    deltax = abs(Xgrille[1] - Xgrille[0])
    deltay = abs(Ygrille[1] - Ygrille[0])
    xmin = min(Xgrille)   
    ymin = min(Ygrille)
    j = int((z[0] - xmin) / deltax)
    i = int((z[1] - ymin ) / deltay) 
    if Xgrille[j]>z[0] :
        print 'erreur dans la localisation'   
    if (j+1 <len(Xgrille) and  Xgrille[j+1]<z[0]):
        print 'erreur dans la localisation'        
    if Ygrille[i]>z[1]:
        print 'erreur dans la localisation'          
    if (i+1 <len(Ygrille) and Ygrille[i+1]<z[1]):
        print 'erreur dans la localisation'   
    return [j,i]

def phin(x,y, p, h, Dk):
    z = numpy.array([x, y])
    u = numpy.dot(Dk, (z -p))/h
    print 'Dk', Dk
    vol = calculdet(Dk)/h
    u1 = u[0]
    u2 = u[1]
    return vol * phi(u1, u2)


def taille_deform_part(D, N):
    Deform_part = numpy.array([0. for i in range(N)])
    for k in range(N):
        a = D[0, k]
        b = D[1, k]   
        c = D[2, k]
        d = D[3, k]
        det_Dk = abs(a*d-b*c)      
        # norme inf de l inverse de Dk
        Deform_part[k] = max(abs(d)+abs(b),abs(c)+abs(a))/det_Dk  
    return Deform_part
      
def volume_particules(D, N, h):
    Vol_part = numpy.array([0. for i in range(N)])
    for k in range(N):
        a = D[0, k]
        b = D[1, k]   
        c = D[2, k]
        d = D[3, k]
        Vol_part[k] = (h*h)/abs(a*d-b*c)      
    return Vol_part

    
    
def diffusion_field(X, N , M , D,hx_remap,hy_remap):
	diffusion_f2py  = diffusion_field_avec_f2py(X, N , M , D, hx_remap, hy_remap)
	#diffusion_python = diffusion_field_avec_python(X, N , M , D, hx_remap, hy_remap)
	#print 'diffusion_f2py : ' , diffusion_f2py
	#print 'diffusion_python : ' , diffusion_python
	return diffusion_f2py
	
def diffusion_field_barenblatt(X, N , M , D,m_barenblatt,hx_remap,hy_remap):
    diffusion_barenblatt_f2py  = calculsfor_rec_modf90.diffusion_field_barenblatt(X, M , D, m_barenblatt, hx_remap, hy_remap, N)
    return diffusion_barenblatt_f2py
    
    
def diffusion_field_avec_python(X, N , M , D,hx_remap,hy_remap):
	diffusion=numpy.zeros([2,N])
	diffusion_num=numpy.zeros([2,N])
	diffusion_denom=numpy.zeros([1,N])
	for k in range(N) :
		if (k % 50 == 0) : print 'k = ',k
		#ADD LINKED LISTS TO NOT LOOP ON EVERY PARTICLES ? OR MOVE TO FOTRAN
		for j in range(N) :
			Dj = D[:, j].reshape((2, 2))
			det_Dj=calculdet(Dj)
			phihn = phi_2d_h( (numpy.dot(Dj,X[:,k]-X[:,j]))[0] , (numpy.dot(Dj,X[:,k]-X[:,j]))[1], hx_remap, hy_remap)
			grad_phihn = numpy.dot(Dj.transpose() , grad_phi_h((numpy.dot(Dj,X[:,k]-X[:,j]))[0] , (numpy.dot(Dj,X[:,k]-X[:,j]))[1], hx_remap, hy_remap))
			diffusion_num[0,k] += M[j] * det_Dj * grad_phihn[0]
			diffusion_num[1,k] += M[j] * det_Dj * grad_phihn[1]
			diffusion_denom[0,k] += M[j] * phihn
			
			
			#a = D[0 , j]
			#b = D[1 , j]   
			#c = D[2 , j]
			#d = D[3 , j]
			#det_Dj = abs(a*d-b*c)
			#phi_x = phi( (a*(X[0,k] - X[0,j]) + b*((X[1,k] - X[1,j]))) / hx_remap )
			#phi_y = phi( (c*(X[0,k] - X[0,j]) + d*((X[1,k] - X[1,j]))) / hy_remap )
			#phi_2D = (phi_x * phi_y) / (hx_remap*hy_remap) 
			
			#derivative_phi_x = (a/hx_remap) * derivative_phi( (a*(X[0,k] - X[0,j]) + b*((X[1,k] - X[1,j])))  / hx_remap ) * phi_y + (c/hy_remap) * derivative_phi( (c*(X[0,k] - X[0,j]) + d*((X[1,k] - X[1,j])))  / hx_remap ) * phi_x
			#derivative_phi_y = (b/hx_remap) * phi_y * derivative_phi( (a*(X[0,k] - X[0,j]) + b*((X[1,k] - X[1,j]))) / hx_remap ) + (d/hy_remap) * phi_x * derivative_phi( (c*(X[0,k] - X[0,j]) + d*((X[1,k] - X[1,j]))) / hy_remap )
			
			#diffusion_num[0,k]+=M[j] * det_Dj/ (hx_remap*hy_remap) * derivative_phi_x
			#diffusion_num[1,k]+=M[j] * det_Dj/ (hx_remap*hy_remap) * derivative_phi_y			
			#diffusion_denom[0,k]+=M[j] * det_Dj/ (hx_remap*hy_remap) * phi_2D
		
		diffusion[:,k]=diffusion_num[:,k] / diffusion_denom[0,k]
	diffusion = sigma*diffusion 
	return diffusion

def diffusion_field_avec_f2py(X, N , M , D,hx_remap,hy_remap):
	#    subroutine diffusion_field(Xpart, M, DD, hx_remap, hy_remap, sigma, diff_field, N_part)
	diffusion= calculsfor_rec_modf90.diffusion_field(X, M , D,hx_remap,hy_remap, sigma,N)		
	#for j in range(N) :	
		#print "j , diffusion(0,j) , diffusion(1,j) = " , j , diffusion[0,j] , diffusion[1,j]
	return diffusion
	
	
	
def number_from_pos(X,N,Xk0,Xk1) :
	pos=-1	
	#print "in number_from_pos : " ,Xk0, Xk1
	for j in range(N) :
		
		if (isclose(X[0,j],Xk0)) & (isclose(X[1,j], Xk1)) :
			pos=j
	#print "pos = " , pos
	return pos

def tab_cardinal_particles(X,N,hx,hy) :
	return tab_cardinal_particles_avec_f2py(X,N,hx,hy)	
	#return tab_cardinal_particles_avec_python(X,N,hx,hy)
	
	
def tab_cardinal_particles_avec_python(X,N,hx,hy) :
	#X = numpy.zeros([2, N])
	print "distance between particles after remapping, hx , hy : " , hy , hy
	tab = {}
	#tab[k] = [N,S,E,W] particles number of particles at North, East, South and West
	#Ntab[k][0] = -1 means there are no particles on North
	for k in range(N) :
		tab[k] = [0,0,0,0]		
		tab[k][0] = number_from_pos(X,N,X[0,k],X[1,k]+hy)
		tab[k][1] = number_from_pos(X,N,X[0,k],X[1,k]-hy)
		tab[k][2] = number_from_pos(X,N,X[0,k]+hx,X[1,k])
		tab[k][3] = number_from_pos(X,N,X[0,k]-hx,X[1,k])
	return tab 

def tab_cardinal_particles_avec_f2py(X,N,hx,hy) :	
	tab_f2py = numpy.zeros([N,4])
	#tab[k,:] = [N,S,E,W] particles number of particles at North, East, South and West	
	tab_f2py = calculsfor_var_modf90.tab_cardinal_particles(X,hx,hy,N)
	
	#change format to python dictionnary : tab[k]=[ind_N,ind_S,ind_E,ind_W]
	tab={}
	for k in range(N) :
		tab[k] = [x for x in tab_f2py[k,:]]						
	return tab 
	
def boundary_particles(tab) :
	boundary_list=[]
	for k in tab.keys() :
		if (tab[k][0] == -1) or (tab[k][1] == -1) or (tab[k][2] == -1) or (tab[k][3] == -1) :
			boundary_list.append(k)
	return boundary_list
			
def solution_barenblatt(Xg,Yg, t) :
	mx = size(Xg)
	my = size(Yg)
	sol = numpy.zeros((mx,my))
	for i in range(mx) :
		for j in range(my) :
			sol[i,j] = barenblatt(Xg[i],Yg[j],t)
	return sol

def compute_solution(Xg,Yg, t) :
	mx = size(Xg)
	my = size(Yg)
	sol = numpy.zeros((mx,my))
	for i in range(mx) :
		for j in range(my) :
			if name_solution == 'barenblatt' :
				sol[i,j] = barenblatt(Xg[i],Yg[j],t)
			elif name_solution == 'SW' : # analytic solution only at t=Tmax
				sol[i,j] = rhoini(Xg[i],Yg[j])
			elif name_solution == 'linear_advection' :
				sol[i,j] = numpy.exp(-alpha*t) * rhoini(Xg[i]*numpy.exp(-alpha*t),Yg[j])
			elif name_solution == 'quadratic_advection' :
				sol[i,j] = (1/(1+t*Xg[i])**2) * rhoini(Xg[i]/(Xg[i]*t+1),Yg[j])
			elif name_solution == 'test_incompressible' :
				sol[i,j] = rhoini(Xg[i]*numpy.exp(-alpha*t), Yg[j]*numpy.exp(alpha*t))
	return sol

def sol_final_SW(Xg,Yg) :
	mx = size(Xg)
	my = size(Yg)
	sol = numpy.zeros((mx,my))
	for i in range(mx) :
		for j in range(my) :
			sol[i,j] = rhoini(Xg[i],Yg[j])
	return sol	

def error(analytic_sol, approx_sol, error_type='L_infinity'): 
	size=analytic_sol.shape[0] * analytic_sol.shape[1]
	vect_analytic_sol = analytic_sol.reshape(size,1)
	vect_approx_sol = approx_sol.reshape(size,1)
	#for i in range(0,analytic_sol.shape[0]) :
		#for j in range(0,analytic_sol.shape[1]) :
			#print "analytic_sol[i][j] , approx_sol[i][j] : " , analytic_sol[i][j] , approx_sol[i][j]
	if error_type == 'L_infinity' :
		error = max(abs(vect_analytic_sol - vect_approx_sol)) / max(abs(vect_analytic_sol))
	elif error_type == 'absolute_error' :
		error = max(abs(vect_analytic_sol - vect_approx_sol))			
	return error

def move_particles(X,N,t,dt) :
    X=move_particles_avec_f2py(X,N,t,dt)        
    return X
    
def move_particles_avec_f2py(X,N,t,dt) :  
	calculsfor_var_modf90.set_flow_params(int(flot), Tmax, alpha)          
	if time_scheme == 'RK2' : 
		time_scheme_number = 2
	elif time_scheme == 'RK4' :
		time_scheme_number = 3
	elif time_scheme =='euler_explicit' :
		time_scheme_number = 1
	else :
		raise ValueError('uncorrect time scheme given')
	Xfor =  calculsfor_var_modf90.move_particles(X, t, dt, time_scheme_number, N)
	return Xfor
    
def move_particles_avec_python(X,N,t,dt) :    
    if time_scheme == 'RK4' : 
        k1=calcule_U_part(X, N, t)
        k2=calcule_U_part(X +(dt/2.)*k1, N, t+(dt/2.))	
        k3=calcule_U_part(X +(dt/2.)*k2, N, t+(dt/2.)	)	
        k4=calcule_U_part(X +dt*k3, N, t+dt )		        
        X[0, :] += (dt/6.) * (k1[0, :] + 2*k2[0, :] + 2*k3[0, :] + k4[0, :])
        X[1, :] += (dt/6.) * (k1[1, :] + 2*k2[1, :] + 2*k3[1, :] + k4[1, :])
    elif time_scheme == 'RK2' : 
        k1=calcule_U_part(X, N, t)
        k2=calcule_U_part(X +(dt/2.)*k1, N, t+(dt/2.))	
        X[0, :] += dt*k2[0,:]
        X[1, :] += dt*k2[1,:]   
    elif time_scheme == 'euler_explicit' :     
		k1 = calcule_U_part(X, N, t)
		X[0, :] += dt*k1[0,:]
		X[1, :] += dt*k1[1,:]
    else : 
		raise ValueError('uncorrect time scheme given')
    return X

def move_particles_diffusion_barenblatt(X, N , M , D,m_barenblatt,hx_remap,hy_remap) :    
	if time_scheme == 'euler_explicit' :
		k1 = - diffusion_field_barenblatt(X, N , M , D,m_barenblatt,hx_remap,hy_remap)  
		X[0, :] += dt*k1[0,:]
		X[1, :] += dt*k1[1,:]	
	elif time_scheme == 'RK2' : 
		k1=-diffusion_field_barenblatt(X, N , M , D,m_barenblatt,hx_remap,hy_remap)  		
		k2=-diffusion_field_barenblatt(X+(dt/2.)*k1, N , M , D,m_barenblatt,hx_remap,hy_remap)
		X[0, :] += dt*k2[0,:]
		X[1, :] += dt*k2[1,:]        	
	elif time_scheme == 'RK4' : 
		k1 = - diffusion_field_barenblatt(X, N , M , D,m_barenblatt,hx_remap,hy_remap)  		
		k2 = - diffusion_field_barenblatt(X + (dt/2.)*k1, N , M , D,m_barenblatt,hx_remap,hy_remap)
		k3 = - diffusion_field_barenblatt(X + (dt/2.)*k2, N , M , D,m_barenblatt,hx_remap,hy_remap)
		k4 = - diffusion_field_barenblatt(X + dt*k3, N , M , D,m_barenblatt,hx_remap,hy_remap)
		X[0, :] += (dt/6.) * (k1[0,:] + 2*k2[0,:] + 2*k3[0,:] + k4[0,:])
		X[1, :] += (dt/6.) * (k1[1,:] + 2*k2[1,:] + 2*k3[1,:] + k4[1,:])
	else : 
		raise ValueError('uncorrect time scheme given')
	
	return X

def errors_reconstruct_ltp_at_particles_for_barenblatt(X, N , M , D,m_barenblatt,hx_remap,hy_remap,t,Xg,Yg) :
    grad_rho_hn=zeros([2, N])
    rho_hn=zeros(N)
    sol_bar=zeros(N)
    grad_sol_bar=zeros([2, N])
    for k in range(0,N) :
        rho_hn[k] = calculsfor_rec_modf90.rec_ltp_at_coordinate(X, M, D, hx_remap, hy_remap,X[:,k], N)
        sol_bar[k]=barenblatt(X[0,k], X[1,k],t)
        grad_rho_hn[:,k] = calculsfor_rec_modf90.rec_ltp_grad_at_coordinate(X, M, D, hx_remap, hy_remap,X[:,k], N)
        grad_sol_bar[:,k]=grad_barenblatt(X[0,k], X[1,k],t)        
        #print "k ,sol_bar[k], rho_hn[k] = " , k ,sol_bar[k], rho_hn[k]
        #print "k, X[:,k] , sol_bar[k], rho_hn[k] , grad_sol_bar[:,k] ,  grad_rho_hn[:,k] = " , k, X[:,k] , sol_bar[k], rho_hn[k] , grad_sol_bar[:,k] ,  grad_rho_hn[:,k]
    
    print "error between sol barenblatt and reconstruct at particle : " ,  max(abs(sol_bar - rho_hn)) / max(abs(sol_bar))
    print "error between grad x sol barenblatt et grad x reconstruct at particle : " ,  max(abs(grad_sol_bar[0,:] - grad_rho_hn[0,:])) / max(abs(grad_sol_bar[0,:]))
    print "error between grad y sol barenblatt et grad y reconstruct at particle : " ,  max(abs(grad_sol_bar[1,:] - grad_rho_hn[1,:])) / max(abs(grad_sol_bar[1,:]))    
    
    mx = size(Xg)
    my = size(Yg)
    grad_x_sol = numpy.zeros((mx,my))
    grad_y_sol = numpy.zeros((mx,my))
    for i in range(mx) :
        for j in range(my) :	
            grad=grad_barenblatt(Xg[i],Yg[j],t)
            grad_x_sol[i,j] = grad[0]
            grad_y_sol[i,j] = grad[1]
    draw_analytical_solution(Xg,Yg,grad_x_sol, 1, t, 'grad_x_sol_barenblatt')            
    draw_analytical_solution(Xg,Yg,grad_y_sol, 1, t, 'grad_y_sol_barenblatt')            
        
    return diffusion_barenblatt_f2py

def errors_reconstruct_ltp_at_particles_for_linear_advection(X, N , M , D,m_barenblatt,hx_remap,hy_remap,t,Xg,Yg) :
    grad_rho_hn=zeros([2, N])    
    rho_hn=zeros(N)
    sol_advec=zeros(N)
    grad_sol_advec=zeros([2, N])    
    for k in range(0,N) :
        rho_hn[k] = calculsfor_rec_modf90.rec_ltp_at_coordinate(X, M, D, hx_remap, hy_remap,X[:,k], N)
        sol_advec[k]=numpy.exp(-alpha*t) * rhoini(X[0,k]*numpy.exp(-alpha*t),X[1,k])
        grad_rho_hn[:,k] = calculsfor_rec_modf90.rec_ltp_grad_at_coordinate(X, M, D, hx_remap, hy_remap,X[:,k], N)
        grad_sol_advec[:,k]=grad_sol_linear_advection(t,X[0,k], X[1,k])                        
        #print "grad_rohini_gaus(X[0,k], X[1,k],t)        : " , X[:,k] , grad_rohini_gaus(X[0,k]*numpy.exp(-t), X[1,k])
        #print "k, X[:,k] , sol_advec[k], rho_hn[k] , grad_sol_advec[:,k] ,  grad_rho_hn[:,k] = " , k, X[:,k] , sol_advec[k], rho_hn[k] , grad_sol_advec[:,k] ,  grad_rho_hn[:,k]
        #print "k, X[:,k] , sol_advec[k], grad_sol_advec[:,k]  = " , k, X[:,k] , sol_advec[k], grad_sol_advec[:,k]
        #print grad_sol_advec[0,k] ,  grad_rho_hn[0,k]
        #print "k, X[:,k] , sol_advec[k], grad_sol_advec[:,k]  = " , k, X[:,k] , sol_advec[k], grad_sol_advec[:,k]
    
    print "error between sol linear advect and reconstruct at particle : " ,  max(abs(sol_advec - rho_hn)) / max(abs(sol_advec))
    print "error between grad x sol linear advect et grad x reconstruct at particle : " ,  max(abs(grad_sol_advec[0,:] - grad_rho_hn[0,:])) / max(abs(grad_sol_advec[0,:]))    
    
    #print "error between grad y sol barenblatt et grad y reconstruct at particle : " ,  max(abs(grad_sol_advec[1,:] - grad_rho_hn[1,:])) / max(abs(grad_sol_advec[1,:]))    
    
    #mx = size(Xg)
    #my = size(Yg)
    #grad_x_sol = numpy.zeros((mx,my))
    #grad_y_sol = numpy.zeros((mx,my))
    #grad=numpy.zeros(2)
    #for i in range(mx) :
        #for j in range(my) :	
            #grad=grad_sol_linear_advection(t,Xg[i],Yg[j])
            #grad_x_sol[i,j] = grad[0]
            #grad_y_sol[i,j] = grad[1]
    #draw_analytical_solution(Xg,Yg,grad_x_sol, 1, t, 'grad_x_sol_lin_advec')            
    #draw_analytical_solution(Xg,Yg,grad_y_sol, 1, t, 'grad_y_sol_lin_advec')            
        
    #return diffusion_barenblatt_f2py
