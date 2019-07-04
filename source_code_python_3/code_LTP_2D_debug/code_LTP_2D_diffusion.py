#!/usr/bin/env python
# ###########################################################################################
# Code LTP 2D
# dernieres modifs 09-10-15 par Frederique
#    - ajout routines pour sous-traitance f2py de reconstruction, initialisation, remapping
# Resolution de \partial_t rho + div(rho u) = 0 avec u donne


#  -*- codin: utf-8 -*-
from __future__ import division
import numpy  #as np
import time
import sys
import os

sys.path.append('./trunk/python_srcs/')


import config
import data_code_LTP_2D

path_to_input_file_name = sys.argv[1]
path_to_output_file_name =  sys.argv[2]
data_code_LTP_2D.initialize_data_parameters(path_to_input_file_name)


import initialisation_2D
import calculs_2D
import shapefunction2D
import graphes2D
import reconstruction_densites_2D
import remapping_2D
import barenblatt_fcts
import diffusion_fcts

# for test, move later
from _calculsfor_f90 import  calculsfor_rec_modf90
from _calculsfor_f90 import  calculsfor_var_modf90

start_code = time.clock()
start_code_bis = time.time()


#=====  Initialisation ====== 
t = config.Tini
npic_part = 0


	
#----- initialisation des positions et poids -----

print 'Initialisation'
if config.redemarrage == 0:
    #X, M = initialize_avec_f2py(Sx, Sy, Nx, Ny)    
    X, M = initialisation_2D.initialize_avec_python(config.Sx, config.Sy, config.Nx, config.Ny, config.rhoini)
    #print 'suppression des particules vides'
    #N, X, M = supprime_particule(X,M)
    N=config.Nini
    #~ if N==config.Nini:
        #~ print 'pas de suppression de particules'
    npic = 0
else:
    X0 = numpy.genfromtxt("out/X0.txt",delimiter=",") 
    X1 = numpy.genfromtxt("out/X1.txt",delimiter=",") 
    M = numpy.genfromtxt("out/M.txt",delimiter=",") 
    N = len(M)
    X = numpy.zeros([2, N])
    t, npic = numpy.genfromtxt("out/data.txt",delimiter=",") 
    #~ print 'redemarrage a t= ', t
    X[0, :] = X0
    X[1, :] = X1

#~ print 'Nombre de particules N=', N 


   
U = numpy.zeros([2, N])
U_barenblatt = numpy.zeros([2, N])
U_ltp = numpy.zeros([2, N])
U_sp = numpy.zeros([2, N])
X_exact = numpy.zeros([2, N])


Mini = M.copy()
numpy.savetxt('out/Mini.txt', Mini)
Xinif = X.copy()
#~ print 'somme des poids', sum(Mini)
Nboucle = 0

#----- Grille de reconstruction -----
[Xgrille, Ygrille] = calculs_2D.make_grid_unif(config.Ix, config.Iy, config.Nmesh)

if (config.D_method == "D_convol") :
    xg_convol_min = config.Sx[0]-abs(config.Sx[1]-config.Sx[0])/2.
    xg_convol_max = config.Sx[1]+abs(config.Sx[1]-config.Sx[0])/2.
    yg_convol_min = config.Sy[0]-abs(config.Sy[1]-config.Sy[0])/2.
    yg_convol_max = config.Sy[1]+abs(config.Sy[1]-config.Sy[0])/2.
    Ix_convol = numpy.array([xg_convol_min, xg_convol_max])
    Iy_convol = numpy.array([yg_convol_min, yg_convol_max])
    delta_x = config.epsilon / 7.
    delta_y = config.epsilon / 7.
    nqx = abs(xg_convol_max - xg_convol_min)/delta_x
    nqy = abs(yg_convol_max - yg_convol_min)/delta_y
    Nmesh_convol = int(round(max(nqx,nqy),0))
    #~ Nmesh_convol = 400
    #~ print "Nmesh_convol = ", Nmesh_convol
    [Xgrid_convol, Ygrid_convol] = calculs_2D.make_grid_unif(Ix_convol, Iy_convol, Nmesh_convol)
    quadrature_pts_convol = calculsfor_var_modf90.quadrature_pts_trapezoidal_method(Xgrid_convol, Ygrid_convol, len(Xgrid_convol), len(Ygrid_convol))
    rho_theory = numpy.zeros([len(Xgrid_convol),len(Ygrid_convol)])
    #~ print quadrature_pts_convol
    #~ raise ValueError("stop")

# ----- Initialisation des supports -----
D = numpy.zeros([4, N]) #matrice de deformation
D_middle = numpy.zeros([4, N])
D_new = numpy.zeros([4, N])
D_old = numpy.zeros([4, N])
D_exact = numpy.zeros([4, N])


D[0,:] = 1.  # coeff (1,1) de Dk^0
D[3,:] = 1.  # coeff (2,2) de Dk^0
D_exact[0,:] = 1.  # coeff (1,1) de Dk^0
D_exact[3,:] = 1.  # coeff (2,2) de Dk^0
D_middle[0,:] = 1.
D_middle[3,:] = 1.
D_new[0,:] = 1.
D_new[3,:] = 1.
D_old[0,:] = 1.
D_old[3,:] = 1.
# D[0,:] -> coef (1,1) de Dkn
# D[1,:] -> coef (1,2) de Dkn
# D[2,:] -> coef (2,1) de Dkn
# D[3,:] -> coef (2,2) de Dkn

J_exact = numpy.zeros([2,2])
inv_J_exact = numpy.zeros([2,2])
Dk_exact_old = numpy.zeros([2,2])
Dk_exact = numpy.zeros([2,2])
Dk_exact_new = numpy.zeros([2,2])
Dk_exact_test = numpy.zeros([2,2])
Save_det_Dk = numpy.zeros(N);
delta_X = numpy.zeros(N);


R_ltp = reconstruction_densites_2D.rec_densite_grille_unif_LTP(X, M, D, Xgrille, Ygrille)
graphes2D.fait_series_dessins(Xgrille, Ygrille, R_ltp, 0, t, 'LTP')    

nqx = len(Xgrille)
nqy = len(Ygrille)
solution = numpy.zeros([nqx, nqy])
for i in range(0,nqx) :
    for j in range(0,nqy) :
        solution[i,j] = diffusion_fcts.diffusion_init(Xgrille[i],Ygrille[j])

graphes2D.draw_analytical_solution(Xgrille,Ygrille,solution, npic, config.Tini, label=config.name_solution)


#relative_error=calculs_2D.error(solution , R_ltp)
#print "relative_error = " , relative_error

#~ raise ValueError("STOP")

if config.D_method == 'implicit' :	
	#tab=tab_cardinal_particles(X,N,hx_remap,hy_remap)
	tab = calculsfor_var_modf90.tab_cardinal_particles(X,config.hx_remap,config.hy_remap,N)	
	
#======  Boucle en temps ===== 
#~ print "time taken to complete all initialisations : " , time.time() - start_code_bis


print 'Start'
while (round(t,13) < (round((config.Tmax),13)) ) : 
#~ while (abs(t-config.Tmax) < 0.00000001)  :    	
    if (config.method == 'PIC') or (config.method == 'LTPIC'):
        U = calcule_U_PIC(X, M, D, Xgrille, Ygrille)
	
    #~ if ((Nboucle % int(Nboucle_total/4)) ==0) : print 't=', t # on affiche t achaque quart du run
    #~ print "t= " , t
    
    # ----- Trace densites reconstruites et enregistrement -----
    if (Nboucle % config.param == 0) :
    #~ if ((Nboucle % int(Nboucle_total/4)) ==0) : 
        
        if config.method == 'LTP':
            start2 = time.clock()
            time1= time.time()
            R_ltp = reconstruction_densites_2D.rec_densite_grille_unif_LTP(X, M, D, Xgrille, Ygrille)            
            #~ print 'temps de reconstruction de rho avec LTP ' , time.time() - time1
            npicf_ltp =  npic
            graphes2D.fait_series_dessins(Xgrille, Ygrille, R_ltp, npicf_ltp, t, 'LTP')    
    
            R_method = R_ltp

        if (config.method =='part') or (config.trace_particules == 'oui'):
            #~ print "trace particles : "
            graphes2D.traceparticules(X, t, npic)

        if (config.method =='SP'):
            start3 = time.clock()
            R_sp = reconstruction_densites_2D.rec_densite_grille_unif_sp(X, M, config.epsilon, Xgrille, Ygrille)
            #~ print 'temps de reconstruction de rho avec SP', time.clock() - start3 
            npicf_sp = npic
            graphes2D.fait_series_dessins(Xgrille, Ygrille, R_sp, npicf_sp, t, 'SP')      
            R_method = R_sp
		        
        if (config.draw_sol == 'oui' ) :
            solution=calculs_2D.compute_solution(Xgrille,Ygrille, t)            
            graphes2D.draw_analytical_solution(Xgrille,Ygrille,solution, npic, t, label=config.name_solution)
            indic_y=int(len(Ygrille)/2)
    
            relative_error=calculs_2D.error(solution , R_method)
            #~ print "\n\nrelative L-infini error with N = " , N, ", at t =  " , t , 'is : '  , relative_error
            
            #------------------------------       

        if (config.save_data == 'oui') :
            numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_solution.txt', solution)
            numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_R_'+str(config.method)+'D_method'+str(config.D_method)+'_mat.txt', R_method)
            numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_Xgrille.txt', Xgrille)
            numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_Ygrille.txt', Ygrille)
        
        #~ print 'Nplot = ', npic
        #~ npic += 1    
    		        
    #__________________________________________________________________#
    #
    ### SP Method : - Compute U(t,X) with gradient_SP
    #                - Move particles with this U(t,X)
    #__________________________________________________________________#
    if (config.method == 'SP') :
        U = - calculsfor_rec_modf90.diffusion_field_barenblatt_sp(X, M, config.m_barenblatt,config.epsilon, N)
        #~ U = calculsfor_rec_modf90.diffusion_field_barenblatt_sp(X, M, config.m_barenblatt,config.hx_remap, N)
        X[0,:] += config.dt * U[0,:]
        X[1,:] += config.dt * U[1,:]
        
        #~ for k in range(0,N) :
            #~ print U[0,k], U[1,k]
        size = U.shape[0]*U.shape[1]
        print "t, error max L infty U = " , t, max(abs(U.reshape(size)))
        #~ myfile=open('./particles_sp.debug','a')
        #~ for k in range(0,N) :
            #~ myfile.write(str(k)+"\t"+str(X[0,k])+"\t"+str(X[1,k])+"\t"+str(U[0,k])+"\t"+str(U[1,k])+'\n')
        #~ myfile.close()			
        Nboucle += 1
        t += config.dt		
        
    if (config.problem == "diffusion") and (config.method == 'SP') :
        if (config.time_scheme == 'euler_explicit') : 
            U = - calculsfor_rec_modf90.diffusion_field_sp(X, M, config.epsilon, N)
            X[0,:] += config.dt * U[0,:]
            X[1,:] += config.dt * U[1,:]
        
            Nboucle += 1
            t += config.dt		
        if (config.time_scheme == 'middle_point') : 
            raise ValueError("TO DO")
            # HERE
    
    if (config.problem == "diffusion") and (config.method == 'LTP') :
        if (config.time_scheme == 'euler_explicit') : 
            start = time.time()
            print "t = " , t
            D_old = D.copy()
            X_old = X.copy()                       
            
            
            tmp = time.time()            
            U = - calculsfor_rec_modf90.diffusion_field_ltp(X, M , D_old, config.hx_remap, config.hy_remap, N) 
            print "time for U = " , time.time() - tmp
            X[0,:] += config.dt * U[0,:]
            X[1,:] += config.dt * U[1,:]
                                    
            # update deformations of particles 
            tmp = time.time()
            update_d_scheme = 1
            indice_max_norm_Dm1, Norm_inf_Dm1, Df2py  = calculsfor_var_modf90.update_d_diffusion(X_old, D_old, M, t, config.dt, update_d_scheme, config.hx_remap, config.hy_remap, N)            
            print "time for D = " , time.time() - tmp
            D = Df2py.copy()
            indice_max_norm_Dm1 = int(indice_max_norm_Dm1 - 1) #so this array makes sense in pythonic world (not sure if useful)                                    
            #~ R_ltp = reconstruction_densites_2D.rec_densite_grille_unif_LTP(X, M, D, Xgrille, Ygrille)            
            #~ R_theo = calculs_2D.compute_solution(Xgrille, Ygrille, t  - config.Tini)
            #~ 
            #~ relative_error=calculs_2D.error(R_theo , R_ltp)
            #~ print "relative_error = " , relative_error            
            t += config.dt
            npic += 1
            Nboucle +=1
            end = time.time()
            print "time for one iteration = " , end - start
            #~ raise ValueError("STOP")
            
            
        if (config.time_scheme == 'middle_point') : 
            raise ValueError("TO DO")
    
    if (config.method == 'analytic_debug') :		    
        if (config.time_scheme == 'euler_explicit') :             
            D_old = D.copy()
            X_old=X.copy()                    
                        
            D, Norm_inf_Dm1, indice_max_norm_Dm1 = calculs_2D.update_d_hessian_analytic(X_old, N, D_old, config.m_barenblatt, t, config.dt, config.hx_remap, config.hy_remap)                                                                
            
            for k in range(0,N) :
                U_barenblatt[:,k] = -config.m_barenblatt * barenblatt_fcts.grad_barenblatt(X_old[0,k], X_old[1,k], t)            
            U=U_barenblatt
            X[0, :] += config.dt * U[0,:]
            X[1, :] += config.dt * U[1,:]
                        
            for kk in range(0,N) :
                #In : X(s) ,s, t. Returns (t/s) ** 1/4 Xs if |Xs|^2 < 16 sqrt(s)  Xs otherwise
                X_exact[:,kk] = barenblatt_fcts.flow_barenblatt(X_old[:,kk], t, t + config.dt)                
                
                D_exact[0,kk] = barenblatt_fcts.jacobian_flow_barenblatt(X_exact[:,kk], t + config.dt, config.Tini)[0][0]
                D_exact[1,kk] = barenblatt_fcts.jacobian_flow_barenblatt(X_exact[:,kk], t + config.dt, config.Tini)[0][1]
                D_exact[2,kk] = barenblatt_fcts.jacobian_flow_barenblatt(X_exact[:,kk], t + config.dt, config.Tini)[1][0]
                D_exact[3,kk] = barenblatt_fcts.jacobian_flow_barenblatt(X_exact[:,kk], t + config.dt, config.Tini)[1][1]                                                
                
            error_particles_positions = calculs_2D.relative_error_norm_max_all_particles(X_exact, X)
            error_particles_deformations = calculs_2D.relative_error_norm_max_all_matrices(D_exact, D)
                        
            list_max_D_and_D_exact = []
            list_max_D_ex = []
            list_max_X_and_X_exact = []
            list_max_X_ex = []
            count_X=0
            count_D=0
            for k in range(0,N) : 
            
                if not( (abs(D_exact[0,k] + 10000) < 0.0001) and (abs(D_exact[1,k] - (-10000)) <= 0.0001) and (abs(D_exact[2,k] - (-10000)) <= 0.0001) and (abs(D_exact[3,k] - (-10000)) <= 0.0001) ) :
                    #~ print D_exact[:,k]
                    absolute_error_D_and_D_ex_k = max( abs(D_exact[0,k] - D[0,k]) + abs(D_exact[1,k] - D[1,k]), abs(D_exact[2,k] - D[2,k]) + abs(D_exact[3,k] - D[3,k]))
                    absolute_error_D_ex = max( abs(D_exact[0,k]) + abs(D_exact[1,k]), abs(D_exact[2,k]) + abs(D_exact[3,k]))
                    list_max_D_and_D_exact.append(absolute_error_D_and_D_ex_k)                
                    list_max_D_ex.append(absolute_error_D_ex)
                    count_D+=1;
                if not( (abs(X_exact[0,k] - (-10000)) <= 0.0001) and (abs(X_exact[1,k] - (-10000)) <= 0.0001) ) :
                    #~ print X_exact[0,k] , X_exact[1,k]
                    absolute_error_X_and_X_ex_k = max (abs(X[:,k] - X_exact[:,k]))
                    absolute_error_X = max (abs(X[:,k]))
                    list_max_X_ex.append(absolute_error_X)
                    list_max_X_and_X_exact.append(absolute_error_X_and_X_ex_k)
                    count_X+=1                                        
                    
            print "number of particles counted for X errors : " , count_X , "\n"
            print "number of particles counted for D errors : " , count_D , "\n"
            
            #~ print error_particles_positions , error_particles_deformations
            
            error_particles_positions = max(list_max_X_and_X_exact) / max(list_max_X_ex)
            error_particles_deformations = max(list_max_D_and_D_exact) / max(list_max_D_ex)
            
            print error_particles_positions , error_particles_deformations
            
            t += config.dt
            Nboucle += 1
            
        elif (config.time_scheme == 'middle_point') :
            print "H3Y"
            D_old = D.copy()
            X_old=X.copy()
            X, D, Norm_inf_Dm1, indice_max_norm_Dm1 = calculs_2D.update_d_and_x_middle_point_gradient_and_hessian_analytic_v2(X_old, M, N, D_old, config.m_barenblatt, t, config.dt, config.hx_remap, config.hy_remap)                            
            #~ X, D, Norm_inf_Dm1, indice_max_norm_Dm1 = calculs_2D.update_d_and_x_middle_point_gradient_and_hessian_analytic(X_old, M, N, D_old, config.m_barenblatt, t, config.dt, config.hx_remap, config.hy_remap)                            
            for kk in range(0,N) :
                #Xs ,s, t. Returns (t/s) ** 1/4 Xs if |Xs|^2 < 16 sqrt(s)  Xs otherwise
                X_exact[:,kk] = barenblatt_fcts.flow_barenblatt(X_old[:,kk], t, t + config.dt)
                                        
                #Xs ,s, t. Returns (t/s) ** 1/4 Id si |Xs|^2 < 16 sqrt(s) Id otherwise
                #We call jacobian_flow_barenblatt( X(tn), tn , config.Tini)
                D_exact[0,kk] = barenblatt_fcts.jacobian_flow_barenblatt(X_exact[:,kk], t + config.dt, config.Tini)[0][0]
                D_exact[1,kk] = barenblatt_fcts.jacobian_flow_barenblatt(X_exact[:,kk], t + config.dt, config.Tini)[0][1]
                D_exact[2,kk] = barenblatt_fcts.jacobian_flow_barenblatt(X_exact[:,kk], t + config.dt, config.Tini)[1][0]
                D_exact[3,kk] = barenblatt_fcts.jacobian_flow_barenblatt(X_exact[:,kk], t + config.dt, config.Tini)[1][1]
                            
            error_particles_positions = calculs_2D.relative_error_norm_max_all_particles(X_exact, X)
            error_particles_deformations = calculs_2D.relative_error_norm_max_all_matrices(D_exact, D)
            #~ error_particles_positions = calculs_2D.absolute_error_norm_max_all_particles(X_exact, X)
            #~ error_particles_deformations = calculs_2D.absolute_error_norm_max_all_matrices(D_exact, D)
            print error_particles_positions , error_particles_deformations
            
            t += config.dt
            Nboucle += 1
#__________________________________________________________________#
    #
    ### LTP with direct Dk computation method (finite differences) :
    #   - Compute flow U(t,X) with LTP reconstruction of gradient 
    #     and move particles
    #   - Compute the new deformation Dkn matrix with jacobian of U(t,X)
    #     by Finite Differences
    #__________________________________________________________________#    
    if (config.D_method == 'implicit') :
        ### barenblatt diffusion reconstruction
        X_old = X.copy()
        D_old = D.copy()
        #~ X=move_particles_diffusion_barenblatt(X, N , M , D,config.m_barenblatt,config.hx_remap,config.hy_remap) 		
        U = - calculsfor_rec_modf90.diffusion_field_barenblatt(X, M , D_old, config.m_barenblatt, config.hx_remap, config.hy_remap, N)
        X[0, :] += config.dt*U[0,:]
        X[1, :] += config.dt*U[1,:]	
        
        
        D, Norm_inf_Dm1, indice_max_norm_Dm1 = calculs_2D.MaJ_D_FD_scheme(X, D_old ,tab, t, config.dt , config.hx_remap ,config.hy_remap)
        #~ print 'max des normes des Dk', Norm_inf_Dm1
        #~ print 'indice_max_norm_Dm1', indice_max_norm_Dm1
        
        if numpy.isinf(Norm_inf_Dm1) | numpy.isnan(Norm_inf_Dm1) : raise ValueError("Found Inf or NaN, Norm_inf_Dm1 = "+str(Norm_inf_Dm1)+"\n")
        
        Nboucle += 1
        t += config.dt
    #__________________________________________________________________#
    #
    ### LTP with incremental Dk computation method by matrix exponential
    #   approximation of Jkn = inv(Dkn):
    #   - Compute flow U(t,X) with LTP reconstruction of gradient 
    #     and move particles
    #   - Compute the new deformation Dkn matrix with jacobian of U(t,X)
    #     by Finite Differences
    #__________________________________________________________________#    
    if (config.D_method == 'explicit') :
        if (config.time_scheme == 'euler_explicit') :     
            # !!! NORMAL VERSION : 
            D_old = D.copy()
            X_old = X.copy()         
            D, Norm_inf_Dm1, indice_max_norm_Dm1 = calculs_2D.MaJ_D_barenblatt(D_old, X_old, M, t, config.dt, config.hx_remap, config.hy_remap)
            print 'max des normes des Dk', Norm_inf_Dm1
            print 'indice_max_norm_Dm1', indice_max_norm_Dm1            
            if numpy.isinf(Norm_inf_Dm1) | numpy.isnan(Norm_inf_Dm1) : raise ValueError("Found Inf or NaN. Norm_inf_Dm1 , t = "+str(Norm_inf_Dm1)+', '+str(t)+"\n")            
            
            U = - calculsfor_rec_modf90.diffusion_field_barenblatt(X, M , D_old, config.m_barenblatt, config.hx_remap, config.hy_remap, N)
            X[0, :] += config.dt * U[0,:]
            X[1, :] += config.dt * U[1,:]										                                               
            t += config.dt
            Nboucle += 1            
                        
            

        if (config.time_scheme == 'middle_point') :            
            X_old = X.copy()            
            D_old = D.copy()
            
            U = - calculsfor_rec_modf90.diffusion_field_barenblatt(X_old, M , D_old, config.m_barenblatt, config.hx_remap, config.hy_remap, N)
            X_middle = X + (config.dt/2.) * U   
            D_middle, D_new,  Norm_inf_Dm1, indice_max_norm_Dm1 = calculs_2D.MaJ_D_barenblatt_mid_scheme(D_old, X, X_middle, M, t, config.dt, config.hx_remap, config.hy_remap)                        
            U_middle = - calculsfor_rec_modf90.diffusion_field_barenblatt(X_middle, M , D_middle, config.m_barenblatt, config.hx_remap, config.hy_remap, N)
            X_new = X + config.dt * U_middle
            
            D = D_new
            X = X_new
            
            print 'max des normes des Dk', Norm_inf_Dm1
            print 'indice_max_norm_Dm1', indice_max_norm_Dm1            
            if numpy.isinf(Norm_inf_Dm1) | numpy.isnan(Norm_inf_Dm1) : raise ValueError("Found Inf or NaN. Norm_inf_Dm1 , t = "+str(Norm_inf_Dm1)+', '+str(t)+"\n")
            
            Nboucle += 1
            t += config.dt
    
    #__________________________________________________________________#
    #
    ### LTP hybrid : DK^n+1 is computed with Hess[\rho ^\epsilon _{LTP}]
    #__________________________________________________________________#
    if (config.D_method == 'hybrid') :        
        if (config.time_scheme == 'euler_explicit') :         
            #~ print "D_method = " , config.D_method
            D_old = D.copy()
            D, Norm_inf_Dm1, indice_max_norm_Dm1 = calculs_2D.MaJ_D_barenblatt(D_old, X, M, t, config.dt, config.epsilon, config.epsilon)
            
            if numpy.isinf(Norm_inf_Dm1) | numpy.isnan(Norm_inf_Dm1) : raise ValueError("Found Inf or NaN. Norm_inf_Dm1 , t = "+str(Norm_inf_Dm1)+', '+str(t)+"\n")
            
            #~ print 'max des normes des Dk', Norm_inf_Dm1
            #~ print 'indice_max_norm_Dm1', indice_max_norm_Dm1        
            
            diff_field = calculsfor_rec_modf90.diffusion_field_barenblatt(X, M , D_old, config.m_barenblatt, config.hx_remap, config.hy_remap, N)
            X[0, :] += config.dt * ( -diff_field[0,:] )
            X[1, :] += config.dt * ( -diff_field[1,:] )
            
            ### QUICK TEST
            #~ print "t = ", t 
            #~ for k in range(0,N) :
                #~ print k, D[0,k] , D[1,k] , D[2,k] , D[3,k]
            
            
            t += config.dt
            Nboucle += 1    
            
        if (config.time_scheme == 'middle_point') :
            #~ print "D_method = " , config.D_method
            #~ print "time_scheme = " , config.time_scheme
            U = - calculs_2D.diffusion_field_barenblatt(X, N , M , D, config.m_barenblatt, config.hx_remap, config.hy_remap)
            X_middle = X + (config.dt/2.) * U
            D_old = D
            
            D_middle, D_new,  Norm_inf_Dm1, indice_max_norm_Dm1 = calculs_2D.MaJ_D_barenblatt_mid_scheme(D_old, X, X_middle, M, t, config.dt, config.epsilon, config.epsilon)
            
            if numpy.isinf(Norm_inf_Dm1) | numpy.isnan(Norm_inf_Dm1) : raise ValueError("Found Inf or NaN, Norm_inf_Dm1 = "+str(Norm_inf_Dm1)+', '+str(t)+"\n")
            
            #~ print 'max des normes des Dk', Norm_inf_Dm1
            #~ print 'indice_max_norm_Dm1', indice_max_norm_Dm1
            
            U_middle = -calculs_2D.diffusion_field_barenblatt(X_middle, N , M , D_middle, config.m_barenblatt, config.hx_remap, config.hy_remap)
            
            X_new = X + config.dt * U_middle
            
            D = D_new
            X = X_new
            Nboucle += 1
            t += config.dt
    
    
    #__________________________________________________________________#
    #
    ### LTP and convolution with other test function for jacobian 
    #__________________________________________________________________#
    if (config.D_method == 'D_convol') :
        X_old = X.copy()
        D_old = D.copy()
        print "epsilon = " , config.epsilon
        #~ D, Norm_inf_Dm1, indice_max_norm_Dm1 = calculs_2D.MaJ_D_convolution_barenblatt(D, X, M, t, config.dt, config.hx_remap, config.hy_remap, config.epsilon, Xgrid_convol, Ygrid_convol, quadrature_pts_convol)
        #~ D, Norm_inf_Dm1, indice_max_norm_Dm1 = calculs_2D.MaJ_D_barenblatt(D_old, X, M, t, config.dt, config.epsilon, config.epsilon)
        D, Norm_inf_Dm1, indice_max_norm_Dm1 = calculs_2D.MaJ_D_barenblatt(D_old, X, M, t, config.dt, config.hx_remap, config.hy_remap)
        print 'max des normes des Dk', Norm_inf_Dm1
        print 'indice_max_norm_Dm1', indice_max_norm_Dm1
        U = - calculsfor_rec_modf90.diffusion_field_barenblatt_convolution(X, M , D_old, config.m_barenblatt, config.hx_remap, config.hy_remap, config.epsilon, Xgrid_convol, Ygrid_convol, quadrature_pts_convol, len(Xgrid_convol), len(Ygrid_convol), N)
        
        X[0, :] += config.dt * U[0,:]
        X[1, :] += config.dt * U[1,:]
        for ii in range(0,N) : U_barenblatt[:,ii] = -2*barenblatt_fcts.grad_barenblatt(X[0,ii],X[1,ii],t)
        U_relative_err = calculs_2D.error(U_barenblatt, U, "L_infinity")
        print "U_relative_err = " , U_relative_err
        
        for kk in range(0,N) :
            print "python : kk , U[0,kk] , U[1,kk] = ", kk , U[0,kk] , U[1,kk] , U_barenblatt[0,kk] , U_barenblatt[1,kk]
        
        rho_ltp_conv_grid = calculsfor_rec_modf90.rec_ltp_on_convolution_grid(X, M, D, Xgrid_convol, Ygrid_convol, config.hx_remap, config.hy_remap, len(Xgrid_convol), len(Ygrid_convol), N)
        for ix in range(0,len(Xgrid_convol)) :
            for jy in range(0,len(Ygrid_convol)) :
                rho_theory[ix,jy] = barenblatt_fcts.barenblatt(Xgrid_convol[ix],Ygrid_convol[jy],t)
                print Xgrid_convol[ix] , Ygrid_convol[jy] , rho_theory[ix,jy] , rho_ltp_conv_grid[ix,jy]
        
        rho_conv_relative_err = calculs_2D.error(rho_theory, rho_ltp_conv_grid,"L_infinity")
        print "rho_conv_relative_err = " , rho_conv_relative_err
        Nboucle += 1
        t += config.dt
        
        if numpy.isinf(Norm_inf_Dm1) | numpy.isnan(Norm_inf_Dm1) : raise ValueError("Found Inf or NaN. Norm_inf_Dm1 , t = "+str(Norm_inf_Dm1)+', '+str(t)+"\n")
        
        
        

    #---- trace forme particule -----
    if  (Nboucle % config.param == 0) and (len(config.Part_suivies)>0) and ((config.method=='LTP') or (config.method == 'part_avec_suivi_vol')):
         #  reconstruction de la forme des particules suivis
         #dessine_forme_particule(D, M, X, h, Part_suivies, npic_part)
         for k in config.Part_suivies:
             Xgk, Ygk, R_k = rec_forme_particule(k, X[:,k], D[:, k], M[k], 100)
             rep2d_contourf(Xgk, Ygk, R_k, npic, t, 'Shape_of_particle_'+str(k))
                #plt.close("all")
         npic_part +=1
	

    
    #---- Mise a jour eventuelle de la grille -----
    config.Ix, config.Iy, modif_grille = calculs_2D.MAJ_grille(X,config.Ix,config.Iy)
    if modif_grille == 1:
        [Xgrille, Ygrille] = calculs_2D.make_grid_unif(config.Ix, config.Iy, config.Nmesh)

    # ---- remapping -----------------------
    if (config.method == 'LTP'):
        if (Norm_inf_Dm1 >config.radius_remap) and (config.indic_remapping == 'yes'):
            print '\nRemapping!\n'
            if D_method == 'implicit' :
                X, M, D  = remapping_avec_python_test(config.hx_remap, config.hy_remap, X, M, D, t) 
                N = len(M)   
                tab = calculsfor_var_modf90.tab_cardinal_particles(X,config.hx_remap,config.hx_remap,N)	
            if D_method =='explicit' :
                X, M, D  = remapping_avec_python_test(config.hx_remap, config.hy_remap, X, M, D, t) # remapping that also returns remapping grid step				
                N = len(M)   
            print 'nouveau nombre de particules', N

#    
#     ---- enregistrement des tableaux -----    
#    ~ if ((Nboucle % int(Nboucle_total/4)) ==0) : 
    if ((Nboucle % config.param) ==0 and (config.method != "LTP")) : 
        if (config.save_data == 'oui') :
            time_tabs=time.time()
            numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_M.txt', M)
            numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_X0.txt', X[0,:])
            numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_X1.txt', X[1,:])
            numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_data.txt', numpy.array([t, npic]))
            numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D0.txt', D[0,:])
            numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D1.txt', D[1,:])
            numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D2.txt', D[2,:])
            numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D3.txt', D[3,:])
            print "time to save data in txt files : " , time.time()-time_tabs
    if (config.problem == "diffusion") and (config.method == "LTP_casier"):
        if (config.time_scheme == 'euler_explicit') : 
            break
            
#===== fin de la boucle en temps ======
#=========================================== DATA FORTRAN FOR LTP METHOD=======================================================
D_old = D.copy()
X_old = X.copy()                

#Nx, Ny, hx, hy, dt, deg_spline, radius_phi, bool_phi_radial 
sizefile = open("trunk/fortran_srcs/size.data",'w')
sizefile.write(str(config.Nx)+'\n')
sizefile.write(str(config.Ny)+'\n')
sizefile.write(str(config.hx_remap)+'\n')
sizefile.write(str(config.hy_remap)+'\n')
sizefile.write(str(config.dt)+'\n')
sizefile.write(str(config.Tini)+'\n')
sizefile.write(str(config.Tmax)+'\n')
sizefile.write(str(1)+'\n')  #update_d_scheme = 1
sizefile.write(str(config.deg_spline)+'\n')
sizefile.write(str(config.radius_phi)+'\n')
sizefile.write(str(config.bool_phi_radial))
sizefile.close()


#NOTEEEE: Lx1 = mesh(1,1), Lx2 = mesh(1,2),  Ly1 = mesh(2,1) , Ly2 = mesh(2,2)
mesh = numpy.array([config.Ix,config.Iy]) 
#print('Im here',mesh)
mesh.T.tofile('trunk/fortran_srcs/mesh.bin')
# Coordinates, Deformation Matrix D, and matrix M
X.T.tofile('trunk/fortran_srcs/coords4fortran.bin')
#print(X)
D.T.tofile('trunk/fortran_srcs/deformmatrix4fortran.bin')  
M.T.tofile('trunk/fortran_srcs/matrixM4fortran.bin')  
#print(numpy.shape(U))   
#print(U[:,1:10])                 
U.T.tofile('trunk/fortran_srcs/velocity4fortran.bin')
#print(M)


#=========================================================================================================================================

if (config.problem == "diffusion") and (config.method == "LTP_casier"):
        if (config.time_scheme == 'euler_explicit') : 
            start = time.time()
#            print "t = " , t   -g -fcheck=all -Wall
                             
            cmd = 'mpif90  -fopenmp trunk/fortran_srcs/tri_casier_method.f90 -o outfile \
                    trunk/fortran_srcs/launchfortran.o \
                    trunk/fortran_srcs/mod_particle_2D.o \
                    trunk/fortran_srcs/MPI_2D_spatialisation.o \
                    trunk/fortran_srcs/MPI_2D_structures.o \
                    trunk/fortran_srcs/calculsfortran_rec.o \
                    trunk/fortran_srcs/calculs_2D.o \
                    trunk/fortran_srcs/calculsfortran_ini.o \
                    trunk/fortran_srcs/pack.o \
                     trunk/fortran_srcs/jacobi_method.o   '
            os.system(str(cmd)) 
#            print ("launch : " , cmd )c
            cmd = 'mpiexec -n 9 ./outfile'
            os.system(str(cmd))
            
            Xread = numpy.fromfile('trunk/fortran_srcs/coords4fortran.bin')
            X = Xread.reshape((config.Nini,2)).T
#            print(X[:,1:10])
            
            
#        if (Norm_inf_Dm1 >config.radius_remap) and (config.indic_remapping == 'yes'):
#            print '\nRemapping!\n'
#            if D_method == 'implicit' :
#                X, M, D  = remapping_avec_python_test(config.hx_remap, config.hy_remap, X, M, D, t) 
#                N = len(M)   
#                tab = calculsfor_var_modf90.tab_cardinal_particles(X,config.hx_remap,config.hx_remap,N)	
#            if D_method =='explicit' :
#                X, M, D  = remapping_avec_python_test(config.hx_remap, config.hy_remap, X, M, D, t) # remapping that also returns remapping grid step				
#                N = len(M)   
#            print 'nouveau nombre de particules', N
#            
#            # Make no sense these 2 following lines       
#            indice_max_norm_Dm1 = 0
#            Norm_inf_Dm1= 0
#            
#            
#            t += config.dt
#            npic += 1
#            Nboucle +=1
#            end = time.time()
#            print "time for one iteration = " , end - start       



if (config.method == 'LTP') or (config.method == 'part_avec_suivi_vol'): 
     start4 = time.clock()     
     #if D_method == 'explicit' : 
		##D, Norm_inf_Dm1, indice_max_norm_Dm1 = MaJ_D(D, X, t, dt)
		#D, Norm_inf_Dm1, indice_max_norm_Dm1 = MaJ_D_barenblatt(D, X, M, t, dt, hx_remap, hy_remap)
 
#===== desssins finaux =====
if config.animmp4 == 1:
    nom ='Density_'+config.method+'_'
    graphes2D.faitlefilm(nom)
       
#traceparticules(X,t,npic)

name = config.method + '_N='+str(N)+'-Nmesh='+str(config.Nmesh) 
name_analytical = config.name_solution + '_N='+str(N)+'-Nmesh='+str(config.Nmesh)

if (config.name_solution == 'barenblatt') :
    solution=calculs_2D.compute_solution(Xgrille,Ygrille, t) 
elif (config.name_solution == 'diffusion') :    
    solution = calculs_2D.compute_solution(Xgrille, Ygrille, t  - config.Tini)
    graphes2D.draw_analytical_solution(Xgrille,Ygrille,solution, npic, t, label=name_analytical)

if (config.method=='SP') :
    R_sp = reconstruction_densites_2D.rec_densite_grille_unif_sp(X, M, config.epsilon, Xgrille, Ygrille) 
    print "\nepsilon  = " , config.epsilon  
    npicf_sp = 2000 + npic
    graphes2D.fait_series_dessins(Xgrille, Ygrille, R_sp, npicf_sp, t, name)
    relative_error=calculs_2D.error(solution , R_sp)
    relative_error_L1=calculs_2D.error_Lp(Xgrille, Ygrille, solution , R_sp, 1)
    relative_error_L2=calculs_2D.error_Lp(Xgrille, Ygrille, solution , R_sp, 2)
    numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_R_sp_mat.txt', R_sp)

if ((config.method=='LTP') or (config.method == 'analytic_debug')):
    R_ltp = reconstruction_densites_2D.rec_densite_grille_unif_LTP(X, M, D, Xgrille, Ygrille)
    npicf_ltp = 3000 + npic
    graphes2D.fait_series_dessins(Xgrille, Ygrille, R_ltp, npicf_ltp, t, name,show=False)
    relative_error=calculs_2D.error(solution , R_ltp)
    relative_error_L1=calculs_2D.error_Lp(Xgrille, Ygrille, solution , R_ltp, 1)
    relative_error_L2=calculs_2D.error_Lp(Xgrille, Ygrille, solution , R_ltp, 2)
#    print(path_to_output_file_name)
    numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_R_ltp_mat.txt', R_ltp)        
    numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D0.txt', D[0,:])
    numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D1.txt', D[1,:])
    numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D2.txt', D[2,:])
    numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D3.txt', D[3,:])
    
if (config.method=='LTP_casier'):
    R_ltp = reconstruction_densites_2D.rec_densite_grille_unif_LTP(X, M, D, Xgrille, Ygrille)
    npicf_ltp = 3000 + npic
    graphes2D.fait_series_dessins(Xgrille, Ygrille, R_ltp, npicf_ltp, t, name,show=False)
    relative_error=calculs_2D.error(solution , R_ltp)
    relative_error_L1=calculs_2D.error_Lp(Xgrille, Ygrille, solution , R_ltp, 1)
    relative_error_L2=calculs_2D.error_Lp(Xgrille, Ygrille, solution , R_ltp, 2)
#    print(path_to_output_file_name)
    numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_R_ltp_mat.txt', R_ltp)        
    numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D0.txt', D[0,:])
    numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D1.txt', D[1,:])
    numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D2.txt', D[2,:])
    numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D3.txt', D[3,:])
    
    #[Xgrille, Ygrille] = make_grid_unif(Ix, Iy, Nmesh_visu3D)
    #R_ltp = rec_densite_grille_unif_LTP(X, M, D, Xgrille, Ygrille)
    #npicf_ltp = 5000 + npic
    #fait_series_dessins(Xgrille, Ygrille, R_ltp, npicf_ltp, t, nom)
    #dessin3d(Xgrille, Ygrille, R_ltp, npicf_ltp, t, label= nom)
    
if (config.draw_sol == 'oui' ) :
    numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_solution_mat.txt', solution)
    graphes2D.draw_analytical_solution(Xgrille,Ygrille,solution, npic, t, label=config.name_solution)

numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_Xgrille.txt', Xgrille)
numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_Ygrille.txt', Ygrille)
numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_M.txt', M)
numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_X0.txt', X[0,:])
numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_X1.txt', X[1,:])
numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_data.txt', numpy.array([t, npic]))
numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D0.txt', D[0,:])
numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D1.txt', D[1,:])
numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D2.txt', D[2,:])
numpy.savetxt(str(path_to_output_file_name)+'_t='+str(t)+'_D3.txt', D[3,:])
	
print "\nnumber of iterations : " , Nboucle
print "\n\nrelative L-infini error with N = " , N, ", at t =  " , t , 'is : '  , relative_error , "\n"
print "relative L1 error with N = " , N, ", at t =  " , t , 'is : '  , relative_error_L1 , "\n"
print "relative L2 error with N = " , N, ", at t =  " , t , 'is : '  , relative_error_L2 , "\n"

if ( not(os.path.isfile(str(path_to_output_file_name))) ) : #if the file does not exist, we create it
    os.system("touch "+str(path_to_output_file_name))
    myfile=open(str(path_to_output_file_name),'a')
    if config.method == 'analytic_debug' :        
        myfile.write("input_file_name\tLinfinity_error\tL1_error\tL2_error\tLinfinity_error_X\tLinfinity_error_D\n")
    else :
        myfile.write("input_file_name\tLinfinity_error\tL1_error\tL2_error\n")
else :
    myfile=open(str(path_to_output_file_name),'a')

if (config.method == 'analytic_debug') :
    myfile.write(str(path_to_input_file_name)+"\t"+str(float(relative_error))+"\t"+str(float(relative_error_L1))+"\t"+str(float(relative_error_L2))+"\t"+str(error_particles_positions)+"\t"+str(error_particles_deformations)+"\n")
else :
    myfile.write(str(path_to_input_file_name)+"\t"+str(float(relative_error))+"\t"+str(float(relative_error_L1))+"\t"+str(float(relative_error_L2))+"\n")

myfile.close()
    
print "\n\nWe're done!"
print 'simulation time: ',  time.time() - start_code_bis



# ===== A FAIRE =======
# faire calcul de U et D par f2py ?
# calcul de D par Differences finies
# fusion et dedoublement de particules

# ==== A DEBUGGER =====
# remapping f2py
# calculsfor_rec_modf90.rec_sp	
# initializepartavecbords            
