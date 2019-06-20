# ###########################################################################################
# Code LTP 2D
# dernieres modifs 09-10-15 par Frederique
#    - ajout routines pour sous-traitance f2py de reconstruction, initialisation, remapping
# Resolution de \partial_t rho + div(rho u) = 0 avec u donne

#  !/usr/bin/env python
#  -*- codin: utf-8 -*-
from __future__ import division
import numpy  #as np
import matplotlib.pyplot as plt
import time


from shapefunction2D import *
from data_code_LTP_2D import *
from initialisation_2D import *
from calculs_2D import *
from graphes2D import *
from reconstruction_densites_2D import *
from remapping_2D import *

# for test, move later
from _calculsfor_f90 import  calculsfor_rec_modf90
from _calculsfor_f90 import  calculsfor_var_modf90

start_code = time.clock()
start_code_bis = time.time()

#=====  Initialisation ====== 
t = Tini
npic_part = 0


	
#----- initialisation des positions et poids -----

print 'Initialisation'
if redemarrage == 0:
    #X, M = initialize_avec_f2py(Sx, Sy, Nx, Ny)
    X, M = initialize_avec_python(Sx, Sy, Nx, Ny)
    #print 'suppression des particules vides'
    #N, X, M = supprime_particule(X,M)
    N=Nini
    if N==Nini:
        print 'pas de suppression de particules'
    npic = 0
else:
    X0 = numpy.genfromtxt("out/X0.txt",delimiter=",") 
    X1 = numpy.genfromtxt("out/X1.txt",delimiter=",") 
    M = numpy.genfromtxt("out/M.txt",delimiter=",") 
    N = len(M)
    X = numpy.zeros([2, N])
    t, npic = numpy.genfromtxt("out/data.txt",delimiter=",") 
    print 'redemarrage a t= ', t
    X[0, :] = X0
    X[1, :] = X1

print 'Nombre de particules N=', N 
   
U = numpy.zeros([2, N])
Mini = M.copy()
numpy.savetxt('out/Mini.txt', Mini)
Xinif = X.copy()
print 'somme des poids', sum(Mini)
Nboucle = 0

#----- Grille de reconstruction -----
[Xgrille, Ygrille] = make_grid_unif(Ix, Iy, Nmesh)


# ----- Initialisation des supports -----

D = numpy.zeros([4, N]) #matrice de deformation
D[0,:] = 1.  # coeff (1,1) de Dk^0
D[3,:] = 1.  # coeff (2,2) de Dk^0
# D[0,:] -> coef (1,1) de Dkn
# D[1,:] -> coef (1,2) de Dkn
# D[2,:] -> coef (2,1) de Dkn
# D[3,:] -> coef (2,2) de Dkn
#  Dk = D[:,k].reshape((2, 2)) -> matrice Dkn
#for k in range(N) :
	#print k , ' , ' , X[0,k] , ' , ' , X[1,k]
### Moved to data_code_LTP_2D
#hx_remap = (Sx[1] - Sx[0]) / float(Nx)  #  distance between particles
#hy_remap = (Sy[1] - Sy[0]) / float(Ny)


print "hx_remap , hy_remap : " , hx_remap , hy_remap
start1=time.clock()

if D_method == 'implicit' :	
	#tab=tab_cardinal_particles(X,N,hx_remap,hy_remap)
	tab = calculsfor_var_modf90.tab_cardinal_particles(X,hx_remap,hx_remap,N)	
	for k in range(0,N) :
	#	print "k , tab[k][:] : " , k , tab[k]
		print "k , X[k][:] : " , k , X[:,k]
	print '\n\ntab_cardinal_particles : ' , tab  , 	'\n\n' 
	print 'temps de calcul du tableau des positions des particules', time.clock() - start1 	,'\n'

	#my_list=boundary_particles(tab)
	#print "\n\nlist of boundary particles : " , my_list ,'\n'

#======  Boucle en temps ===== 
print 'Debut de la boucle en temps'
print "time taken to complete all initialisations : " , time.time() - start_code_bis
#while (round(t,13) < (round((Tmax -dt),13)) ) :
while (round(t,13) < (round((Tmax),13)) ) :    

	#copy X to compare the explosions in the particles moving (see end of the loop)	        
    X_old = X.copy()
    #D_old = D.copy()
    if (method == 'PIC') or (method == 'LTPIC'):
        U = calcule_U_PIC(X, M, D, Xgrille, Ygrille)
	
    # INUTILE PUISQUE RECALCULE AU MOMENT DE BOUGER LES PARTICULES DANS LA FONCTION move_particles
    #else:
    #    U = calcule_U_part(X, N, t)

    # ----- Trace densites reconstruites et enregistrement -----
    if (Nboucle % param == 0) :
		
        if method == 'LTP':
            start2 = time.clock()
            time1= time.time()
            R_ltp = rec_densite_grille_unif_LTP(X, M, D, Xgrille, Ygrille)            
            #for i in range(0,R_ltp.shape[0]) :
			#	for j in range(0,R_ltp.shape[1]) :
			#		print "R_ltp[i,j] = " , R_ltp[i,j]
            #print 'temps de reconstruction de rho avec LTP', time.clock() - start2 
            #for k in range(0,N) : 
				#R_ltp_part=calculsfor_rec_modf90.rec_ltp_at_coordinate(X, M, D, hx_remap, hy_remap, X[:,k], N)
				#sol = barenblatt(X[0,k],X[1,k],t)
				#print  " k , sol , R_ltp_part " , k , sol , R_ltp_part 
            print 'temps de reconstruction de rho avec LTP ' , time.time() - time1
            npicf_ltp =  npic
            fait_series_dessins(Xgrille, Ygrille, R_ltp, npicf_ltp, t, 'LTP')    
            #if Nboucle == 0 :    
				#solution=compute_solution(Xgrille,Ygrille, t) #= initial profile because reverse swirling velocity field respects u(0,x) = u(Tmax,x)
				#print "\n\nerror(solution, R_ltp) at t =  " , t , 'is : ' , error(solution ,R_ltp)
                        
            #indic_y=int(len(Ygrille)/2)
            #rep1d_plot(Xgrille, indic_y, R_ltp, npic, t, label='Density_LTP_1D')    	

        if (method =='part') or (trace_particules == 'oui'):
            traceparticules(X, t, npic)

        if (method =='SP'):
            start3 = time.clock()
            R_sp = rec_densite_grille_unif_sp(X, M, epsilon, Xgrille, Ygrille)
            print 'temps de reconstruction de rho avec SP', time.clock() - start3 
            npicf_sp = npic
            fait_series_dessins(Xgrille, Ygrille, R_sp, npicf_sp, t, 'SP')      
                  
        if (draw_sol == 'oui' ) :
            solution=compute_solution(Xgrille,Ygrille, t)            
            draw_analytical_solution(Xgrille,Ygrille,solution, npic, t, label=name_solution)
            indic_y=int(len(Ygrille)/2)
            rep1d_plot(Xgrille, indic_y, solution, npic, t,name_solution+'_1D')
            relative_error=error(solution , R_ltp)
            print "\n\nrelative L-infini error with N = " , N, ", at t =  " , t , 'is : '  , relative_error
            
        #------------------------------       
        print 'Nplot = ', npic

        npic += 1    
        plt.close("all")
		
    print 't=', t
    
    if method == 'SP' :		
		k=318
		grad_f = calculsfor_rec_modf90.rec_sp_grad_at_coordinate(X, M, epsilon, X[:,k], N)
		grad_sol = grad_barenblatt(X[0,k], X[1,k],t)
		print " k , grad_sol , grad_f : " , k , grad_sol , grad_f
		time_part=time.time()
		#X=move_particles(X,N,t,dt) # update X with time_scheme chosen in data_code_LTP_2D.py
		diff_field = calculsfor_rec_modf90.diffusion_field_barenblatt_sp(X, M, m_barenblatt,epsilon, N)
		X[0,:] += dt * (-diff_field[0,:])
		X[1,:] += dt * (-diff_field[1,:])
		print "time taken to move the particles X : " , time.time() - time_part								
		Nboucle += 1
		t += dt		
    
    #__________________________________________________________________#
    #
    ### MOVE PARTICLES (CASE D by IMPLICIT method = finite differences)
    #__________________________________________________________________#
    if D_method == 'implicit' :        
		### barenblatt with gradient reconstructed with SP
		time_part=time.time()
		#diff_field = calculsfor_rec_modf90.diffusion_field_barenblatt_sp(X, M, m_barenblatt,epsilon, N)
		#X[0,:] += dt * (-diff_field[0,:])
		#X[1,:] += dt * (-diff_field[1,:])
				
		### barenblatt
		X=move_particles_diffusion_barenblatt(X, N , M , D,m_barenblatt,hx_remap,hy_remap) 
		#for k in range(0,N) : 
			#print "k, X[0,k], X[1,k] : " , k, X[0,k] , X[1,k]
			
		### advection
        #X=move_particles(X,N,t,dt) # update X with time_scheme chosen in data_code_LTP_2D.py
		print "time taken to move the particles X : " , time.time() - time_part
		Nboucle += 1
		t += dt
		
	#__________________________________________________________________#
    #
    # MOVE PARTICLES (CASE D by EXPLICIT method Dkn+1 = Dkn * inv(Jtn,tn+1)
    #				  and Jk computed by euler implicit scheme)
    #__________________________________________________________________#
    if ((D_method == 'explicit') & (D_method_scheme == 'euler_implicit')) :
		time_part=time.time()
		### advection case
        #X=move_particles(X,N,t,dt) # update X with time_scheme chosen in data_code_LTP_2D.py
		
		###barenblatt diffusion reconstruction
		#X=move_particles_diffusion_barenblatt(X, N , M , D,m_barenblatt,hx_remap,hy_remap)        
		
		#diff_field = calculsfor_rec_modf90.diffusion_field_barenblatt_sp(X, M, m_barenblatt,epsilon, N)
		#X[0,:] += dt * (-diff_field[0,:])
		#X[1,:] += dt * (-diff_field[1,:])
		
		### exact barenblatt for test
		#for i in range(0,N) :
			#grad_sol = grad_barenblatt(X[0,i], X[1,i],t)
			#X[0,i] += -2*dt*grad_sol[0]	#case m=2
			#X[1,i] += -2*dt*grad_sol[1] #case m=2
		print "time taken to move the particles X : " , time.time() - time_part
		Nboucle += 1
		t += dt
		k=318
		grad_f = calculsfor_rec_modf90.rec_ltp_grad_at_coordinate(X, M, D, hx_remap, hy_remap, X[:,k], N)				
		grad_sol = grad_barenblatt(X[0,k], X[1,k],t)
		print "X[:,318] : " , X[:,318]
		print " k , det(D[:,k]), grad_sol, grad_f : " , k , calculdet(D[:,k].reshape((2, 2))), grad_sol, grad_f
	#__________________________________________________________________#
    #
	### UPDATE D
    #__________________________________________________________________#
    #----- Mise a jour matrice de deformation -----
    if (method == 'LTP') or (method == 'part_avec_suivi_vol'): 
        start4 = time.time()
        if D_method == 'implicit' : 
			D, Norm_inf_Dm1, indice_max_norm_Dm1 = MaJ_D_FD_scheme(X, D ,tab, t, dt , hx_remap ,hy_remap) 
		
        if D_method == 'explicit' : 
			#D, Norm_inf_Dm1, indice_max_norm_Dm1 = MaJ_D(D, X, t, dt)	
			D_old = D.copy()		
			D, Norm_inf_Dm1, indice_max_norm_Dm1 = MaJ_D_barenblatt(D, X, M, t, dt, hx_remap, hy_remap)			
			#D, Norm_inf_Dm1, indice_max_norm_Dm1 = MaJ_D_barenblatt_analytic(D, X, M, t, dt)
							        
        print 'temps de calcul de D', time.time() - start4
        print 'max des normes des Dk', Norm_inf_Dm1
        print 'indice_max_norm_Dm1', indice_max_norm_Dm1
         
         
    #__________________________________________________________________#
    #
    # MOVE PARTICLES (CASE D by EXPLICIT method Dkn+1 = Dkn * inv(Jtn,tn+1))
    #__________________________________________________________________#
    if ((D_method == 'explicit') & (D_method_scheme != 'euler_implicit')) :
		time_part=time.time()
		### advection case
        #X=move_particles(X,N,t,dt) # update X with time_scheme chosen in data_code_LTP_2D.py
		
		###barenblatt diffusion reconstruction
		X=move_particles_diffusion_barenblatt(X, N , M , D_old,m_barenblatt,hx_remap,hy_remap)        
		
		### exact barenblatt for test
		#for i in range(0,N) :
			#grad_sol = grad_barenblatt(X[0,i], X[1,i],t)
			#X[0,i] += -2*dt*grad_sol[0]	#case m=2
			#X[1,i] += -2*dt*grad_sol[1] #case m=2
		print "time taken to move the particles X : " , time.time() - time_part
		Nboucle += 1
		t += dt
		k=318
		grad_f = calculsfor_rec_modf90.rec_ltp_grad_at_coordinate(X, M, D, hx_remap, hy_remap, X[:,k], N)				
		grad_sol = grad_barenblatt(X[0,k], X[1,k],t)
		print "X[:,318] : " , X[:,318]
		print " k , det(D[:,k]), grad_sol, grad_f : " , k , calculdet(D[:,k].reshape((2, 2))), grad_sol, grad_f
		         
    #---- trace forme particule -----
    if  (Nboucle % param == 0) and (len(Part_suivies)>0) and ((method=='LTP') or (method == 'part_avec_suivi_vol')):
         #  reconstruction de la forme des particules suivis
         #dessine_forme_particule(D, M, X, h, Part_suivies, npic_part)
         for k in Part_suivies:
             Xgk, Ygk, R_k = rec_forme_particule(k, X[:,k], D[:, k], M[k], 100)
             rep2d_contourf(Xgk, Ygk, R_k, npic, t, 'Shape_of_particle_'+str(k))
                #plt.close("all")
         npic_part +=1
	

    
    #---- Mise a jour eventuelle de la grille -----
    Ix, Iy, modif_grille = MAJ_grille(X,Ix,Iy)    
    if modif_grille == 1:
        [Xgrille, Ygrille] = make_grid_unif(Ix, Iy, Nmesh)

    # ---- remapping -----------------------
    if (method == 'LTP'):
        if (Norm_inf_Dm1 >radius_remap) and (indic_remapping == 'oui'):
			print '\nRemapping!\n'
			if D_method == 'implicit' :
				X, M, D  = remapping_avec_python_test(hx_remap, hy_remap, X, M, D, t) 
				N = len(M)   
				tab = calculsfor_var_modf90.tab_cardinal_particles(X,hx_remap,hx_remap,N)	
			if D_method =='explicit' :
				X, M, D  = remapping_avec_python_test(hx_remap, hy_remap, X, M, D, t) # remapping that also returns remapping grid step				
				N = len(M)   
			print 'nouveau nombre de particules', N
	
	#X_new_test=np.zeros([2,N])
	#print "\ngradient of non-exploding particles : \n" 
	#for k in range(0,N) :
		#grad_f = calculsfor_rec_modf90.rec_ltp_grad_at_coordinate(X_old, M, D, hx_remap, hy_remap, X_old[:,k], N)
		#X_new_test[:,k] = X_old[:,k] - 2* dt * grad_f[:]				
		#grad_sol = grad_barenblatt(X_old[0,k], X_old[1,k],t-dt)
		#print k , grad_sol, grad_f
		
	#______________________________________________________________________#
	#
	#		PARTICLES MOVING EXPLOSION TEST
	#______________________________________________________________________#
	### We test if particles moving explode (related to dt)
    if (explosion == "yes") :
        X_diff = X_old - X
        explo_list=[]
        explo = False
        for k in range(0,N) :
            if ((abs(X_diff[0,k]) >= dx) | (abs(X_diff[1,k]) >= dy)) :
                explo_list.append([k, X_old[0,k], X_old[1,k], X[0,k], X[1,k], X_diff[0,k], X_diff[1,k], dx, dy])
                explo = True
        if (explo) :
            myfile=open('./out/'+str(outfile_name),'a')
            myfile.write('explosion of particles moving at t = '+str(t)+'\n')
            myfile.close()			
            N_explo = int(len(explo_list))
            X_explo_old = np.zeros([2,N_explo])
            X_explo = np.zeros([2,N_explo])
            Xini_explo = np.zeros([2,N_explo])
            for i in range(0,N_explo) : 
                X_explo_old[0,i] = explo_list[i][1]
                X_explo_old[1,i] = explo_list[i][2]
                X_explo[0,i] = explo_list[i][3]
                X_explo[1,i] = explo_list[i][4]
                Xini_explo[0,i] = Xinif[0,explo_list[i][0]]
                Xini_explo[1,i] = Xinif[1,explo_list[i][0]]
			### draw position of particles who will explode at INTIAL TIME
			#traceparticules(Xini_explo, Tini, 10000-1+npic,show=True)
			### draw position of particles BEFORE they explode
			#traceparticules(X_explo_old, t -dt, 10000+npic,show=True)					
			### draw position of particles AFTER they explode
			#traceparticules(X_explo, t, 10000+npic+1,show=True)
            R_ltp_vec=np.zeros(len(explo_list))
            R_ltp_old_vec=np.zeros(len(explo_list))
            fexplo=open('./out/explosion.txt','w')
            fexplo.write('k, X_old[0,k], X_old[1,k], X[0,k], X[1,k], X_diff[0,k], X_diff[1,k], dx, dy, R_ltp_old[k], R_ltp[k]\n')
            for i in range(0,int(len(explo_list))) :
                k = explo_list[i][0]				
                #R_ltp_old_vec[i] = calculsfor_rec_modf90.rec_ltp_at_coordinate(X_old, M, D_old, hx_remap, hy_remap, X_old[:,k], int(len(X_old[0,:])))
                #R_ltp_vec[i] = calculsfor_rec_modf90.rec_ltp_at_coordinate(X, M, D, hx_remap, hy_remap, X[:,k], int(len(X[0,:])))
                explo_list[i].append(R_ltp_old_vec[i])
                explo_list[i].append(R_ltp_vec[i])
                fexplo.write(str(explo_list[i])[1:-1]+'\n')
                print k , Xinif[:,k] , barenblatt(Xinif[0,k],Xinif[1,k],Tini)
            fexplo.close()
            if (method=='LTP'):
                R_ltp_old = rec_densite_grille_unif_LTP(X_old, M, D_old, Xgrille, Ygrille)	
                fait_series_dessins(Xgrille, Ygrille, R_ltp_old, npicf_ltp+10000, t-dt, 'LTP', show=True)
                R_ltp = rec_densite_grille_unif_LTP(X, M, D, Xgrille, Ygrille)	
                fait_series_dessins(Xgrille, Ygrille, R_ltp, npicf_ltp+1+10000, t, 'LTP', show=True)                        
            if (method=='SP'):
                R_sp_old = rec_densite_grille_unif_sp(X_old, M, epsilon, Xgrille, Ygrille)				
                R_sp = rec_densite_grille_unif_sp(X, M, epsilon, Xgrille, Ygrille)				
                npicf_sp = npic				
                fait_series_dessins(Xgrille, Ygrille, R_sp_old, 30000 + npicf_sp , t-dt, 'SP')
                fait_series_dessins(Xgrille, Ygrille, R_sp, 30000 + npicf_sp +1, t, 'SP')
            traceparticules(X_old, t-dt, 20000+npic-1,show=True)
            traceparticules(X, t-dt, 20000+npic,show=True)
			#print "\ngradient and positions of exploding particles : \n" 
			#for i in range(0,int(len(explo_list))) :
				#k = explo_list[i][0]
				#grad_f = calculsfor_rec_modf90.rec_ltp_grad_at_coordinate(X_old, M, D, hx_remap, hy_remap, X_old[:,k], N,1)
				#grad_sol = grad_barenblatt(X_old[0,k], X_old[1,k],t-dt)
				#print "\nk , grad_sol, grad_f : " , k , grad_sol, grad_f , "\n"
			
##			X_new_test=np.zeros([2,N])
			#print "\ngradient of non-exploding particles : \n" 
			#for k in range(0,N) :
				#grad_f = calculsfor_rec_modf90.rec_ltp_grad_at_coordinate(X_old, M, D, hx_remap, hy_remap, X_old[:,k], N)
				##X_new_test[:,k] = X_old[:,k] - 2* dt * grad_f[:]				
				#grad_sol = grad_barenblatt(X_old[0,k], X_old[1,k],t)
				#print k , grad_sol, grad_f
			
			### TEST GRADIENT DENSITY
			#X_test=move_particles_diffusion_barenblatt(X_old, N , M , D ,m_barenblatt,hx_remap,hy_remap)			
			#for i in range(0,int(len(explo_list))) :
				#k = explo_list[i][0]
				#print k, X_old[:,k], X[:,k], X_test[:,k], R_ltp_old, R_ltp
			    
            raise ValueError('Explosion in particles moving see file in ./out/explosion.txt')

    # ---- enregistrement des tableaux -----
    if save_data == 'oui' :
        time_tabs=time.time()
        numpy.savetxt('out/M.txt', M)
        numpy.savetxt('out/X0.txt', X[0,:])
        numpy.savetxt('out/X1.txt', X[1,:])
        numpy.savetxt('out/data.txt', numpy.array([t, npic]))
        #numpy.savetxt('out/Dtranspose.txt', numpy.transpose(D))
        numpy.savetxt('out/D0.txt', D[0,:])
        numpy.savetxt('out/D1.txt', D[1,:])
        numpy.savetxt('out/D2.txt', D[2,:])
        numpy.savetxt('out/D3.txt', D[3,:])
        print "time to save data in txt files : " , time.time()-time_tabs
        
#===== fin de la boucle en temps ======



if (method == 'LTP') or (method == 'part_avec_suivi_vol'): 
     start4 = time.clock()     
     #if D_method == 'explicit' : 
		##D, Norm_inf_Dm1, indice_max_norm_Dm1 = MaJ_D(D, X, t, dt)
		#D, Norm_inf_Dm1, indice_max_norm_Dm1 = MaJ_D_barenblatt(D, X, M, t, dt, hx_remap, hy_remap)
 
#===== desssins finaux =====
if animmp4 == 1:
    nom ='Density_'+method+'_'
    faitlefilm(nom)
       
#traceparticules(X,t,npic)

nom = method+'N='+str(N)+'-Nmesh='+str(Nmesh) 

#solution=sol_final_SW(Xgrille,Ygrille) #= initial profile because reverse swirling velocity field respects u(0,x) = u(Tmax,x)
#draw_analytical_solution(Xgrille,Ygrille,solution, npic, t, label='Solution_SW')

solution=compute_solution(Xgrille,Ygrille, t) 


if method=='SP':
    R_sp = rec_densite_grille_unif_sp(X, M, epsilon, Xgrille, Ygrille) 
    print "\nepsilon  = " , epsilon  
    npicf_sp = 2000 + npic
    fait_series_dessins(Xgrille, Ygrille, R_sp, npicf_sp, t, nom)
    relative_error=error(solution , R_sp)

	
if method=='LTP':
    R_ltp = rec_densite_grille_unif_LTP(X, M, D, Xgrille, Ygrille)
    npicf_ltp = 3000 + npic
    fait_series_dessins(Xgrille, Ygrille, R_ltp, npicf_ltp, t, nom,show=True)
    relative_error=error(solution , R_ltp)
    
    #[Xgrille, Ygrille] = make_grid_unif(Ix, Iy, Nmesh_visu3D)
    #R_ltp = rec_densite_grille_unif_LTP(X, M, D, Xgrille, Ygrille)
    #npicf_ltp = 5000 + npic
    #fait_series_dessins(Xgrille, Ygrille, R_ltp, npicf_ltp, t, nom)
    #dessin3d(Xgrille, Ygrille, R_ltp, npicf_ltp, t, label= nom)
    
if (draw_sol == 'oui' ) :
	#solution=solution_barenblatt(Xgrille,Ygrille, t) #Barenblatt profile
	#solution=sol_final_SW(Xgrille,Ygrille) #= initial profile because reverse swirling velocity field respects u(0,x) = u(Tmax,x)
	
	draw_analytical_solution(Xgrille,Ygrille,solution, npic, t, label=name_solution)
	#len_Ygrille=len(Ygrille)
	
	#indic_y=int(len_Ygrille/2)
	#print "int(len_Ygrille/2) , y = " , indic_y , Ygrille[indic_y]
	#rep1d_plot(Xgrille, indic_y, solution, npic, t, label='Solution_linear_advection_1D')    	
	#rep1d_plot(Xgrille, indic_y, R_ltp, npic, t, label='Density_LTP_1D')    	

print "\nnumber of iterations : " , Nboucle

print "\n\nrelative L-infini error with N = " , N, ", at t =  " , t , 'is : '  , relative_error
myfile=open('./out/'+str(outfile_name),'a')
myfile.write('relative L-infini error with N = '+str(N)+' at t = '+str(t)+' is : '+str(float(relative_error))+'\n')
myfile.close()
print "C'est fini !"
print 'temps de simulation: ',  time.clock() - start_code

print 'temps de simulation bis: ',  time.time() - start_code_bis
# ===== A FAIRE =======
# faire calcul de U et D par f2py ?
# calcul de D par Differences finies
# fusion et dedoublement de particules

# ==== A DEBUGGER =====
# remapping f2py
# calculsfor_rec_modf90.rec_sp	
# initializepartavecbords
