# derniere modif le 09-10-15 par Frederique

from __future__ import division
import numpy  #as np
import scipy  #as sp

import math  # ou : from math import *
import time
import sys
import csv

import config
import barenblatt_fcts
import diffusion_fcts
import shapefunction2D

from _calculsfor_f90 import  calculsfor_rec_modf90


"""
#______________________________________________________________________#
#
#			 EQUATION SOLVED: BARENBLATT EQUATION
# $\partial _t f - \nabla \cdot (mf^{m-1} \nabla f) = 0$
# m is a parameter (we mainly consider the case m=2)
#______________________________________________________________________#
"""


"""
#______________________________________________________________________#

        /!\ /!\ /!\ /!\ /!\ /!\        
        /!\ /!\ WARNING /!\ /!\
        /!\ /!\ /!\ /!\ /!\ /!\
        
	-Comments with (??) are things that must be CHECKED
    -Comments with (!!) are things that must be CHANGED
    -We DID NOT do the case redemarrage == 1 in this new version of data (2016/12/01)

#______________________________________________________________________#
"""

def initialize_data_parameters(file_name) : 	
	# we loop through file named 'file_name' and set some config variables  
	# after the loop through the file, the other config variables are set 
	# accordingly
    config.redemarrage = 0
    # ------------- START READING FILE -------------------
    with open(file_name, 'rb') as f:
        reader = csv.reader(f, delimiter='=')
        for row in reader:
            #~ print row
			#  ---- Set initial time, time step and final time  ----
			# !!! WARNING !!! Tmax must be such as (Tmax - Tini) is a modulo of dt 
			# otherwise set : Tmax = Tini+dt*round((Tmax-Tini)/dt,0) ?
            if (row[0] == 'Tini') : 	config.Tini = float(row[1])
            if (row[0] == 'Tmax') : 	config.Tmax = float(row[1])			
            if (row[0] == 'dt') : 		config.dt = float(row[1])
			
			# ----- Time discretisation -----
			# time_scheme = 'euler_explicit'
			# time_scheme = 'middle_point'
            if (row[0] == 'time_scheme') : 		config.time_scheme = str(row[1])
			
			#  ---- Distance to origin of initial domain borders  ----
            if (row[0] == 'Lx1') : 		config.Lx1 = float(row[1])
            if (row[0] == 'Lx2') : 		config.Lx2 = float(row[1])
            if (row[0] == 'Ly1') : 		config.Ly1 = float(row[1])
            if (row[0] == 'Ly2') : 		config.Ly2 = float(row[1])
									
            if (row[0] == 'problem') : config.problem = str(row[1])
            if (row[0] == 'm_barenblatt') : config.m_barenblatt = float(row[1])
            if (row[0] == 'd') : config.d = float(row[1])
            if (row[0] == 'C') : config.C = float(row[1])						            					

			##-----  Numerical method choice ----- 
			# method = 'part_Prand' : particulaire avec positions initiales aleatoires
			# method = 'part_Mrand' : particulaire avec masse aleatoire
			# method = 'part' : particulaire avec discretisation de rhoini
			# method = 'part_avec_suivi_vol' : on ne recontruit pas les densites mais on deforme les particules
			# method = 'SP' : smooth particule, particules de taille epsilon
			# method = 'LTP' : deformation des particules
			# methode = 'PIC' : SP avec reconstruction PIC champ d'advection
			# methode = 'LTPIC' : LTP avec reconstruction PIC champ d'advection
			#LTP_convol not implemented yet 
			# method = 'LTP_convol' : LTP with velocity field U evaluated as  G_\epsilon * \rho_{LTP} D_method can still be either hybrid, explicit or implicit
									 # with G_epsilon --> grad(dirac)
			#method =  'SP' # 'LTP' #			
			#~ method =  'SP'
			# method = 'analytic_debug'
            if (row[0] == 'method') : 
				config.method = str(row[1])
				if config.method =='analytic_debug' : config.D_method = None
			
			# ----- Choice of method to compute the matrices Dk ----- 
			# D_method = 'explicit' # one updates particle positions Xk AFTER  Dk matrices (by "explicit" expression of Jacobian of flow)
			# D_method = 'implicit' # one updates particle positions Xk BEFORE Dk matrices. (by Finite Differences)
			# D_method = 'none' # in case of SP method
			# D_method = 'hybrid' # in case of LTP method where we compute explicitly matrices Dkn 
			# D_method = 'convol' # compute D_k^{n+1} as D_k^{n+1} = D_k^{n} * (J_k^n)^{-1}
								  # with (J_k^n)^{-1} = exp(-\delta t * J_U)
								  # where J_U is jacobian of flow U and computed as J_U = H_\epsilon * \rho_{LTP}
								  # with H_\epsilon -> Hess(dirac)
			#but with a larger shape in approximation of jacobian of flow (replace hx ,hy with a epsilon > h = max(hx,hy))

            if (row[0] == 'D_method') : 
				config.D_method = str(row[1])			
				# D_method should be None but just in case there is still 
				# D_method in csv file we set it to None if we are in SP method
				if (config.method) == 'SP' : config.D_method = None

			# indice remapping :
			# 'yes' = we remapp particles based on a criterion
			# 'no' = no remapping
            if (row[0] == 'indic_remapping') : config.indic_remapping = str(row[1])
            if (row[0] == 'l_remap') : config.l_remap = float(row[1]) 				#value set before : l_remap = 6.
            if (row[0] == 'radius_remap') : config.radius_remap = float(row[1]) 	#value set before : radius_remap = 2.


			#----- Method choice to compute weights wk ----- 
			# met_calcul_poids = 'ponctuelle'
			# met_calcul_poids = 'quasi_interp' (methode Martin)
			# met_calcul_poids = 'quadrature'
            if (row[0] == 'met_calcul_poids') : config.met_calcul_poids = str(row[1]) 	#value set before met_calcul_poids = 'quasi_interp' 

			# set particles on border of domain (oui ou non)		
            if (row[0] == 'part_sur_bords') : config.part_sur_bords=str(row[1]).strip()


			#----- Shape function ----- 
			# spline of order != 3 are not implemented...
			# exponential shape functions are imlemented BUT (!!) : 
			# we did NOT changed the initialization of the weights accordingly ... -_______-
			### TO BE MODIFIED : REALLY SHITTY (!!)
            if (row[0] == 'bool_phi_radial') : config.bool_phi_radial = int(row[1]) # value set before bool_phi_radial = 1
            if (row[0] == 'deg_spline') : config.deg_spline = int(row[1]) # value set before deg_spline = 3
    
			#  ---- Number of particules in x and y directions ----			
            if (row[0] == 'Nx') : 		config.Nx = int(row[1])
            if (row[0] == 'Ny') : 		config.Ny = int(row[1])
			
			# -------  set shape of particle in case of LTP-hybrid or SP method  -------
            if (row[0] == 'epsilon') : 		config.epsilon = float(row[1])
    
			# ----- Explosion of particles position moving -----
			# explosion = "yes" or "no"# study explosion or not
            if (row[0] == 'explosion') : 
                config.explosion = str(row[1])
				### criterion to know if particles epxloded in X and/or Y direction
				### basically : if particles moved more than dx in x one time step then it exploded (resp Y)
                config.dx = config.hx_remap
                config.dy = config.hy_remap

    # ------------- END READING FILE -------------------

	#  ---- Interval of initial domain in x and y directions  ----
    
    config.Sx = numpy.array([config.Lx1, config.Lx2])
    config.Sy = numpy.array([config.Ly1, config.Ly2])
	
	#  ---- Grid for vizualization----
    config.xgmin = config.Sx[0]-abs(config.Sx[1]-config.Sx[0])/5
    config.xgmax = config.Sx[1]+abs(config.Sx[1]-config.Sx[0])/5
    config.ygmin = config.Sy[0]-abs(config.Sy[1]-config.Sy[0])/5
    config.ygmax = config.Sy[1]+abs(config.Sy[1]-config.Sy[0])/5
    config.Ix = numpy.array([config.xgmin, config.xgmax])
    config.Iy = numpy.array([config.ygmin, config.ygmax])    
    config.Lx = config.xgmax - config.xgmin
    config.Ly = config.ygmax - config.ygmin 
    
    #  ---- Barenblatt variables ----
    config.alpha = config.d/(config.d*(config.m_barenblatt-1)+2)    
    config.beta = config.alpha / config.d
    config.k = config.alpha * (config.m_barenblatt-1) / (2*config.m_barenblatt*config.d)
    
    #  ---- Initial density (non adimensionned)----
    config.initialisation = 'deterministe' #'aleatoire'
    
    if (config.problem == 'barenblatt') : 
        config.rhoini = barenblatt_fcts.barenblatt_init 
        config.name_solution = 'barenblatt' # 'barenblatt' is the name of the function defined above
    elif (config.problem == 'diffusion') : 
        config.rhoini = diffusion_fcts.diffusion_init         
        config.name_solution = 'diffusion'
    elif (config.problem == 'convection_diffusion') : 
        raise ValueError("problem = " , config.problem, "not implemented yet \n")        
        
        
    # -----  Particles list that one wants to trace ----- 
    config.Part_suivies = []

    if (config.bool_phi_radial == 1) : # phi is radial ( so far <=> phi is a spline
        calculsfor_rec_modf90.set_phi_radial(config.bool_phi_radial)    
        calculsfor_rec_modf90.set_degre_spline(config.deg_spline)        
        if (config.deg_spline==1) :
            config.phi = shapefunction2D.b1 
            config.radius_phi = 1 # 2 # rayon du support de phi
            #~ config.radius_phi=1.
        elif (config.deg_spline == 3) :
            config.phi = shapefunction2D.b3 
            config.derivative_phi = shapefunction2D.derivative_b3
            config.phi_2d = shapefunction2D.b3spline2d
            config.phi_2d_h = shapefunction2D.b3spline2d_h
            config.grad_phi = shapefunction2D.grad_b3spline2d
            config.grad_phi_h = shapefunction2D.grad_b3spline2d_h
            config.radius_phi = 2 # 2 # rayon du support de phi
        else :
            print "ERROR degre_spline value must be either 1 or 3. Current degre_spline = " , deg_spline

    elif (config.bool_phi_radial == 0) : # not a radial function
        # phi(x,y) = exp(-1/(1-x**2-y**2)) if (x**2 + y**2) < 1 , 
        # phi(x,y) = 0 otherwise 
        calculsfor_rec_modf90.set_phi_radial(config.bool_phi_radial)
        calculsfor_rec_modf90.set_degre_spline(config.deg_spline)
        
        config.phi = shapefunction2D.exp_compact #radius_phi is declared when we call that function
        config.radius_phi = 1.
    elif (config.bool_phi_radial == 2) : # not a radial function
        # phi(x,y) = exp(-x**2-y**2)/2) / 2 pi 
        calculsfor_rec_modf90.set_phi_radial(config.bool_phi_radial)
        calculsfor_rec_modf90.set_degre_spline(config.deg_spline)        
        
        config.phi = shapefunction2D.gaussian
        config.radius_phi = 3.

    # total number of particles
    config.Nini = config.Nx * config.Ny  
    
    # initial distance between particules :			
    # /max(1,Nx-1) if we set particles at border
    # /(Nx+1) if we do not set particles at border
    if  config.method == 'part_avec_suivi_vol':
        config.hx = 0.06
        config.hy = 0.06
    else:
        config.hx = config.Lx / max(1, config.Nx - 1)  # pour visu mettre 0.06 avec N=4
        config.hy = config.Ly / max(1, config.Ny - 1)
        
    if( config.part_sur_bords == 'non') :
        config.hx_remap = (config.Sx[1] - config.Sx[0]) / (config.Nx)
        config.hy_remap = (config.Sy[1] - config.Sy[0]) / (config.Ny)        
    else :         
        config.hx_remap = (config.Sx[1] - config.Sx[0]) / (config.Nx-1)
        config.hy_remap = (config.Sy[1] - config.Sy[0]) / (config.Ny-1)		
    
    print "config.hx_remap , config.hy_remap = " , config.hx_remap , config.hy_remap
    ### Define epsilon (shape function "size") 
    #~ h = max((Sx[1] - Sx[0]) / float(Nx-1) , (Sy[1] - Sy[0]) / float(Ny-1))
    config.h = max(config.hx_remap,config.hy_remap)

    ### D_method_scheme should disappear in new version (!!)
    if (config.time_scheme == 'euler_explicit') :
        config.D_method_scheme = 'euler_explicit'			
    if (config.time_scheme == 'RK2') :
        config.D_method_scheme = 'RK2'
    if (config.time_scheme == 'middle_point') :
        config.D_method_scheme = 'middle_point'


    #  ---- Reconstruction parameter at final time ------
    config.Nmesh = 50  # parameter of reconstructed function (500 initially)
    config.Nmesh_visu3D = 50
    config.phi_grid = shapefunction2D.b3 # pour la quasi interpolation (utilise ? a verifier) (??)
    config.radius_phi_grid = 2   #(??)
    #  m=7 #nb de points pour la quadrature de W*varphi
    #  l=3*N+1 #number of points for the numerical quadrature


    #  ======================= output data parameters   ======================================
    config.animgif = 0
    config.animmp4 = 0
    config.position_particles = 1
    #~ exportdata = 1
    #~ config.param = -1	   # we plot every time step
    config.param = 100000
    config.pictini = 0
    config.pictfin = 0
    #~ config.trace_particules = 'oui'
    config.trace_particules = 'non'
    #~ config.draw_sol = 'oui'        
    config.draw_sol = 'non'        
    #~ draw_sol = 'non'   
    #~ config.save_data = 'oui'
    config.save_data = 'non'
