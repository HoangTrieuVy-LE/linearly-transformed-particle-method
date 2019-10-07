# derniere modif le 09-10-15 par Frederique

from __future__ import division
import numpy  #as np
import scipy  #as sp
import pylab  #as plt
import math  # ou : from math import *
import time
from shapefunction2D import *
from initial_densities import *

from _calculsfor_f90 import  calculsfor_rec_modf90

redemarrage = 0

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
#______________________________________________________________________#
"""

### Initial and final times 
# Tini must be strictly greater than 0, 
# since at T=0 the density function is a Dirac and we can't implement it 
Tini=0.05	
Tmax = 0.052


#______________________________________________________________________#
#
# 				THEORETICAL BARENBLATT DENSITIES
#______________________________________________________________________#

## order m of PME equation $\partial _t f - \nabla \cdot (mf^{m-1} \nabla f) = 0$
m_barenblatt = 2
d=2 # spatial dimension 


def barenblatt(x,y,t=Tini) :	 
    """
    Theoretical profile of Barenblatt density
    """   
    alpha = d/(d*(m_barenblatt-1)+2)
    C=1. #arbitrary constant
    beta = alpha / d
    k = alpha * (m_barenblatt-1) / (2*m_barenblatt*d)
    sol_barenblatt = 1/(t**alpha) * (max( C - k *(x**2 + y**2) * (1/t**(2*beta)) , 0))**(1/(m_barenblatt-1))		    
    return sol_barenblatt

def grad_barenblatt(x,y,t=Tini) :	
    """
	Theoretical profile of gradient of Barenblatt density
    """  
    grad_sol_barenblatt=numpy.zeros(2) 
    alpha = d/(d*(m_barenblatt-1)+2)
    C=1. #arbitrary constant
    beta = alpha / d
    k = alpha * (m_barenblatt-1) / (2*m_barenblatt*d)
    if (C - k *(x**2 + y**2) * (1/t**(2*beta)) > 0):
		if (m_barenblatt == 2) :
			grad_sol_barenblatt[0] = 1/(t**alpha) * (-2*k*x*(1/t**(2*beta))) 
			grad_sol_barenblatt[1] = 1/(t**alpha) * (-2*k*y*(1/t**(2*beta))) 
		else :
			grad_sol_barenblatt[0] = 1/(t**alpha) * (-2*k*x*(1/t**(2*beta))) * ( C - k *(x**2 + y**2) * (1/t**(2*beta)))**(1/(m_barenblatt-1) - 1)
			grad_sol_barenblatt[1] = 1/(t**alpha) * (-2*k*y*(1/t**(2*beta))) * ( C - k *(x**2 + y**2) * (1/t**(2*beta)))**(1/(m_barenblatt-1) - 1)
    else : 
        grad_sol_barenblatt[:] = 0.
    return grad_sol_barenblatt

def hess_barenblatt(x,y,t=Tini) :	
    """
	Theoretical profile of Hessian Barenblatt density
    """  
    hess_sol_barenblatt=numpy.zeros([2,2])
    alpha = d/(d*(m_barenblatt-1)+2)
    C=1. #arbitrary constant
    beta = alpha / d
    k = alpha * (m_barenblatt-1) / (2*m_barenblatt*d)
    if (C - k *(x**2 + y**2) * (1/t**(2*beta)) > 0):
		if (m_barenblatt == 2) :
			hess_sol_barenblatt[0,0] = 1/(t**alpha) * (-2*k*(1/t**(2*beta))) 
			hess_sol_barenblatt[0,1] = 0. 
			hess_sol_barenblatt[1,0] = 0.
			hess_sol_barenblatt[0,0] = 1/(t**alpha) * (-2*k*(1/t**(2*beta))) 
		else :
			# A MODIFIER POUR CAS GENERAL
			hess_sol_barenblatt[0,0] = 1/(t**alpha) * (-2*k*(1/t**(2*beta))) 
			hess_sol_barenblatt[0,1] = 0. 
			hess_sol_barenblatt[1,0] = 0.
			hess_sol_barenblatt[0,0] = 1/(t**alpha) * (-2*k*(1/t**(2*beta))) 
    else : 
        hess_sol_barenblatt[:,:] = 0.
    return hess_sol_barenblatt

     

#______________________________________________________________________#
#
###					MODEL DATA
#______________________________________________________________________#

#  ---- Initial density (non adimensionned)----
initialisation = 'deterministe' #'aleatoire'
rhoini = barenblatt # 'barenblatt' is the name of the function defined above
name_solution = 'barenblatt' # 'barenblatt' is the name of the function defined above

#Size of domain in X and Y directions
#Lx1, Lx2, Ly1, Ly2 the min/max length in X and Y directions defined in file 'initial_densities.py' 
Sx = numpy.array([Lx1, Lx2])
Sy = numpy.array([Ly1, Ly2])


#______________________________________________________________________#
#
###					velocity field
#______________________________________________________________________#

flot = 'barenblatt'
### NO VELOCITY FIELD IN BARENBLATT EQUATION 
    
#______________________________________________________________________#
#
###					NUMERICAL DATA
#______________________________________________________________________#

### ----- Grid size for vizualization (??) ----- 
xgmin = Sx[0]-abs(Sx[1]-Sx[0])/5
xgmax = Sx[1]+abs(Sx[1]-Sx[0])/5
ygmin = Sy[0]-abs(Sy[1]-Sy[0])/5
ygmax = Sy[1]+abs(Sy[1]-Sy[0])/5
Ix = numpy.array([xgmin, xgmax])
Iy = numpy.array([ygmin, ygmax])
Lx = xgmax - xgmin
Ly = ygmax - ygmin 


### -----  Particles list that one wants to trace ----- 
Part_suivies = []

### -----  Numerical method choice ----- 
# method = 'part_Prand' : particulaire avec positions initiales aleatoires
# method = 'part_Mrand' : particulaire avec masse aleatoire
# method = 'part' : particulaire avec discretisation de rhoini
# method = 'part_avec_suivi_vol' : on ne recontruit pas les densites mais on deforme les particules
# method = 'SP' : smooth particule, particules de taille epsilon
# method = 'LTP' : deformation des particules
# methode = 'PIC' : SP avec reconstruction PIC champ d'advection
# methode = 'LTPIC' : LTP avec reconstruction PIC champ d'advection
#method =  'SP' # 'LTP' #
#~ method =  'LTP'
method =  'SP'


# 'oui' = we remapp particles based on a criterion
# 'non' = no remapping
indic_remapping = 'non' #'oui' 

#----- parameters criterion for remapping (??) ----- 
l_remap = 6
radius_remap = 2

### ----- Choice of method to compute the matrices Dk ----- 
# D_method = 'explicit' # one updates particle positions Xk AFTER  Dk matrices (by "explicit" expression of Jacobian of flow)
# D_method = 'implicit' # one updates particle positions Xk BEFORE Dk matrices. (by Finite Differences)
# D_method = 'none' # in case of SP method

if (method == 'SP') :
	D_method = 'none'
else :
	#~ D_method = 'explicit'
	D_method = 'hybrid'


#----- Method choice to compute weights wk ----- 
# met_calcul_poids = 'ponctuelle'
# met_calcul_poids = 'quasi_interp' (methode Martin)
# met_calcul_poids = 'quadrature'
#~ met_calcul_poids =  'ponctuelle'
#~ met_calcul_poids =  'quadrature'
met_calcul_poids = 'quasi_interp' 

# set particles on border of domain (yes or no)
part_sur_bords = 'non'
#part_sur_bords = 'oui'


#----- Shape function ----- 
# spline of order != 3 are not implemented...
# exponential shape functions are imlemented BUT (!!) : 
# we did NOT changed the initialization of the weights accordingly ... -_______-
### TO BE MODIFIED : REALLY SHITTY (!!)
bool_phi_radial = 1
deg_spline = 3

if (bool_phi_radial == 1) : # phi is radial ( so far <=> phi is a spline
    calculsfor_rec_modf90.set_phi_radial(bool_phi_radial)    
    calculsfor_rec_modf90.set_degre_spline(deg_spline)
    if (deg_spline==1) :
        phi = b1 
        radius_phi = 1 # 2 # rayon du support de phi
        radius_phi=1.
    elif (deg_spline == 3) :
        phi = b3 
        derivative_phi = derivative_b3
        phi_2d = b3spline2d
        phi_2d_h = b3spline2d_h
        grad_phi = grad_b3spline2d
        grad_phi_h = grad_b3spline2d_h
        radius_phi = 2 # 2 # rayon du support de phi
    else :
        print "ERROR degre_spline value must be either 1 or 3. Current degre_spline = " , deg_spline

else : # not a radial function
    # phi(x,y) = exp(-1/(1-x**2-y**2)) if (x**2 + y**2) < 1 , 
    # phi(x,y) = 0 otherwise 
    calculsfor_rec_modf90.set_phi_radial(bool_phi_radial)
    calculsfor_rec_modf90.set_degre_spline(deg_spline)
    
    phi = exp_compact #radius_phi is declared when we call that function
    radius_phi = 1.
    


#----- Particles ----- CHECK THE ALL PART (??) (!!)
Nx = 50 # initial number of numerical particles in X direction
Ny = 50 # initial number of numerical particles in Y direction
#~ Nx = 100 # initial number of numerical particles in X direction
#~ Ny = 100 # initial number of numerical particles in Y direction
Nini = Nx * Ny  #  total number of particlesdt
# initial distance between particules :
# ici : on ne met pas de particules au bord
# /max(1,Nx-1) si on  met  des particules au bord
# /(Nx+1) si on ne met pas des particules au bord
if  method == 'part_avec_suivi_vol':
    hx = 0.06
    hy = 0.06
else:
    hx = Lx / max(1, Nx - 1)  # pour visu mettre 0.06 avec N=4
    hy = Ly / max(1, Ny - 1)

if( part_sur_bords == 'non') :
	hx_remap = (Sx[1] - Sx[0]) / (Nx)
	hy_remap = (Sy[1] - Sy[0]) / (Ny)
else : 
	hx_remap = (Sx[1] - Sx[0]) / (Nx-1)
	hy_remap = (Sy[1] - Sy[0]) / (Ny-1)


### Define epsilon (shape function "size") 
#~ h = max((Sx[1] - Sx[0]) / float(Nx-1) , (Sy[1] - Sy[0]) / float(Ny-1))
h = max(hx_remap,hy_remap)

epsilon = 0.0106062528555


# ----- Time discretisation -----
dt = dt_TEMPLATE #1e-0001   #Initial time step
time_scheme = 'euler_explicit'
#time_scheme = 'RK4'
#~ time_scheme = 'middle_point'

### D_method_scheme should disappear in new version (!!)
if (time_scheme == 'euler_explicit') :
	D_method_scheme = 'euler_explicit'

if (time_scheme == 'RK2') :
	D_method_scheme = 'RK2'
if (time_scheme == 'middle_point') :
	D_method_scheme = 'middle_point'
    
# ----- Explosion of particles position moving -----
#~ explosion = "yes" # study explosion or not
explosion="no"
### criterion to know if particles epxloded in X and/or Y direction
### basically : if particles moved more than dx in x one time step then it exploded (resp Y)
dx = hx_remap
dy = hy_remap

#  ---- Reconstruction parameter ------
Nmesh = 500  # parameter of reconstructed function (500 initially)
Nmesh_visu3D = 50
phi_grid = b3 # pour la quasi interpolation (utilise ? a verifier) (??)
radius_phi_grid = 2   #(??)
#  m=7 #nb de points pour la quadrature de W*varphi
#  l=3*N+1 #number of points for the numerical quadrature


#  ======================= output data   ======================================
animgif = 0
animmp4 = 0
position_particles = 1
#~ exportdata = 1
#~ param = -1	   # we plot every time step
param = 200
pictini = 0
pictfin = 0
trace_particules = 'oui'
#trace_particules = 'non'
draw_sol = 'oui'        
#~ draw_sol = 'non'   
save_data = 'non'

# ======================= write data of run   ======================================
myfile=open('./out/'+str(outfile_name),'a')
myfile.write("Tini = "+str(Tini)+'\n')
myfile.write("Tmax = "+str(Tmax)+'\n')
myfile.write("flot = "+str(flot)+'\n')
myfile.write("rhoini = "+str(rhoini)+'\n')
myfile.write("name_solution = "+str(name_solution)+'\n')
myfile.write("time_scheme = "+str(time_scheme)+'\n')
myfile.write("dt = "+str(dt)+'\n')
myfile.write("Nx = "+str(Nx)+'\n')
myfile.write("Ny = "+str(Ny)+'\n')
myfile.write("method = "+str(method)+'\n')
if ( (method == 'SP') | (D_method == 'hybrid') ) : 
	myfile.write("epsilon = "+str(epsilon)+'\n')

myfile.write("D_method = "+str(D_method)+'\n')
myfile.write("indic_remapping = "+str(indic_remapping)+'\n')
myfile.close()

#function to test equalities of floats (maybe move it somewhere else on the code)     
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)            
