#!/usr/bin/env python
# - * - codin: utf-8 -*-
# derniere modif le 09-10-15 par Frederique
from __future__ import division
import numpy  # as np
import scipy  # as sp
from scipy import integrate
import math  # ou : from math import *
import time
from scipy import *
import os
from random import *
from scipy.integrate import dblquad


#from pylab import *

#~ from data_code_LTP_2D import *

import config

#from _calculsfor_f90 import  calculsfor_ini_modf90  # Appel au module FORTRAN 90 compile

       
#===== Appel de l initialisation =====
def initialize_avec_python(sx, sy, nx, ny, rhoini): 
    N = nx * ny
    X = numpy.zeros([2, N])
    M = numpy.array([0. for i in range(N)])    
    mat =  init_poids_et_position(rhoini, sx, sy, nx, ny)
    X[0, :] = mat[0, :]
    X[1, :] = mat[1, :]
    M = numpy.array(mat[2, :])    
    return X, M

def initialize_avec_f2py(sx, sy, nx, ny):
    #initialisation sans mettre de particules sur les bords,
    #calcul des position et poids deterministes
    #utilisation de rhoini2 avec r=1  ----> a generaliser
    calculsfor_ini_modf90.set_num_rhoini(num_rhoini, r, Lx1, Lx2, Ly1, Ly2)
    N_part = nx*ny
    Domaine = numpy.zeros((2,2))
    Domaine[0,:] = sx
    Domaine[1,:] = sy
    calculsfor_ini_modf90.alloc_vect(N_part)
    calculsfor_ini_modf90.initialize_spline(Domaine, nx, ny, deg_spline)
    Xp = calculsfor_ini_modf90.xpart
    X = Xp.copy()
    Mp = calculsfor_ini_modf90.mpart
    M = Mp.copy()
    calculsfor_ini_modf90.dealloc_vect()
    return X, M

#==== Calcul approche integrale avec methode trapeze =====
def int2d(f, x1, x2, y1, y2, l):
#  l=nb de points pour la quadrature
    dx = (x2 - x1) / float(l - 1)
    dy = (y2 - y1) / float(l - 1)
    xsi1 = [x1 + (2 * i + 1.) / 2. * dx for i in range(l - 1)]
    xsi2 = [y1 + (2 * i + 1.) / 2 * dy for i in range(l - 1)]
    fj = numpy.array([0. for i in range(l - 1)])
    for i in range(l - 1):
        fi = numpy.array([f(xsi1[i], xsi2[j]) for j in range(l - 1)])
        fj[i] = sum(fi)
    return sum(fj) * dx * dy
    
def calculmasse(f, sx,sy,l):
    return int2d(f, sx[0], sx[1], sy[0], sy[1], l)


#===== Calcul des positions initiales des particules =====
def init_poids_et_position(f, sx, sy, nx, ny):
    #----- cas ou on discretise aleatoirement sans rhoini -----
    if config.method == 'part_Mrand':
        return initialize_part_Mrand(sx, sy, nx, ny)
    elif config.method == 'part_Prand':
        return initialize_part_Prand(sx, sy, nx, ny) 
    #----- autres cas -----
    elif config.part_sur_bords == 'non':        
        return initializepartssbords(f, sx, sy, nx, ny)
    else:  
        return initializepartavecbords(f, sx, sy, nx, ny)

#----- Pas de particules sur les bords -----
def initializepartssbords(f, sx, sy, nx, ny):  #pas de particules sur les bords du domaine
    hx = (sx[1] - sx[0]) / float(nx)  #  distance between particles
    hy = (sy[1] - sy[0]) / float(ny)
    z = zeros([3, nx * ny])	
    k = 0
    for i in range(ny):		
        if ny == 1:
            y = (sy[1] + sy[0]) / 2.
        else:
            y = sy[1] - (i + 0.5) * hy
        for j in range(nx):
            if nx == 1:
                x = (sx[1] + sx[0]) / 2.
            else:
                x = sx[0] + (j + 0.5) * hx
            z[0, k] = x
            z[1, k] = y
            if (config.met_calcul_poids == 'quasi_interp') and (config.deg_spline ==3):
                #t0 = f(x, y) * (8/6)**2
                #t01 = -8/36 *(f(x-hx, y)+f(x+hx, y)+f(x, y+hy)+f(x, y-hy))
                #t11 =1/36* (f(x-hx, y-hy)+ f(x-hx, y+hy)+f(x+hx, y-hy)+f(x+hx, y+hy) )
                ##z[2, k] = hx * hy *(t0+t01+t11)  
                #z[2, k] = hx * hy * max(t0+t01+t11,0) 
                ### RAJOUT DES '.' , e.g. : 8. / 6.
                t0 = f(x, y) * (8./6.)**2
                t01 = -8./36. *(f(x-hx, y)+f(x+hx, y)+f(x, y+hy)+f(x, y-hy))
                t11 =1./36.* (f(x-hx, y-hy)+ f(x-hx, y+hy)+f(x+hx, y-hy)+f(x+hx, y+hy) )
                #z[2, k] = hx * hy *(t0+t01+t11)  
                z[2, k] = hx * hy * max(t0+t01+t11,0) 
            elif  config.met_calcul_poids =='quadrature':
                z[2, k] = int2d(f,x-hx/2.,x+hx/2.,y-hy/2.,y+hy/2.,10)
                #~ print "f (x,y) = " , f(x,y)
                #~ print "int2d with dblquad = " , dblquad(f, y-hy/2., y+hy/2., lambda x: x-hx/2., lambda x: x+hx/2)[0]
                #~ print "int2d with trapezes = " , int2d(f,x-hx/2.,x+hx/2.,y-hy/2.,y+hy/2.,10)
            else: # inclus (met == 'formulespline') and (deg_spline ==1) 
				#config.met_calcul_poids = 'ponctuelle'                
                z[2, k] = hx * hy * f(x,y)
            k += 1
            
    return z

#----- Des particules sur les bords ----- > A DEBUGGER
# met_calcul_poids = 'ponctuelle' ou 'integration' mais pas 'quasi_interp'
def initializepartavecbords(f, sx, sy, nx, ny):  #particules sur les bords du domaine
    methode='ponctuelle'  
    hx = (sx[1] - sx[0]) / float(min(1, nx-1))  #  distance between particles
    hy = (sy[1] - sy[0]) / float(min(1, ny-1))
    z = zeros([3, nx * ny])
    k = 0
    for i in range(ny):
        if ny == 1:
            y = (sy[1] + sy[0]) / 2.
        else:
            y = sy[1] - i * hy
        for j in range(nx):
            if nx == 1:
                x = (sx[1] + sx[0]) / 2.
            else:
                x = sx[0] + j * hx
            z[0, k] = x
            z[1, k] = y
            # initialisation des poids avec valeur ponctuelle
            if (i > 0) and (i < config.Nx) and (j > 0) and (j < config.Nx):
                z[2, k] = f(x, y) * hx * hy
            elif ((i == 1) or (i == config.Nx)) and ((j == 1) or (j == config.Nx)):
                z[2, k] = f(x, y) * hx * hy * 7. / 4.
            else:
                if  config.met_calcul_poids == 'ponctuelle':
                    #print  'initialisation des poids ponctuelle'
                    z[2, k] == f(x, y) * hx * hy * 3. / 4.                        
                else: # initialisation des poids avec  quadrature
                    z[2, k] = int2d(f,x-hx/2.,x+hx/2,y-hy/2.,y+hy/2.,100)
            k += 1
    return z

# ----- Pour des tests => ne sert plus -----
def initialize_part_Mrand(sx, sy, nx, ny):
    #ici on met des particules sur les bords aussi
    hx = (sx[1] - sx[0]) / float(min(1, nx-1))  #  distance between particles
    hy = (sy[1] - sy[0]) / float(min(1, ny-1))
    z = zeros([3, nx * ny])
    k = 0
    for i in range(ny):
        if ny == 1:
            y = (sy[1] + sy[0]) / 2.
        else:
            y = sy[1] - i * hy
        for j in range(nx):
            if nx == 1:
                x = (sx[1] + sx[0]) / 2.
            else:
                x = sx[0] + j * hx
            z[0, k] = x
            z[1, k] = y
            z[2, k] = random()
            k += 1
    somme_poids = sum(z[2,:])
    z[2,:] = z[2,:]/somme_poids
    return z
    
# ----- Initialisation aleatoire des particules, pour des tests visuels -----
def initialize_part_prand(sx, sy, nx, ny):
    z = zeros([3, nx * ny])
    k = 0
    for i in range(ny):
        for j in range(nx):
            z[0, k] = sx[0]+random()*(sx[1] - sx[0]) 
            z[1, k] = sy[0]+random()*(sy[1] - sy[0]) 
            z[2, k] = 1./(nx*ny)
            k +=1
    return z

#===== Suppression des particules dont le poids est quasi-nul =====    
def supprime_particule(X,M):
    N = len(M)
    poid_tot = sum(M)
    tol = poid_tot/N*0.01 
    Nnew = N
    M1 = numpy.array([0. for i in range(N)])
    X1 =zeros((2,N))
    knew = 0
    for k in range(N):
        if abs(M[k])> tol:
            M1[knew] = M[k]
            X1[:,knew] = X[:,k]
            knew += 1            
        else:
            Nnew += -1

    Xnew = X1[:,0:Nnew]
    Mnew = M1[0:Nnew]
    return Nnew, Xnew, Mnew

    
