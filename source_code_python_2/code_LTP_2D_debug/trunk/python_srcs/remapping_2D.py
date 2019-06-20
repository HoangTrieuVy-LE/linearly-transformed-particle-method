# -*- coding: utf-8 -*-
# derniere modif le 09-10-15 par Frederique
from __future__ import division
import numpy  #as np

import config

#~ from shapefunction2D import *
#~ from initialisation_2D import *
#~ from _calculsfor_f90 import  calculsfor_ini_modf90

#~ from data_code_LTP_2D import *

#~ from reconstruction_densites_2D import *
#~ from graphes2D import *


def remapping_avec_python(Nxnew, Nynew, Xold, Mold, Dold,t):
    l = l_remap
    Nold = len(Mold)
    Domaine = numpy.zeros((2,2))
    Domaine[0,0] = min(Xold[0,:])
    Domaine[0,1] = max(Xold[0,:])
    Domaine[1,0] = min(Xold[1,:])
    Domaine[1,1] = max(Xold[1,:])
    N_part = Nxnew * Nynew    
    xpart = numpy.zeros((2,N_part))
    mpart = numpy.array([0. for i in range(N_part)])
    xmin = Domaine[0,0]
    xmax = Domaine[0,1]
    ymin = Domaine[1,0] 
    ymax = Domaine[1,1]
    hx = (xmax - xmin)/float(Nxnew-1)  # distance between particles
    hy = (ymax - ymin)/float(Nynew-1)
    deltax = hx/(l+1)
    deltay = hy/(l+1)
    NgrilleX = Nxnew +l*(Nxnew-1)+2
    NgrilleY = Nynew +l*(Nynew-1)+2 
    Xg = numpy.array([xmin-deltax + (k)*deltax for k in range(NgrilleX)])
    Yg = numpy.array([ymin-deltay + k*deltay for k in range(NgrilleY)])
    ro = rec_densite_grille_unif_LTP(Xold, Mold, Dold, Xg, Yg)
    k = 0
    for j in range(Nynew):
        y =  ymin + hy*(j)
        indice_Yg_y = (j)*(l+1)+2-1  
        for i in range(Nxnew):
            x = xmin + (i )*hx
            indice_Xg_x = (i)*(l+1)+2-1
            xpart[0,k] = x
            xpart[1,k] = y
            if (deg_spline ==1):
                 z = ro[indice_Xg_x, indice_Yg_y]
                 mpart[k] =  hx * hy * z 
            else: 
                z00 = ro[indice_Xg_x, indice_Yg_y]
                zm0 = ro[indice_Xg_x-1, indice_Yg_y]
                zp0 = ro[indice_Xg_x+1, indice_Yg_y]
                z0m = ro[indice_Xg_x, indice_Yg_y-1]
                z0p = ro[indice_Xg_x, indice_Yg_y+1]
                zpp = ro[indice_Xg_x+1, indice_Yg_y+1]
                zmm = ro[indice_Xg_x-1, indice_Yg_y-1]
                zpm = ro[indice_Xg_x+1, indice_Yg_y-1]
                zmp = ro[indice_Xg_x-1, indice_Yg_y+1]
                t00 = z00 * (8./6.)**2
                t01 = -8./36. *(zm0+zp0+z0m+z0p)
                t11 =1./36.* (zmm+zpp+zmp+zpm)
                mpart[k] = hx * hy * max(t00+t01+t11,0.)
            k = k+1
    print 'coucou'
    fait_series_dessins(Xg, Yg, ro, int(t), t, 'graphe_remapping')
    print 'somme ancien poids', sum(Mold)
    print 'somme nouveau poids', sum(mpart)
    #N_part_new, X, M = supprime_particule(xpart,mpart)
    N_part_new = N_part
    X = xpart
    M = mpart
    D = numpy.zeros([4, N_part_new]) #matrice de deformation
    D[0,:] = 1.  # coeff (1,1) de Dk^0
    D[3,:] = 1.       
    return X, M, D , hx , hy    

#We change the window of remapping, so that we keep a grid of same initial size hx_remap , hy_remap
def remapping_avec_python_test(hx,hy, Xold, Mold, Dold,t):
    l = l_remap
    Nold = len(Mold)
    Domaine = numpy.zeros((2,2))
    Domaine[0,0] = min(Xold[0,:])
    Domaine[0,1] = max(Xold[0,:])
    Domaine[1,0] = min(Xold[1,:])
    Domaine[1,1] = max(Xold[1,:])
    xmin = Domaine[0,0]
    xmax = Domaine[0,1]
    ymin = Domaine[1,0] 
    ymax = Domaine[1,1]
    
    Nxnew=int((xmax-xmin)/hx + 1)
    eps_x=(hx*Nx-(xmax-xmin))/2
    xmin=xmin-eps_x
    xmax=xmax+eps_x
    
    Nynew=int((ymax-ymin)/hy + 1)
    eps_y=(hy*Ny-(ymax-ymin))/2
    ymin=ymin-eps_y
    ymax=ymax+eps_y
    
    N_part=Nxnew*Nynew
    xpart = numpy.zeros((2,N_part))
    mpart = numpy.array([0. for i in range(N_part)])
    
    #hx = (xmax - xmin)/float(Nxnew-1)  # distance between particles
    #hy = (ymax - ymin)/float(Nynew-1)
    deltax = hx/(l+1)
    deltay = hy/(l+1)
    NgrilleX = Nxnew +l*(Nxnew-1)+2
    NgrilleY = Nynew +l*(Nynew-1)+2 
    Xg = numpy.array([xmin-deltax + (k)*deltax for k in range(NgrilleX)])
    Yg = numpy.array([ymin-deltay + k*deltay for k in range(NgrilleY)])
    ro = rec_densite_grille_unif_LTP(Xold, Mold, Dold, Xg, Yg)
    k = 0
    for j in range(Nynew):
        y =  ymin + hy*(j)
        indice_Yg_y = (j)*(l+1)+2-1  
        for i in range(Nxnew):
            x = xmin + (i )*hx
            indice_Xg_x = (i)*(l+1)+2-1
            xpart[0,k] = x
            xpart[1,k] = y
            if (deg_spline ==1):
                 z = ro[indice_Xg_x, indice_Yg_y]
                 mpart[k] =  hx * hy * z 
            else: 
                z00 = ro[indice_Xg_x, indice_Yg_y]
                zm0 = ro[indice_Xg_x-1, indice_Yg_y]
                zp0 = ro[indice_Xg_x+1, indice_Yg_y]
                z0m = ro[indice_Xg_x, indice_Yg_y-1]
                z0p = ro[indice_Xg_x, indice_Yg_y+1]
                zpp = ro[indice_Xg_x+1, indice_Yg_y+1]
                zmm = ro[indice_Xg_x-1, indice_Yg_y-1]
                zpm = ro[indice_Xg_x+1, indice_Yg_y-1]
                zmp = ro[indice_Xg_x-1, indice_Yg_y+1]
                t00 = z00 * (8./6.)**2
                t01 = -8./36. *(zm0+zp0+z0m+z0p)
                t11 =1./36.* (zmm+zpp+zmp+zpm)
                mpart[k] = hx * hy * max(t00+t01+t11,0.)
            k = k+1
    print 'coucou'
    fait_series_dessins(Xg, Yg, ro, int(t), t, 'graphe_remapping')
    print 'somme ancien poids', sum(Mold)
    print 'somme nouveau poids', sum(mpart)
    #N_part_new, X, M = supprime_particule(xpart,mpart)
    N_part_new = N_part
    X = xpart
    M = mpart
    D = numpy.zeros([4, N_part_new]) #matrice de deformation
    D[0,:] = 1.  # coeff (1,1) de Dk^0
    D[3,:] = 1.       
    return X, M, D      

def remapping_avec_f2py(Nxnew, Nynew, Xold, Mold, Dold): # A DEBUGGER
    Nold = len(Mold)
    Domaine = numpy.zeros((2,2))
    Domaine[0,0] = min(Xold[0,:])
    Domaine[0,1] = max(Xold[0,:])
    Domaine[1,0] = min(Xold[1,:])
    Domaine[1,1] = max(Xold[1,:])
    N_part = Nxnew * Nynew
    calculsfor_ini_modf90.alloc_vect(N_part)
    calculsfor_ini_modf90.set_param_remap(Domaine,  deg_spline,  l_remap)
    print 'Mold', Mold
    print 'len(Mold)', len(Mold)
    calculsfor_ini_modf90.remap(Xold, Mold, Dold, Nxnew, Nynew, Nold)
    Xp = calculsfor_ini_modf90.xpart
    X = Xp.copy()
    Mp = calculsfor_ini_modf90.mpart
    M = Mp.copy()
    calculsfor_ini_modf90.dealloc_vect()
    N_part_new, X, M = supprime_particule(X,M)
    D = numpy.zeros([4, N_part_new]) #matrice de deformation
    D[0,:] = 1.  # coeff (1,1) de Dk^0
    D[3,:] = 1. 
    return X, M, D 

#===== Methode ou on ne change pas la position des particules -> voir si ca marche =====
def remapping_test_vieux(N, X, M, D, hx, hy):
# on garde la meme position des particules, on ne change que le poids et la  forme
# methode brutale, a ameliorer (super couteux)
    Mnew = numpy.array([0. for i in range(N)])    
    for i in range(N):
        xi = X[0,i]
        yi = X[1,i]    
        for k in range(N):
            Dk = D[:, k].reshape((2, 2))
            Yk = numpy.dot(Dk,numpy.array([X[0,k]-xi,X[1,k]-yi]) )                
            phix = phi(Yk[0]/hx)
            phiy = phi(Yk[1]/hy) 
            Mnew[i] += M[k] * phix * phiy
            
    Dnew = numpy.zeros([4, N]) #matrice de deformation
    Dnew[0,:] = 1.  # coeff (1,1) de Dk^0
    Dnew[3,:] = 1.  # coeff (2,2) de Dk^0 
    return Mnew, Dnew        

