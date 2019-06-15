# -*- coding: utf-8 -*-
# derniere modif le 09-10-15 par Frederique

from __future__ import division
import numpy  #as np
#~ from scipy import *  #as sp
import scipy
import config

#~ from shapefunction2D import *
#~ from data_code_LTP_2D  import *
#~ from calculs_2D import *

from _calculsfor_f90 import  calculsfor_rec_modf90 # Appel au module FORTRAN 90 compile

#----- Cas de la methode SP -----
def rec_densite_grille_unif_sp(X, M, epsilon, Xg, Yg):
    #return rec_densite_grille_unif_sp_avec_f2py(X, M, epsilon, Xg, Yg)
    return rec_densite_grille_unif_sp_avec_python(X, M, epsilon, Xg, Yg)

def rec_densite_grille_unif_sp_avec_f2py(X, M, epsilon, Xg, Yg):
    # derniere modif 07-10-15
    # interface avec fortran
    NgrilleX = len(Xg)
    NgrilleY = len(Yg)
    N = len(M)
    calculsfor_rec_modf90.alloc_ro(NgrilleX, NgrilleY)
    calculsfor_rec_modf90.set_degre_spline(deg_spline)
    calculsfor_rec_modf90.set_hx_hy(hx, hy)
    calculsfor_rec_modf90.set_grille(Xg, Yg, NgrilleX, NgrilleY)
    #--- routine f2py specifique : ---> ne marche pas, a debugguer    
    #calculsfor_rec_modf90.rec_sp(X, M, epsilon, N)  
    #-----
    #----- utiisation de celle faite pour reconstruction LTP (un peu plus lente)
    D = numpy.zeros((4,N))
    D[0, :] = 1.
    D[3, :] = 1.
    calculsfor_rec_modf90.rec_ltp(X, M, D, N)  
    #-----
    Ro = calculsfor_rec_modf90.ro
    Rho = Ro.copy()
    calculsfor_rec_modf90.dealloc_ro()
    return  Rho
    
def rec_densite_grille_unif_sp_avec_python(X, M, epsilon, Xg, Yg):
    # version 04/08/15, avec les x qui correspond au premier indice de la matrice
    # methode avec depot
    deltax = abs(Xg[1] - Xg[0])
    deltay = abs(Yg[1] - Yg[0])    
    xmin = min(Xg)
    ymin = min(Yg)
    mx = scipy.size(Xg)
    my = scipy.size(Yg)
    R = numpy.zeros((mx,my))
    N = scipy.size(M) 
    qx = int(config.radius_phi * epsilon / deltax) + 1
    qy = int(config.radius_phi * epsilon / deltay) + 1
    for k in range(N):
        ik = int((X[0,k] - xmin) / deltax)
        jk = int((X[1,k] - ymin ) / deltay)
        depotx = numpy.zeros((mx,1))
        depoty = numpy.zeros((1,my))        
        for i in range(max(0,ik- qx), min(ik+qx+1, mx)): 
            #depotx[0,j] = phi((X[0,k]-Xg[j]) / epsilon)
            depotx[i,0] = config.phi((X[0,k]-Xg[i]) / epsilon)
        for j in range(max(jk-qy,0), min(jk+qy+1,my)):
            #depoty[i,0] = phi((X[1,k] -Yg[i]) / epsilon) 
            depoty[0,j] = config.phi((X[1,k] -Yg[j]) / epsilon) 
        Matricedepot =   depotx * depoty
        R += M[k] / (epsilon * epsilon) * Matricedepot
    return  R
 
#----- Cas de la methode LTP -----  
def rec_densite_grille_unif_LTP(X, M, D, Xg, Yg):
    return rec_densite_grille_unif_LTP_avec_f2py(X, M, D, Xg, Yg)
    #return rec_densite_grille_unif_LTP_avec_python(X, M, D, Xg, Yg)
    #return rec_densite_grille_unif_LTP_avec_f2py_phi_support_non_compact(X, M, D, Xg, Yg)
    #~ return rec_densite_grille_unif_LTP_avec_f2py_v2(X, M, D, Xg, Yg)
 
def rec_densite_grille_unif_LTP_avec_f2py(X, M, D, Xg, Yg):
    # derniere modif 07-10-15
    # interface avec fortran
    NgrilleX = len(Xg)
    NgrilleY = len(Yg)
    N = len(M)
    calculsfor_rec_modf90.alloc_ro(NgrilleX, NgrilleY)
    calculsfor_rec_modf90.set_phi_radial(config.bool_phi_radial)
    calculsfor_rec_modf90.set_degre_spline(config.deg_spline)
    calculsfor_rec_modf90.set_hx_hy(config.hx_remap, config.hy_remap)
    calculsfor_rec_modf90.set_grille(Xg, Yg, NgrilleX, NgrilleY)
    calculsfor_rec_modf90.rec_ltp(X, M, D, N)
    #calculsfor_rec_modf90.rec_sansdepot(X, M, D, N)
    Ro = calculsfor_rec_modf90.ro
    Rho = Ro.copy()
    calculsfor_rec_modf90.dealloc_ro()
    return  Rho

def rec_densite_grille_unif_LTP_avec_f2py_v2(X, M, D, Xg, Yg):
    # derniere modif 07-10-15
    # interface avec fortran
    NgrilleX = len(Xg)
    NgrilleY = len(Yg)
    N = len(M)
    calculsfor_rec_modf90.alloc_ro(NgrilleX, NgrilleY)
    calculsfor_rec_modf90.set_phi_radial(bool_phi_radial)
    calculsfor_rec_modf90.set_degre_spline(deg_spline)
    calculsfor_rec_modf90.set_hx_hy(hx, hy)
    calculsfor_rec_modf90.set_grille(Xg, Yg, NgrilleX, NgrilleY)
    calculsfor_rec_modf90.rec_ltp_v2(X, M, D, N)
    #calculsfor_rec_modf90.rec_sansdepot(X, M, D, N)
    Ro = calculsfor_rec_modf90.ro
    Rho = Ro.copy()
    calculsfor_rec_modf90.dealloc_ro()
    return  Rho



def rec_densite_grille_unif_LTP_avec_f2py_phi_support_non_compact(X, M, D, Xg, Yg):
    # derniere modif 07-10-15
    # interface avec fortran
    NgrilleX = len(Xg)
    NgrilleY = len(Yg)
    N = len(M)
    calculsfor_rec_modf90.alloc_ro(NgrilleX, NgrilleY)
    calculsfor_rec_modf90.set_degre_spline(deg_spline)
    calculsfor_rec_modf90.set_hx_hy(hx, hy)
    calculsfor_rec_modf90.set_grille(Xg, Yg, NgrilleX, NgrilleY)
    calculsfor_rec_modf90.rec_ltp_phi_support_non_compact(X, M, D, N)
    #calculsfor_rec_modf90.rec_sansdepot(X, M, D, N)
    Ro = calculsfor_rec_modf90.ro
    Rho = Ro.copy()
    calculsfor_rec_modf90.dealloc_ro()
    return  Rho


def rec_densite_grille_unif_LTP_avec_python(X, M, D, Xg, Yg):
    # derniere modif 06-08-15
    # methode avec depot
    xmin = min(Xg)
    ymin = min(Yg)
    mx = size(Xg)
    my = size(Yg)
    deltax = abs(Xg[1] - Xg[0])
    deltay = abs(Yg[1] - Yg[0]) 
    R = numpy.zeros((mx,my))
    N = size(M) 
    # q = int(radius_phi * epsilon / delta) + 1
    for k in range(N):
        a = D[0 , k]
        b = D[1 , k]   
        c = D[2 , k]
        d = D[3 , k]
        #print 'Dk', D[:,k]   
        det_Dk = abs(a*d-b*c)
        if det_Dk ==0:
            print ('probleme de Dk non inversible')
        rayonkx =  hx /det_Dk*max(abs(d)+abs(b),abs(c)+abs(a))
        rayonky = hy/det_Dk*max(abs(d)+abs(b),abs(c)+abs(a))
        qkx = int(radius_phi * rayonkx / deltax) + 1
        qky = int(radius_phi * rayonky / deltay) + 1
        ik = int((X[0,k] - xmin) / deltax)
        jk = int((X[1,k] - ymin ) / deltay)
        for j in range(max(0,jk- qky), min(jk+qky+1, my)): 
            yj = Yg[j]          
            for i in range(max(ik-qkx,0), min(ik+qkx+1,mx)): 
                xi = Xg[i]
                phix = phi( (a*(X[0,k]-xi)+b*(X[1,k] - yj)) / hx)
                phjy = phi( (c*(X[0,k]-xi)+d*(X[1,k] - yj)) / hy)
                R[i,j] += M[k]*det_Dk/(hx*hy) * phix * phjy
    print ('maxR', R.max())
    return  R

#----- Pour visu de la forme des particules   
def rec_forme_particule(k, Xk, Dk, Mk, Nmesh):
    Sxk = numpy.array([Xk[0]-5*hx, Xk[0]+5*hx])
    Syk = numpy.array([Xk[1]-5*hy, Xk[1]+5*hy])                
    [Xgk, Ygk] = make_grid_unif(Sxk, Syk, Nmesh)
    Dk = Dk.reshape((4, 1))   
    Xk = Xk.reshape((2,1))
    Mk = numpy.array([Mk])
    Z = rec_densite_grille_unif_LTP(Xk, Mk, Dk, Xgk, Ygk)
    return Xgk, Ygk, Z
 
   
#==================================================================
# fonctions non utilisees



def rec_1D_sp(Y, M, epsilon, Yg):
    # version 04/08/15, avec les x qui correspond au premier indice de la matrice
    # methode avec depot
    deltay = abs(Yg[1] - Yg[0])    
    ymin = min(Yg)
    my = scipy.size(Yg)
    R = numpy.array([0. for i in range(my)])
    N = scipy.size(M) 
    qy = int(radius_phi * epsilon / deltay) + 1
    for k in range(N):
        jk = int((Y[k] - ymin ) / deltay)
        for j in range(max(jk-qy,0), min(jk+qy+1,my)):
            #depoty[i,0] = phi((X[1,k] -Yg[i]) / epsilon) 
            R[j] += M[k] / (epsilon * epsilon) *phi((Y[k] -Yg[j]) / epsilon) 
    return  R

    
def rec_densite_grille_unif_LTP2(X, M, D, Xg, Yg):
    # derniere modif 10-09-15
    # ne marche pas (mauvaise reconstruction)
    eps_precis = 1e-2
    xmin = min(Xg)
    ymin = min(Yg)
    mx = size(Xg)
    my = size(Yg)
    print ('mx', mx)
    print ('my', my)
    deltax = abs(Xg[1] - Xg[0])
    deltay = abs(Yg[1] - Yg[0])    #print 'delta', delta
    R = numpy.zeros((mx,my))
    N = size(M) 
    # q = int(radius_phi * epsilon / delta) + 1
    for k in range(N):
        a = D[0 , k]
        b = D[1 , k]   
        c = D[2 , k]
        d = D[3 , k]
        det_Dk = abs(a*d-b*c)
        xk = X[0,k]
        yk = X[1,k]
        ik = int((xk - xmin) / deltax)
        jk = int((yk - ymin ) / deltay)        
        phix = phi( (a*(xk-Xg[ik])+b*(yk - Yg[jk])) / hx)
        # --- premier quadran
        i = 0
        while (abs(phix)>eps_precis) and (ik+i<mx):
            j = 0
            phjy = phi( (c*(xk-Xg[ik+i])+d*(yk - Yg[jk])) / hy) 
            while  (abs(phjy)>eps_precis) and (jk+j<my):
                phjy = phi( (c*(xk-Xg[ik+i])+d*(yk - Yg[jk+j])) / hy) 
                phix = phi( (a*(xk-Xg[ik+i])+b*(yk - Yg[jk+j])) / hx)    
                R[i+ik,j+jk] += M[k]*det_Dk/(hx*hy) * phix * phjy               
                j += 1
                #print 'j+jk', j+jk            
            i += 1    
            #print 'i+ik', i+ik
        # --- 
        phix = phi( (a*(xk-Xg[ik])+b*(yk - Yg[jk])) / hx)
        i = 0
        while (abs(phix)>eps_precis) and (ik-i>=0):
            j = 0
            phjy = phi( (c*(xk-Xg[ik+i])+d*(yk - Yg[jk])) / h) 
            while  (abs(phjy)>eps_precis) and (jk+j<my):
                phjy = phi( (c*(xk-Xg[ik+i])+d*(yk - Yg[jk+j])) / hy) 
                phix = phi( (a*(xk-Xg[ik+i])+b*(yk - Yg[jk+j])) / hx)                    
                R[i,j] += M[k]*det_Dk/(hx*hy) * phix * phjy
                j += 1
                #print 'j+jk', j+jk               
            i -= 1 
            #print 'i+ik', i+ik
        # --- 
        phix = phi( (a*(xk-Xg[ik])+b*(yk - Yg[jk])) / hx)
        i = 0
        while (abs(phix)>eps_precis) and (ik+i<mx):
            j = 0
            phjy = phi( (c*(xk-Xg[ik+i])+d*(yk - Yg[jk])) / hy) 
            while  (abs(phjy)>eps_precis) and (jk-j>=0):
                phjy = phi( (c*(xk-Xg[ik+i])+d*(yk - Yg[jk+j])) / hy) 
                phix = phi( (a*(xk-Xg[ik+i])+b*(yk - Yg[jk+j])) / hx) 
                R[i,j] += M[k]*det_Dk/(hx*hy) * phix * phjy
                j -= 1
                #print 'j+jk', j+jk               
            i += 1     
            #print 'i+ik', i+ik
        # --- 
        phix = phi( (a*(xk-Xg[ik])+b*(yk - Yg[jk])) / hx)
        i = 0
        while (abs(phix)>eps_precis) and (ik-i>=0):
            j = 0
            phjy = phi( (c*(xk-Xg[ik+i])+d*(yk - Yg[jk])) / hy) 
            while  (abs(phjy)>eps_precis) and (jk-j>=0):
                phjy = phi( (c*(xk-Xg[ik+i])+d*(yk - Yg[jk+j])) / hy) 
                phix = phi( (a*(xk-Xg[ik+i])+b*(yk - Yg[jk+j])) / hx)   
                R[i,j] += M[k]*det_Dk/(hx*hy) * phix * phjy
                j -= 1
                #print 'j+jk', j+jk             
            i -= 1     
            #print 'i+ik', i+ik                                   
    return  R
    
def rec_densite_grille_unif_sp_sansdepot(X, M, epsilon, Xg, Yg):
    # methode basique avec double boucle (pour comparaison)
    # modifie le 04-08-15
    mx = len(Xg)
    my = len(Yg)
    R = numpy.zeros((mx, my))
    N = int(len(M))
    for i in range(mx):    
        xi = Xg[i]
        for j in range(my):  
            yj =  Yg[j]
            for k in range(N):
                phix = phi((X[0,k]-xi) / epsilon)
                phiy = phi((X[1,k]-yj) / epsilon)
                R[i,j] += M[k]/(epsilon*epsilon) * phix * phiy
    return R

def rec_densite_grille_unif_LTP_sansdepot(X, M, D, Xg, Yg):
    # methode basique avec double boucle (pour comparaison)
    # modifie le 10-09-15
    mx = len(Xg)
    my = len(Yg)
    R = numpy.zeros((mx, my))
    N = int(len(M))
    for i in range(mx):    
        xi = Xg[i]
        for j in range(my):  
            yj =  Yg[j]
            for k in range(N):
                Dk = D[:, k].reshape((2, 2))
                Yk = numpy.dot(Dk,numpy.array([X[0,k]-xi,X[1,k]-yj]) )                
                phix = phi(Yk[0]/hx)
                phiy = phi(Yk[1]/hy) 

                R[i,j] += M[k]*det_Dk/(hx*hy) * phix * phiy
    return R

def depose_sur_part_sp(V, X, Xg, Yg): #la meme chose que l'on soit en PIC ou LTPIC
    # V 2d-array de taille mx*my  
    # modif le 5/08/15
    # a revoir
    mx = size(Xg)
    my = size(Yg)
    xmin = min(Xg)
    ymin = min(Yg)    
    deltax = abs(Xg[1] - Xg[0])
    deltay = abs(Yg[1] - Yg[0])    
    Vp = 0*X         # numpy.array([0. for x in X])
    N = size(X[0,:]) # pas tres propre, a modifier
    # q = floor(radius_phi_grid * epsilon / delta) + 1
    for k in range(N):
        ik = int((X[0, k] - xmin) / deltax)
        binf_i = max(0, ik - radius_phi_grid - 1)
        bsup_i = min(mx, ik + radius_phi_grid + 2)
        jk = int((X[1, k] - ymin ) / deltay)
        binf_j = max(0, jk - radius_phi_grid - 1)
        bsup_i = min(my, jk + radius_phi_grid + 2)
        for j in range(binf_j,bsup_j):
            yj = Y[j]
            contrib_j=phi_grid((X[1,k]-yj) / deltay)
            for i in range(binf_i,bsup_i): 
                xi= X[i]
                contrib_i=phi_grid((X[0,k]- xi) / deltax)
                Vp[k] += V[i,j] * contrib_i * contrib_j
    return Vp
    
    
def densite_part(X, M):
    Xg = list(X[0, :])
    Xg.sort()
    Yg = list(X[1, :])
    Yg.sort()
    N = len(M)
    R = numpy.zeros((N, N))
    for k in range(N):
        i = Xg.index(X[0, k])
        j = Yg.index(X[1, k])
        R[i,j]+=M[k]
    
    Xgrille = numpy.array(Xg)
    Ygrille = numpy.array(Yg)    
    return Xgrille, Ygrille, R
    
def rec_densite_grille_unif_sp0(X, M, epsilon, Xg, Yg):
    # version 04/08/15, avec les x qui correspond au premier indice de la matrice
    # methode avec depot
    deltax = abs(Xg[1] - Xg[0])
    deltay = abs(Yg[1] - Yg[0])    
    xmin = min(Xg)
    ymin = min(Yg)
    mx = scipy.size(Xg)
    my = scipy.size(Yg)
    R = numpy.zeros((mx,my))
    N = scipy.size(M) 
    qx = int(radius_phi * epsilon / deltax) + 1
    qy = int(radius_phi * epsilon / deltay) + 1
    for k in range(N):
        ik = int((X[0,k] - xmin) / deltax)
        jk = int((X[1,k] - ymin ) / deltay)
        depotx = numpy.zeros((mx,1))
        depoty = numpy.zeros((1,my))        
        for i in range(max(0,ik- qx), min(ik+qx+1, mx)): 
            #depotx[0,j] = phi((X[0,k]-Xg[j]) / epsilon)
            depotx[i,0] = phi((X[0,k]-Xg[i]) / epsilon)
        for j in range(max(jk-qy,0), min(jk+qy+1,my)):
            #depoty[i,0] = phi((X[1,k] -Yg[i]) / epsilon) 
            depoty[0,j] = phi((X[1,k] -Yg[j]) / epsilon) 
        Matricedepot =   depotx * depoty
        R += M[k] / (epsilon * epsilon) * Matricedepot
    return  R

def rec_densite_grille_unif_sp0(X, M, epsilon, Xg, Yg):
    # version 04/08/15, avec les x qui correspond au premier indice de la matrice
    # methode avec depot
    deltax = abs(Xg[1] - Xg[0])
    deltay = abs(Yg[1] - Yg[0])    
    xmin = min(Xg)
    ymin = min(Yg)
    mx = scipy.size(Xg)
    my = scipy.size(Yg)
    R = numpy.zeros((mx,my))
    N = scipy.size(M) 
    qx = int(radius_phi * epsilon / deltax) + 1
    qy = int(radius_phi * epsilon / deltay) + 1
    for k in range(N):
        ik = int((X[0,k] - xmin) / deltax)
        jk = int((X[1,k] - ymin ) / deltay)
        depotx = numpy.zeros((mx,1))
        depoty = numpy.zeros((1,my))        
        for i in range(max(0,ik- qx), min(ik+qx+1, mx)): 
            #depotx[0,j] = phi((X[0,k]-Xg[j]) / epsilon)
            depotx[i,0] = phi((X[0,k]-Xg[i]) / epsilon)
        for j in range(max(jk-qy,0), min(jk+qy+1,my)):
            #depoty[i,0] = phi((X[1,k] -Yg[i]) / epsilon) 
            depoty[0,j] = phi((X[1,k] -Yg[j]) / epsilon) 
        Matricedepot =   depotx * depoty
        R += M[k] / (epsilon * epsilon) * Matricedepot
    return  R
     
