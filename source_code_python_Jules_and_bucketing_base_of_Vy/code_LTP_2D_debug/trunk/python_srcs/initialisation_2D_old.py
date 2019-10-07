#!/usr/bin/env python
# - * - codin: utf-8 -*-
from __future__ import division
import numpy  # as np
import scipy  # as sp
import pylab   # as plt
import matplotlib  
import matplotlib.pyplot as plt
from scipy import integrate
import math  # ou : from math import *
import time
from pylab import *
from scipy import *
import os
from random import *

#from _calculsfor_f90 import  calculsfor_ini_modf90

#  ============= Possible initial data ================================
# rem : a mettre dans le fichier data_code_LTP_2D ?       
k = 0.5
Lx1 = -1 
Lx2 = 1  
Ly1 = -1  
Ly2 = 1  

def rhoini1(x, y):
    if Lx1<x<Lx2:
        ax =1+0.01*cos(k*x)
    else:
        ax =0.
    #ax = math.exp(-30 * (x - 0.5) * (x - 0.5)) + 2 * math.exp(-50 * (x + 0.3) * (x + 0.3))
    if Ly1<y<Ly2:
        ay = math.exp(-y*y/2)/sqrt(2*3.14)
    else:       
        ay = 0.
    return ax * ay
    
def rhoini2(x, y, r=1.):
    n =sqrt(x*x+y*y)
    if n<=r:
        return 1.
    else:
        return 0.       
            
    
    
#  ============= Functions for Initialization ==========================

#  Adimensionnement de la densite initiale
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

#  Position initiales des particules
def init_poids_et_position(f, sx, sy, nx, ny):
    initialisation = 'deterministe'
    if initialisation == 'aleatoire':
        if method == 'part_Mrand':
            return initialize_part_Mrand(sx, sy, nx, ny)
        elif method == 'part_Prand':
            return initialize_part_Prand(sx, sy, nx, ny)        
    else:
        #return initializepartssbords(f, sx, sy, nx, ny,'formulespline')
        return initializepartssbords(f, sx, sy, nx, ny,'quadrature')

def initializepartssbords(f, sx, sy, nx, ny, met):  #pas de particules sur les bords du domaine
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
            # initialisation avec Formule Martin
            if met == 'formulespline':
                t0 = f(x, y) * (8/6)**2
                t01 = -8/36 *(f(x-hx, y)+f(x+hx, y)+f(x, y+hy)+f(x+hx, y-hy))
                t11 =1/36* (f(x-hx, y-hy)+f(x+hx, y+hy))
                #z[2, k] = hx * hy *(t0+t01+t11)  
                z[2, k] = hx * hy * max(t0+t01+t11,0) 
            elif met =='quadrature':
                z[2, k] = int2d(f,x-hx/2.,x+hx/2,y-hy/2.,y+hy/2.,10)
            else:
                z[2, k] = hx * hy * f(x,y)
            k += 1
    return z

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
            if (i > 0) and (i < Nx) and (j > 0) and (j < Nx):
                z[2, k] = f(x, y) * hx * hy
            elif ((i == 1) or (i == Nx)) and ((j == 1) or (j == Nx)):
                z[2, k] = f(x, y) * hx * hy * 7. / 4.
            else:
                if methode=='ponctuelle':
                    print  'initialisation des poids ponctuelle'
                    z[2, k] == f(x, y) * hx * hy * 3. / 4.             
            # initialisation des poids avec integrale
                else:
                    print 'initialisation des poids avec integrale'
                    z[2, k] = int2d(f,x-hx/2.,x+hx/2,y-hy/2.,y+hy/2.,100)
            k += 1
    return z

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
    
def initialize_part_Prand(sx, sy, nx, ny):
    z = zeros([3, nx * ny])
    k = 0
    for i in range(ny):
        for j in range(nx):
            z[0, k] = sx[0]+random()*(sx[1] - sx[0]) 
            z[1, k] = sy[0]+random()*(sy[1] - sy[0]) 
            z[2, k] = 1./(nx*ny)
            k +=1
    return z
    
def supprime_particule(X,M):
    N = len(M)
    poid_tot = sum(M)
    tol = poid_tot/N*0.01
    #print 'tol', tol
    Nnew = N
    #print 'M avant suppression particules', M
    #print 'X avant suppression particules', X
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
    #if knew != Nnew:
    #    print 'probleme'
    #print 'M apres suppression particules', Mnew
    #print 'X apres suppression particules', Xnew
    return Nnew, Xnew, Mnew

#def remapping(R,Xg, Yg):
    
