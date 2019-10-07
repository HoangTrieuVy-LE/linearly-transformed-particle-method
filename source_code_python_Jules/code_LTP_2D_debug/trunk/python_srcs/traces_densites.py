#  !/usr/bin/env python
#  - * - codin: utf-8 -*-
from __future__ import division
import numpy  #as np
import matplotlib.pyplot as plt


from calculs_2D import *
from graphes2D import *
from reconstruction_densites_2D import *
from remapping_2D import *

Nmesh = 2000
X0 = numpy.genfromtxt("out/X0.txt",delimiter=",") 
X1 = numpy.genfromtxt("out/X1.txt",delimiter=",") 
M = numpy.genfromtxt("out/M.txt",delimiter=",") 
D0 = numpy.genfromtxt("out/D0.txt",delimiter=",") 
D1 = numpy.genfromtxt("out/D1.txt",delimiter=",") 
D2 = numpy.genfromtxt("out/D2.txt",delimiter=",") 
D3 = numpy.genfromtxt("out/D3.txt",delimiter=",") 

N = len(M)
D = numpy.zeros([4, N]) #matrice de deformation
D[0,:] = D0  
D[1,:] = D1  
D[2,:] = D2  
D[3,:] = D3  
X = numpy.zeros([2, N])
t, npic = numpy.genfromtxt("out/data.txt",delimiter=",") 
X[0, :] = X0
X[1, :] = X1

Sx = numpy.array([-1., 1.])
Sy = numpy.array([-1., 1.])
Ix, Iy, modif_grille = MAJ_grille(X,Sx,Sy)

[Xgrille, Ygrille] = make_grid_unif(Ix, Iy, Nmesh)
R_ltp = rec_densite_grille_unif_LTP(X, M, D, Xgrille, Ygrille)
nom = 'N='+str(N)+'-Nmesh='+str(Nmesh) 
fait_series_dessins(Xgrille, Ygrille, R_ltp, npic, t, nom )


