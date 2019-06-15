# ###########################
#  !/usr/bin/env python
#  - * - codin: utf-8 -*-
from __future__ import division
import numpy  #as np
import scipy  #as sp
import pylab  #as plt
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import integrate
import math  # ou : from math import *
import time
from pylab import *
from scipy import *
from scipy import integrate
from mpl_toolkits.mplot3d import axes3d 
from scipy.interpolate import *
import os
from matplotlib.colors import LightSource
#from matplotlib import anima

#~ from data_code_LTP_2D import *
import config

#~ from calculs_2D import *
#~ from shapefunction2D import b3, b3spline2d
#~ from reconstruction_densites_2D import *

def traceparticules(X,t,npic,show=False):
	figure(0)
	clf()
	plt.axis([config.Ix[0], config.Ix[1], config.Iy[0], config.Iy[1]])
	plt.plot(X[0, :], X[1, :], 'or', ms=2, mfc='k')
	plt.title('Positions of particles at time t=' + str(t) + '.')
	if config.animgif == 1:
		filename = 'out/fichierTemp' + str('%02d' % npic) + '.pdf'
		savefig(filename)
	if config.position_particles == 1:
		#~ print "trace particles"
		filename = 'out/particules' +str('%04d'%npic) + '.png'
		savefig(filename)
	#~ print 'Trace'
	if (show):
		plt.show()		
	else:
		plt.close(0)
    
def traceparticules_avec_vitesse(X, t, U, npic):
    figure(0)
    clf()
    plt.axis([Ix[0], Ix[1], Iy[0], Iy[1]])
    plt.plot(X[0, :], X[1, :], 'or', ms=2, mfc='k')
    plt.title('Positions of particles at time t=' + str(t) + '.')
    X1, Y1= meshgrid(X[0,:], X[1,:])
    DX1, DY1 =  U[0,:], U[1,:]
    M= (hypot(DX1,DY1))
    M [ M== 0] =1.
    DX1 /= M
    DY1 /= M
    plt.quiver (X1, Y1 ,DX1,DY1,M)
    
    if animgif == 1:
        filename = 'out/fichierTemp' + str('%02d' % npic) + '.pdf'
        savefig(filename)
    elif position_particles == 1:
        filename = 'out/particules' +str('%04d'%npic) + '.png'
        savefig(filename)
    print ('Trace')
    plt.close(0)


def faitlefilm(nom):
    if animgif == 1:
        print ('\nanim gif')
        if os.path.exists('./out/Animation.gif') :
            choice=raw_input("file Animation.gif already exists, do you want to erase it ?[y,n]")
            if ( (choice == "y")  or (choice == "yes") ) : 
                print ("computing Animation.gif\n")
                cmd = 'convert ./out/fichierTemp*.pdf ./out/Animation.gif'
                os.system(cmd)
                os.system('rm ./out*.pdf')
            else :
                print ("Animation.gif is not changed\n")
        else :
            print ("computing Animation.gif\n")
            cmd = 'convert ./out/fichierTemp*.pdf ./out/Animation.gif'
            os.system(cmd)
            os.system('rm ./out*.pdf')
            
    if animmp4 == 1:
        print ('anim mp4')
        #   os.system('ffmpeg -r 10 -b 1800 -i %04d.jpg movie.mp4')
        os.system('avconv -r 10 -b 1800 -i out/'+nom+'%04d.png movie.mp4')
        #   os.system('ffmpeg  -i out/%04d.jpg  -r 10 movie.mp4')
        #   os.system('ffmpeg -framerate 1/5 -i %04d.jpg  -r 10 out.mp4')
        #   dans l'exemple precedent chaque image dure 5 secondes
        #   os.system('rm *.jpg')
        
#def fait_series_dessins(Xgrille, Ygrille, R, npic, t, label):
def rep1d_plot(x,y_indic,M,npic,t,label='Density_1D') :
	#X, Y = np.meshgrid(x, y)   
	fig= plt.figure(7)    
	ax = fig.add_subplot(1, 1, 1)
	plt.plot(x,M[:,y_indic])
	#fig.colorbar(surf, shrink=0.5, aspect=10)
	plt.title(label+' at time t='+str(t))
	filename = 'out/'+label+'_'+str('%04d'%npic) +'.png'
	savefig(filename)
	plt.close(7)

def dessin3d(x, y, M, npic, t, label='_'):
    # prend un temps de temps
    X, Y = np.meshgrid(x, y)
    #if truc=='fonction':
    #    vect_f = vectorize(M)
    #    Z = vect_f(X, Y)
    #else:
    #    Z = M
        
    zmax = M.max()
    
    fig = plt.figure(1)
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    surf = ax.plot_surface(X, Y, M, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_zlim3d(0, zmax*1.1)
    fig.colorbar(surf, shrink=0.5, aspect=10)
    plt.title('Density at time t='+str(t))
    filename = 'out/3d' +str('%04d'%npic) + label+'.png'
    savefig(filename)
    plt.close(1)

def rep2d_imshow(M, npic, t):
    fig= plt.figure(2)
    ax1 = fig.add_subplot(1, 1, 1)
    cax = ax1.imshow(numpy.transpose(M), cmap=cm.gnuplot) #cm.YlOrRd) #coolwarm)
    plt.xticks([]); plt.yticks([])
    fig.colorbar(cax, shrink=0.5, aspect=10)
    xlabel('Density at time t='+str(t))
    filename = 'out/imshow' +str('%04d'%npic) + '.png'
    savefig(filename)
    plt.close(2)
    
def rep2d_contourf(x, y, M, npic, t, label='Density',show=False):	
    X, Y = np.meshgrid(x, y)    
    fig= plt.figure(3)    
    ax = fig.add_subplot(1, 1, 1)
    surf = plt.contourf(X, Y, M.transpose() )
    fig.colorbar(surf, shrink=0.5, aspect=10)
    plt.title(label+' at time t='+str(t))
    filename = 'out/'+label+'_'+str('%04d'%npic) +'.png'
    savefig(filename)
	
#    if abs(Tmax/2.-t) < 0.00001 :
#        plt.show()
    if (show):
        plt.show()
    else:
        plt.close(3)        
    
def rep3d_proj(x, y, M, npic, t):
    # prend un temps de temps    
    X, Y = np.meshgrid(x, y)    
    fig= plt.figure(4)
    ax = fig.add_subplot(1, 1, 1)
    ax = fig.gca(projection='3d')  
    ax.view_init(elev=90., azim=0)
    surf = ax.plot_surface(X, Y, M, rstride=1, cstride=1, cmap=cm.YlOrRd, linewidth=0, antialiased=False)
    plt.title('Density at time t='+str(t))
    filename = 'out/3d_proj' +str('%04d'%npic) + '.png'
    savefig(filename)
    plt.close(4)
        
def rep3d_avec_contourf_cote(x, y, M, npic, t, label='_' ):
    M = M.transpose()
    X, Y = np.meshgrid(x, y) 
    zmax = M.max()
    
    fig = plt.figure(5)
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, M, rstride=8, cstride=8, alpha=0.3)
    cset = ax.contourf(X, Y, M, zdir='z', offset=0, cmap=cm.coolwarm)
    cset = ax.contourf(X, Y, M, zdir='x', offset=-1.5, cmap=cm.coolwarm)
    cset = ax.contourf(X, Y, M, zdir='y', offset=1.5, cmap=cm.coolwarm)

    ax.set_xlabel('X')
    ax.set_xlim(min(x)-0.1,  max(x)+0.1)
    ax.set_ylabel('Y')
    ax.set_ylim(min(y)-0.1,  max(y)+0.1)
    ax.set_zlabel('Z')
    ax.set_zlim(0, zmax)
    
    plt.title('Density at time t='+str(t))
    filename = 'out/3dcotesf' +str('%04d'%npic) + label +'.png'
    savefig(filename)
    plt.close(5)
    
def rep3d_contour_cote(x, y, M, npic, t, label='_' ):
    M = M.transpose()
    X, Y = np.meshgrid(x, y) 
    zmax = M.max()   
    fig = plt.figure(6)
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    #ax.plot_wireframe(XX, YY, M, rstride=10, cstride=10)
    cset = ax.contour(X, Y, M, zdir='x', offset=-1.5, cmap=cm.coolwarm)
    cset = ax.contour(X, Y, M, zdir='y', offset=1.5, cmap=cm.coolwarm)
    cset = ax.contour(X, Y, M, zdir='z', offset=0, cmap=cm.coolwarm)
    
    plt.title('Density at time t='+str(t))
    filename = 'out/3dcotes' +str('%04d'%npic) + label + '.png'
    savefig(filename)   
    plt.close(6)
    
def rep3d_mixte(x, y, M, npic, t, label):    
    X, Y = np.meshgrid(x, y) 
    
    fig = plt.figure(figsize=plt.figaspect(0.5))
    #---- First subplot
    ax = fig.add_subplot(1, 2, 1)
    surf = plt.contourf(X, Y, M.transpose() )
    fig.colorbar(surf, shrink=0.5, aspect=10) 
                
    #---- Second subplot
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    ax.plot_wireframe(X, Y, M.transpose(), rstride=10, cstride=10)
    cset = ax.contour(X, Y, M.transpose(), zdir='x', offset=-1.5, cmap=cm.coolwarm)
    cset = ax.contour(X, Y, M.transpose(), zdir='y', offset=1.5, cmap=cm.coolwarm)
    cset = ax.contour(X, Y, M.transpose(), zdir='z', offset=0, cmap=cm.coolwarm)    
    plt.title('Density at time t='+str(t))
    filename = 'out/mixtes' +str('%04d'%npic) +'Density'+label+ '.png'
    savefig(filename)
    #plt.close("all")
    
def fait_series_dessins(Xgrille, Ygrille, R, npic, t, label, show=False):
    nom = 'Density_'+label
    #dessin3d(Xgrille, Ygrille, R, npic, t, nom)
    #rep3d_avec_contourf_cote(Xgrille, Ygrille, R, npic, t, nom)
    #rep3d_contour_cote(Xgrille, Ygrille, R, npic, t, nom )
    #rep3d_mixte(Xgrille, Ygrille, R, npic, t, nom )
    rep2d_contourf(Xgrille, Ygrille, R, npic, t, nom, show)
    #plt.close("all")  
    
def dessine_forme_particule(D, M, X, h, liste_indice_part, npic):
    #phi = b3spline2d
 
    for k in liste_indice_part:
        #phin_vect = vectorize(phin)
        dx = 0.01
        Xg = np.arange(X[0,k]-0.5, X[0,k]+0.5, dx)
        Yg = np.arange(X[1,k]-0.5, X[1,k]+0.5, dx)
        
        # Dk_0 = D[:, k] #.reshape((4, 1))
        # print 'Dk_0', Dk_0
        Dk = D[:, k].reshape((4, 1))
        #print 'Dk', Dk
        Xk = X[:,k].reshape((2,1))
        Mk = numpy.array([1])
                
        Z = rec_densite_grille_unif_LTP(Xk, Mk, Dk, Xg, Yg)
        
        XX, YY = np.meshgrid(Xg, Yg)       
        fig= plt.figure(8)    
        ax = fig.add_subplot(1, 1, 1)        
        
        surf = plt.contourf(XX, YY, Z.transpose())
        fig.colorbar(surf, shrink=0.5, aspect=10)        
        
        #fig = plt.figure()
        #ax = fig.gca(projection='3d')
        #ax.view_init(elev=90., azim=0)
        #surf = ax.plot_surface(XX, YY, Z, rstride=1, cstride=1, cmap=cm.YlOrRd, linewidth=0, antialiased=False)
        #ax.set_zlim(-1.01, 1.01)
        #ax.zaxis.set_major_locator(LinearLocator(10)) # ??
        #ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f')) # ??
        #fig.colorbar(surf, shrink=0.5, aspect=5)
        filename = 'out/forme_part'+str(k)+'_'+str('%04d'%npic) + '.jpg'
        savefig(filename)
        plt.close(8)

def draw_analytical_solution(x,y,sol, npic, t, label='Solution') : 
    X, Y = np.meshgrid(x, y)    
    fig= plt.figure(10)    
    ax = fig.add_subplot(1, 1, 1)
    surf = plt.contourf(X, Y, sol.transpose() )
    fig.colorbar(surf, shrink=0.5, aspect=10)
    plt.title(label+' at time t='+str(t))
    filename = 'out/'+label+'_'+str('%04d'%npic) +'.png'
    savefig(filename)
    plt.close(10)
	#    plt.show(True)
    

def draw_U(X, U,  t, npic, N, label = "U_norm") :
    U_norm=numpy.zeros([N,1])
    for k in range(0,N) :
        U_norm[k] = numpy.sqrt( U[0,k]**2 + U[1,k]**2 ) 
    fig= plt.figure(10)    
    #~ ax = fig.add_subplot(1, 1, 1)	
    #~ for k in range(0,N) :
        #~ surf = plt.contourf(X[0,k], X[1,k], U_norm[k])    
        #~ fig.colorbar(surf, shrink=0.5, aspect=10)
      #~ 
    cm = plt.cm.get_cmap('RdYlBu')
    #~ xy = range(20)
    #~ z = xy
    sc = plt.scatter(X[0,:], X[1,:], c=U_norm[:], vmin=0, vmax=20, s=35, cmap=cm)
    plt.colorbar(sc)    
    plt.title(label+' at time t='+str(t))
    filename = 'out/'+label+'_'+str('%04d'%npic) +'.png'
    savefig(filename)
    #~ plt.show()
    plt.close(10)

