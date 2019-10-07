import numpy

import config
import calculs_2D
from _calculsfor_f90 import  calculsfor_var_modf90
#~ from data_code_LTP_2D import *
						
#______________________________________________________________________#
#
# 				THEORETICAL FUNCTIONS OF DIFFUSION SOLUTION
#______________________________________________________________________#


def diffusion_init(x,y) :	 
    """
    Theoretical profile of initial function for diffusion problem
    if function is axisymetric then the global solution will be its convolution
    by heat kernel
    """       
    #~ K = (1 / (4 * numpy.pi * config.Tini)**2) * numpy.exp( - (x**2 + y**2) / (4*config.Tini) )
    
    
    init  = (1 / (4 * numpy.pi * config.Tini)) * numpy.exp( - (x **2 +y**2) / (4*config.Tini) )
    
    return init

def heat_kernel(X,Y,t) :	 
    """
    Theoretical heat kernel expression
    """       
    
    K = (1 / (4 * numpy.pi * t)) * numpy.exp( - ((X[0] - Y[0])**2 + (X[1] - Y[1])**2) / (4*t) )

    return K

def diffusion_sol(x,y,t) :	 
    """
    Analytical solution to the global diffusion problem in one point 
    We first compute quadrature points for trapeze method
    then we compute convolution on (x,y) coordinate we test function config.phi_2d
    """       
    quadrature_pts , Xgrid_convol, Ygrid_convol = calculs_2D.get_quadrature_pts_convol()
    nqx = quadrature_pts.shape[0]
    nqy = quadrature_pts.shape[1]
    # warning dx and dy might slightly different  (upt to 0.1% difference)
    # than prescribed previously
    dx = abs(Xgrid_convol[1] - Xgrid_convol[0])
    dy = abs(Ygrid_convol[1] - Ygrid_convol[0])
    #~ print "dx , dy  = ", dx , dy
    Y = numpy.zeros(2)
    # get current position 
    X = numpy.zeros(2)
    X[0] = x
    X[1] = y
    sol = 0.
    for i in range(0,nqx) :
        for j in range(0,nqy) :
            # get grid point coordinate
            #~ print Xgrid_convol[i] , Xgrid_convol[j] , quadrature_pts[i,j]
            Y[0] = Xgrid_convol[i]
            Y[1] = Ygrid_convol[j]
            # compute analytic solution on point (x,y)
            # that's a convolution               
            sol = sol + heat_kernel(X,Y,t) * diffusion_init(Y[0], Y[1]) * (1/4.) * dx * dy * quadrature_pts[i,j]  
    return sol 
    

def diffusion_sol_analytic(x,y,t) :	 
    """
    Only if the initial density is a gaussian : because convolution of gaussian is a gaussian.
    """       
    
    sol = 1/(4 * numpy.pi *(config.Tini + t)) * numpy.exp( - (x**2 + y**2) / (4*(t+config.Tini)) )
    #~ sol = 1/(2 * numpy.pi *(config.Tini + t)) * numpy.exp( - (x**2 + y**2) / (2*(t+config.Tini)) )
    
    return sol
