import numpy

import config
#~ from data_code_LTP_2D import *
						
#______________________________________________________________________#
#
# 				THEORETICAL BARENBLATT DENSITIES
#______________________________________________________________________#




## order m of PME equation $\partial _t f - \nabla \cdot (mf^{m-1} \nabla f) = 0$

# necessary to create a function barenblatt_init that takes only 2 arguments because if we consider that
# rhoini(x,y) = barenblatt(x,y,t=config.Tini) like before, config.Tini is initialized with its value 
# from config.py and NOT from data_code_LTP_2D.initialize_data_parameters
# Not sure why but it probably has to do with the way python "compiles" 
# Maybe change that in a later perspective... (not really a big deal)
def barenblatt_init(x,y) :	 
    """
    Theoretical profile of initial Barenblatt (t = config.Tini) density
    """       
    
    if (config.C - config.k *(x**2 + y**2) * (1/config.Tini**(2*config.beta)) >= 0):
        sol_barenblatt = 1/(config.Tini**config.alpha) * (max( config.C - config.k *(x**2 + y**2) * (1./config.Tini**(2.*config.beta)) , 0))**(1/(config.m_barenblatt-1))
    else :
        sol_barenblatt = 0.
    return sol_barenblatt


def barenblatt(x,y,t) :	 
    """
    Theoretical profile of Barenblatt density
    """       
    if (config.C - config.k *(x**2 + y**2) * (1/t**(2*config.beta)) >= 0):
        sol_barenblatt = 1/(t**config.alpha) * (max( config.C - config.k *(x**2 + y**2) * (1/t**(2*config.beta)) , 0))**(1/(config.m_barenblatt-1))		    
    else :
        sol_barenblatt = 0.
        
    return sol_barenblatt

def grad_barenblatt(x,y,t) :	
    """
	Theoretical profile of gradient of Barenblatt density
    """  
    grad_sol_barenblatt=numpy.zeros(2) 
    if (config.C - config.k *(x**2 + y**2) * (1/t**(2*config.beta)) >= 0):
        if (config.m_barenblatt == 2) :
            grad_sol_barenblatt[0] = 1/(t**config.alpha) * (-2*config.k*x*(1/t**(2*config.beta))) 
            grad_sol_barenblatt[1] = 1/(t**config.alpha) * (-2*config.k*y*(1/t**(2*config.beta))) 
        else :
            grad_sol_barenblatt[0] = 1/(t**config.alpha) * (-2*k*x*(1/t**(2*config.beta))) * ( C - k *(x**2 + y**2) * (1/t**(2*config.beta)))**(1/(config.m_barenblatt-1) - 1)
            grad_sol_barenblatt[1] = 1/(t**config.alpha) * (-2*k*y*(1/t**(2*config.beta))) * ( C - k *(x**2 + y**2) * (1/t**(2*config.beta)))**(1/(config.m_barenblatt-1) - 1)
    else : 
        grad_sol_barenblatt[:] = 0.
    return grad_sol_barenblatt

def hess_barenblatt(x,y,t) :	
    """
	Theoretical profile of Hessian Barenblatt density
    """  
    hess_sol_barenblatt=numpy.zeros([2,2])
    
    if (config.C - config.k *(x**2 + y**2) * (1/t**(2*config.beta)) >= 0):
        if (config.m_barenblatt == 2) :
            hess_sol_barenblatt[0,0] = 1/(t**config.alpha) * (-2*config.k*(1/t**(2*config.beta)))
            hess_sol_barenblatt[0,1] = 0. 
            hess_sol_barenblatt[1,0] = 0.
            hess_sol_barenblatt[1,1] = 1/(t**config.alpha) * (-2*config.k*(1/t**(2*config.beta))) 
        else :
			# A MODIFIER POUR CAS GENERAL
            hess_sol_barenblatt[0,0] = 1/(t**config.alpha) * (-2*config.k*(1/t**(2*config.beta))) 
            hess_sol_barenblatt[0,1] = 0. 
            hess_sol_barenblatt[1,0] = 0.
            hess_sol_barenblatt[0,0] = 1/(t**config.alpha) * (-2*config.k*(1/t**(2*config.beta))) 
    else : 
        hess_sol_barenblatt[:,:] = 0.
    return hess_sol_barenblatt

# compute F^{s,t} (x,y) = (t/s)^1/4  * X_(s) (e.g. for m=d=2)
def flow_barenblatt(Xs, s, t) :
    X_exact = numpy.zeros(2)
    if (config.m_barenblatt == 2) :        
        if ( ( Xs[0]**2 + Xs[1]**2) <= 16 * s**0.5 ) :
            X_exact = (t/s)**(1./4.) * Xs            
        else :  # i.e. : ( ( Xs[0]**2 + Xs[1]**2) > 16 * numpy.sqrt(s) )
            if 16 * numpy.sqrt(t) <= ( Xs[0]**2 + Xs[1]**2) :
                X_exact[:] = Xs
            else :
                X_exact = -10000.
                #~ raise ValueError("no solution for values of s, t, Xs : "+str(s)+", "+str(t)+", "+str(Xs))
    else :
        raise ValueError("m_barenblatt != 2 not implemented yet. Current m_barenblatt value : "+str(config.m_barenblatt))
    return X_exact

#~ def inv_flow_barenblatt(Xs, s, t) :
    #~ X_exact = numpy.zeros(2)
    #~ if (config.m_barenblatt == 2) :                
        #~ if ((Xs[0]**2 + Xs[1]**2) <= 16 * numpy.sqrt(t))   :
                #~ X_exact = (s/t)**(1./4.) * Xs
        #~ if ((Xs[0]**2 + Xs[1]**2) >= 16 * numpy.sqrt(s))   :                        
            #~ X_exact[:] = Xs
    #~ else :
        #~ raise ValueError("m_barenblatt != 2 not implemented yet. Current m_barenblatt value : "+str(config.m_barenblatt))
    #~ return X_exact
#~ 
def jacobian_flow_barenblatt(Xs, s, t) :
    J_st = numpy.zeros([2,2])
    Id = numpy.zeros([2,2])
    Id[0,0] = 1.
    Id[1,1] = 1.
    if (config.m_barenblatt == 2) :                
        if not( (abs(Xs[0] - (-10000) ) <= 0.01) and (abs(Xs[1] - (-10000) ) <= 0.01) ) :        
            if ( ( Xs[0]**2 + Xs[1]**2) <= 16 * s**0.5 ) :
                J_st = (t/s)**(1./4.) * Id            
            else :  # i.e. : ( ( Xs[0]**2 + Xs[1]**2) > 16 * numpy.sqrt(s) )
                if (16 * numpy.sqrt(t) <= ( Xs[0]**2 + Xs[1]**2)) :
                    J_st[:] = Id
                else :
                    J_st[0,1] = -10000.
                    J_st[1,0] = -10000.
                    #~ raise ValueError("no solution for values of s, t, Xs : "+str(s)+", "+str(t)+", "+str(Xs))        if (config.C - config.k *(Xs[0]**2 + Xs[1]**2) * (1/s**(2*config.beta)) >= 0):            
        else :            
            J_st = -10000. * Id
            J_st[0,1] = -10000.
            J_st[1,0] = -10000.            
    else :
        raise ValueError("m_barenblatt != 2 not implemented yet. Current m_barenblatt value : "+str(config.m_barenblatt))
    return J_st

#!!!! OLD VERSION	
def Dk_exact_barenblatt(Xini, t) :
    Dk = numpy.zeros([2,2])
    Id = numpy.zeros([2,2])
    Id[0,0] = 1.
    Id[1,1] = 1.
    if (config.m_barenblatt == 2) :        
        #~ if (Xini[0]**2 + Xini[1]**2 < 16*numpy.sqrt(t)) :
        if (config.C - config.k *(Xini[0]**2 + Xini[1]**2) * (1/t**(2*config.beta)) >= 0):
            Dk = (config.Tini/t)**(1./4.) * Id
        else :
            #~ Dk[:] = 0.
            Dk[:] = Id
    else :
        raise ValueError("m_barenblatt != 2 not implemented yet. Current m_barenblatt value : "+str(config.m_barenblatt))
    return Dk
