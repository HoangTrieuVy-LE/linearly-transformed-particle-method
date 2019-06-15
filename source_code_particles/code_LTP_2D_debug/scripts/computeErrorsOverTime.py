#!/usr/bin/env python
# ###########################################################################################

import csv
import os
import sys
import glob
import numpy  as np
import matplotlib.pyplot as plt

sys.path.append('../../trunk/python_srcs/')
sys.path.append('../../trunk/fortran_srcs/')

from calculs_2D import *
from graphes2D import *
from reconstruction_densites_2D import *
from remapping_2D import *
from _calculsfor_f90 import  calculsfor_var_modf90
from _calculsfor_f90 import  calculsfor_rec_modf90

#modele is ~/path/to/dir/inputs/data_LTP_0001.out
#then we look for data files at every time steps that is in directory
input_file_name = sys.argv[1]
out_file_name = sys.argv[2]

inputs_files = glob.glob(str(input_file_name)+'*.txt')
time_steps=[]
bool_density_reconstruted = 0
list_files_name_R_ltp = []
list_files_name_solution = []
bool_data_txt  = 0
for fname in inputs_files :
    if ("solution" in fname) : 
        list_files_name_solution.append(fname)
        bool_density_reconstruted = 1
    if (("R_LTP" in fname) or ("R_ltp" in fname)) :
        list_files_name_R_ltp.append(fname)
        bool_density_reconstruted = 1
    if "data.txt" in fname : 
        bool_data_txt = 1
        with open(fname) as f :
            content = f.readlines()
            content = [x.strip() for x in content]			
            time_steps.append(round(float(content[0]),8))

if (bool_data_txt == 0) :    
    os.system("ls "+str(input_file_name)+"*solution* | cut -d'=' -f3 | cut -d'_' -f1 > tmp")
    with open("tmp") as f :
        time_steps = f.readlines()
        time_steps = [float(x.strip()) for x in time_steps]			

print time_steps
        
extension = os.path.splitext(out_file_name)[1] # either .eps or .png	or 'svg'
path_to_project_dir = "/home/jdichamp/These_MEGEP/Code_Particle_Methods/Code_2D/LTP_2D_convection_diffusion/LTP_2D_barenblatt_OMP/" 
with open(input_file_name) as f :
    content = f.readlines()
    content = [x.strip() for x in content]    
data_name_in_out_file = str(content[1].replace('\t', ',').split(',')[0])

#~ raise ValueError("stop")

data_file_name = str(path_to_project_dir)+str(data_name_in_out_file)

Nmesh = 2000
npic=0
list_errors=[]
list_epsilon=[]
for t in time_steps :    
    if (bool_density_reconstruted == 1) :        
        for l1 in list_files_name_R_ltp  :            
            if (str(t) in str(l1)) :
                file_name_R_ltp = l1
            for l2 in list_files_name_solution  :              
                if (str(t) in str(l2)) :
                    file_name_solution = l2        
        R_ltp = np.genfromtxt(str(file_name_R_ltp), delimiter=" ")         
        solution = np.genfromtxt(str(file_name_solution), delimiter=" ")
    else  :
        file_name = str(input_file_name)+"_t="+str(t)
        print file_name
        X0 = np.genfromtxt(str(file_name)+"_X0.txt",delimiter=",") 
        X1 = np.genfromtxt(str(file_name)+"_X1.txt",delimiter=",") 
        M  = np.genfromtxt(str(file_name)+"_M.txt",delimiter=",") 
        D0 = np.genfromtxt(str(file_name)+"_D0.txt",delimiter=",") 
#  
#  name: inconnu
#  @param
#  @return
#  
        D1 = np.genfromtxt(str(file_name)+"_D1.txt",delimiter=",") 
        D2 = np.genfromtxt(str(file_name)+"_D2.txt",delimiter=",") 
        D3 = np.genfromtxt(str(file_name)+"_D3.txt",delimiter=",") 
        #~ 
        #~ 
        N = len(M)
        D = np.zeros([4, N]) #matrice de deformation
        D[0,:] = D0  
        D[1,:] = D1  
        D[2,:] = D2  
        D[3,:] = D3  
        X = np.zeros([2, N])
        npic +=1
        X[0, :] = X0
        X[1, :] = X1
        
        Sx = np.array([-2., 2.])
        Sy = np.array([-2., 2.])
        Ix, Iy, modif_grille = MAJ_grille(X,Sx,Sy)
        [Xgrille, Ygrille] = make_grid_unif(Ix, Iy, Nmesh)
        NgrilleX = len(Xgrille)
        NgrilleY = len(Ygrille)    
        calculsfor_rec_modf90.alloc_ro(NgrilleX, NgrilleY)
        calculsfor_rec_modf90.set_phi_radial(1)
        calculsfor_rec_modf90.set_degre_spline(3)
        calculsfor_rec_modf90.set_hx_hy(0.08, 0.08)
        calculsfor_rec_modf90.set_grille(Xgrille, Ygrille, NgrilleX, NgrilleY)
        calculsfor_rec_modf90.rec_ltp(X, M, D, N)
        Ro = calculsfor_rec_modf90.ro
        R_ltp = Ro.copy()
        calculsfor_rec_modf90.dealloc_ro()
        nom = 'N='+str(N)+'-Nmesh='+str(Nmesh) 
        #    #~ fait_series_dessins(Xgrille, Ygrille, R_ltp, npic, t, nom )
        solution=compute_solution(Xgrille,Ygrille, t) 
             
    relative_error=error(solution , R_ltp)
    list_errors.append(relative_error)
        
    with open(data_file_name, 'rb') as f:        
        reader = csv.reader(f, delimiter='=')
        for row in reader:
            if (row[0] == 'epsilon') : 		epsilon = float(row[1])
    list_epsilon.append(epsilon)
    
        
print list_epsilon
print list_errors        
        
### PLOT PARAMETERS
params = {'legend.fontsize': 40,
	          'legend.linewidth': 2}
matplotlib.rcParams.update({'font.size': 50})
pylab.rcParams.update(params)
plt.rc('text', usetex=True)

plt.legend(loc=2)
plt.plot(time_steps , list_errors,   linewidth = 2)
    
figure = plt.gcf() # get current figure
figure.set_size_inches(20, 15)
plt.savefig(out_file_name, format=str(extension[1:]), dpi=100)
    
plt.show()

