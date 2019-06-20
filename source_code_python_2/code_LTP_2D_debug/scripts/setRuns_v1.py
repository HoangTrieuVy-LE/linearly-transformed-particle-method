import os
import sys
import numpy as np

inputs_dir = sys.argv[1]

if (not(os.path.isdir(str(inputs_dir)))) :
    os.makedirs(str(inputs_dir))

nb_methods = len(sys.argv[2:])
methods_list  = sys.argv[2:]

if (methods_list == []) :
	raise ValueError('wrong number of arguments = '+str(int(nb_methods+1))+'. \n1 or 2 arguments after inputs_dir argument : first method name and second method name (LTP or SP)')

### TO BE DONE IN THE END
#Nx_arr=np.array([50, 75, 100, 150, 200, 300, 500, 1000])
#dt_arr=np.array([0.01, 0.005, 0.001, 0.0005, 0.0001])

Nx_arr=np.array([50, 100, 150, 200, 250, 500])
dt_arr=np.array([0.001, 0.0001])

delta_no_CFL_arr = np.array([-20, -15, -10, -5, -1., -0.5, -0.1, -0.05, -0.01, 0.,  0.01, 0.05, 0.1, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 7.5, 10., 15, 20.])
q_arr =  np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.])

count=0.


for method in methods_list :

    #______________________________________________________________________#
    #
    #   method=LTP , D_method=hybrid
    #______________________________________________________________________#
    
    if (method == 'LTP') : 
    	print "TOTAL NUMBER OF RUNS FOR LTP, D_method : hybrid = " , len(Nx_arr)*len(dt_arr)*len(delta_no_CFL_arr)
        count += len(Nx_arr) *len(dt_arr)*len(delta_no_CFL_arr)
        i=0.
        for Nx in Nx_arr : 
            for dt in dt_arr :
                h = 4./Nx
                CFL = dt / h
                delta_arr = CFL * delta_no_CFL_arr
                epsilon_arr=(1+delta_arr)*h
                for epsilon in epsilon_arr :            
                    fileName=str(inputs_dir)+'/data_LTP_'+str('%04d'%i)+'.csv'
                    os.system("touch "+str(fileName))
                    myInFile=open(fileName,'w')
                    myInFile.write('Lx1=-2.\n')
                    myInFile.write('Lx2=2.\n')
                    myInFile.write('Ly1=-2.\n')
                    myInFile.write('Ly2=2.\n')
                    myInFile.write('Tini=0.05\n')
                    myInFile.write('time_scheme=euler_explicit\n')
                    myInFile.write('m=2\n')
                    myInFile.write('d=2\n')
                    myInFile.write('C=1.\n')
                    myInFile.write('indic_remapping=no\n')
                    myInFile.write('l_remap=6\n')
                    myInFile.write('radius_remap=2\n')
                    myInFile.write('met_calcul_poids=quasi_interp\n')
                    myInFile.write('part_sur_bords=non\n')
                    myInFile.write('bool_phi_radial=1\n')
                    myInFile.write('deg_spline=3\n')
                    myInFile.write('explosion=no\n')
                    myInFile.write('Tmax=4.\n')
                    myInFile.write('method=LTP\n')
                    myInFile.write('D_method=hybrid\n')
                    myInFile.write('Nx='+str(Nx)+'\n')
                    myInFile.write('Ny='+str(Nx)+'\n')
                    myInFile.write('dt='+str(dt)+'\n')
                    myInFile.write('epsilon='+str(epsilon)+'\n')
                    myInFile.close()
                    i+=1 #incremente file number
        
    #______________________________________________________________________#
    #
    #   method=SP , D_method=None
    #______________________________________________________________________#    
    
    if (method == 'SP') : 
        print "TOTAL NUMBER OF RUNS FOR SP, D_method : None = " , len(Nx_arr)*len(dt_arr)*len(q_arr)
        epsilon_arr = h**q_arr
        count += len(Nx_arr) *len(dt_arr)*len(q_arr)
        i=0.
        for Nx in Nx_arr : 
            for dt in dt_arr :
                for epsilon in epsilon_arr :
                    fileName=str(inputs_dir)+'/data_SP_'+str('%04d'%i)+'.csv'
                    os.system("touch "+str(fileName))
                    myInFile=open(fileName,'w')
                    myInFile.write('Lx1=-2.\n')
                    myInFile.write('Lx2=2.\n')
                    myInFile.write('Ly1=-2.\n')
                    myInFile.write('Ly2=2.\n')
                    myInFile.write('Tini=0.05\n')
                    myInFile.write('time_scheme=euler_explicit\n')
                    myInFile.write('m=2\n')
                    myInFile.write('d=2\n')
                    myInFile.write('C=1.\n')
                    myInFile.write('indic_remapping=no\n')
                    myInFile.write('l_remap=6\n')
                    myInFile.write('radius_remap=2\n')
                    myInFile.write('met_calcul_poids=quasi_interp\n')
                    myInFile.write('part_sur_bords=non\n')
                    myInFile.write('bool_phi_radial=1\n')
                    myInFile.write('deg_spline=3\n')
                    myInFile.write('explosion=no\n')
                    myInFile.write('Tmax=4.\n')
                    myInFile.write('method=SP.\n')
                    myInFile.write('D_method=None.\n')
                    myInFile.write('Nx='+str(Nx)+'.\n')
                    myInFile.write('Ny='+str(Nx)+'.\n')
                    myInFile.write('dt='+str(dt)+'.\n')
                    myInFile.write('epsilon='+str(epsilon)+'\n')
                    myInFile.close()
                    i+=1 #incremente file number
                    
                    
print "TOTAL NUMBER OF FILES CREATED = " , count , "\n\n"
