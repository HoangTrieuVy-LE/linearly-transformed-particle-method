import os
import sys
import numpy as np


parser_file = sys.argv[1]
inputs_dir = sys.argv[2]
method = sys.argv[3]
if (not(os.path.isdir(str(inputs_dir)))) :
    os.makedirs(str(inputs_dir))

#~ nb_methods = len(sys.argv[3:])
#~ methods_list  = sys.argv[3:]
#~ 
#~ if (methods_list == []) :
	#~ raise ValueError('wrong number of arguments = '+str(int(nb_methods+1))+'. \n1 or 2 arguments after inputs_dir argument : first method name and second method name (LTP or SP)')


### TO BE DONE IN THE END
#Nx_arr=np.array([50, 75, 100, 150, 200, 300, 500, 1000])
#dt_arr=np.array([0.01, 0.005, 0.001, 0.0005, 0.0001])


parser = np.genfromtxt(fname=str(parser_file), names=True)
list_Nx=parser['Nx']
list_Ny=parser['Ny']
list_dt=parser['dt']
list_epsilon=parser['epsilon']
list_parameters = []
for i in range(0,len(list_Nx)) :
	dic_parameters = {}
	dic_parameters['Nx'] = int(list_Nx[i])
	dic_parameters['Ny'] = int(list_Ny[i])
	dic_parameters['dt'] = list_dt[i]
	dic_parameters['epsilon'] = list_epsilon[i]
	list_parameters.append(dic_parameters)

count=0.
print list_parameters

i=0.
for params in list_parameters :            
    Nx = params['Nx']
    Ny = params['Ny']
    dt = params['dt']        
    epsilon = params['epsilon']
    print Nx , Ny ,  dt , epsilon
    fileName=str(inputs_dir)+'/data_LTP_'+str('%04d'%i)+'.csv'
    os.system("touch "+str(fileName))
    myInFile=open(fileName,'w')
    myInFile.write('Lx1=-2.\n')
    myInFile.write('Lx2=2.\n')
    myInFile.write('Ly1=-2.\n')
    myInFile.write('Ly2=2.\n')
    myInFile.write('Tini=0.05\n')
    myInFile.write('Tmax='+str(dt+0.05)+'\n')
    myInFile.write('dt='+str(dt)+'\n')
    myInFile.write('time_scheme=middle_point\n')
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
    myInFile.write('method='+str(method)+'\n')
    myInFile.write('D_method=hybrid\n')
    myInFile.write('Nx='+str(Nx)+'\n')
    myInFile.write('Ny='+str(Ny)+'\n')            
    myInFile.write('epsilon='+str(epsilon)+'\n')
    myInFile.close()
    i+=1 #incremente file number
    
    
