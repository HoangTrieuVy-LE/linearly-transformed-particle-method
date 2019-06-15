import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
import pylab

import numpy as np
from scipy import stats

if (str(sys.argv[1]) != '-SeveralErrors') :
    print "first param : -SeveralErrors\nsecond param : option (dt, epsilon or N)\nthird param : outfile.ext (ext = png or svg)\nother params : errors*.txt"

### Errors for several files 
# arg 1 = -SeveralErrors 
# arg 2 = option (error against N, dt, or epsilon
# arg > 2 = all files names
elif (str(sys.argv[1]) == '-SeveralErrors') :		
    
    params = {'legend.fontsize': 40,
              'legend.linewidth': 2}
    matplotlib.rcParams.update({'font.size': 50})
    pylab.rcParams.update(params)
    plt.rc('text', usetex=True)
    linestyles = ['-.', '-', '--', ':',' ']    
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)] 
             
    for i in range(len(tableau20)):  
        r, g, b = tableau20[i]  
        tableau20[i] = (r / 255., g / 255., b / 255.) 

    option=sys.argv[2]
    outFileName = sys.argv[3]
    nb_of_files = len(sys.argv[4:])
    print "\n Number of files from which errors will be computed = " , nb_of_files , "\n"
    print "option choosen : " , option
    print "output file name  : " , outFileName
    
    extension = os.path.splitext(outFileName)[1] # either .eps or .png	or 'svg'
    
    label_list = [os.path.splitext(f)[0] for f in sys.argv[4:]]
    label_list = [label.split('/')[-1].split('.')[0].replace("_"," ") for label in label_list]
    print label_list
    #~ label_list = ["$"+label.split('_')[0]+"_{"+label.split('_')[1]+"}$" for label in label_list]
    
    plt.figure(10)
    if (option == 'epsilon') :					
        for i in range(0,int(nb_of_files)) :			
            inFileName=sys.argv[4+i]
            print "input file name = " , inFileName
            
            data=np.genfromtxt(fname=inFileName, names = True)			
            epsilon_list = data['epsilon']
            error_list = data['Linfinity_error']
            
            #~ delta_list = 
            print "all errors : " , error_list
            
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('$\epsilon$')
            plt.ylabel('$L^{\infty}$ error')
            plt.plot(epsilon_list , error_list, label = '$'+str(label_list[i])+'$', linewidth = 2)
            plt.legend(loc=1)
    
    if (option == 'N') :			
        plt.legend(loc=2)
        for i in range(0,int(nb_of_files)) :
            inFileName=sys.argv[3+i]
            data=np.genfromtxt(fname=inFileName, names = True)			
            N_list = data['N']
            error_list = data['Linfinity_error']
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('$N$')
            plt.ylabel('$L^{\infty}$ error')
            plt.plot(N_list , error_list, label = '$'+str(label_list[i])+'$',  linewidth = 2)

    if (option == 'dt') :
        plt.legend(loc=2)
        for i in range(0,int(nb_of_files)) :
            ls = linestyles[i]
            inFileName=sys.argv[4+i]
            data=np.genfromtxt(fname=inFileName, names = True)
            dt_list = data['dt']
            error_X_list = data['Linfinity_error_X']
            error_D_list = data['Linfinity_error_D']
            indices_sorted = np.argsort(dt_list)
            dt_list_sorted = dt_list[indices_sorted]          #[9:len(dt_list)]
            error_X_list_sorted = error_X_list[indices_sorted]#[9:len(dt_list)]
            error_D_list_sorted = error_D_list[indices_sorted]#[9:len(dt_list)]
            #~ print dt_list_sorted
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('$\Delta t$')
            plt.ylabel('$L^{\infty}$ error')
            
            slope_X, intercept, r_value, p_value, std_err = stats.linregress(np.log(dt_list_sorted), np.log(error_X_list_sorted))
            slope_D, intercept, r_value, p_value, std_err = stats.linregress(np.log(dt_list_sorted), np.log(error_D_list_sorted))            
            print 'slope(dt ,',label_list[i] ,' X) = ' , slope_X
            print 'slope(dt ,',label_list[i] ,' D) = ' , slope_D
            plt.plot(dt_list_sorted , error_X_list_sorted , label = '$\mbox{'+str(label_list[i])+' X}$',  linestyle = linestyles[0], linewidth = 2, color = tableau20[2*i], marker ='x', markersize = 10)
            plt.plot(dt_list_sorted , error_D_list_sorted , label = '$\mbox{'+str(label_list[i])+' D}$',  linestyle = linestyles[1], linewidth = 2, color = tableau20[2*i], marker ='x', markersize = 10)
            
            #~ plt.plot(dt_list_sorted , error_X_list_sorted / dt_list_sorted, label = '$\mbox{'+str(label_list[i])+' X}$',  linestyle = linestyles[0], linewidth = 2, color = tableau20[i+3], marker ='x', markersize = 10)
            #~ plt.plot(dt_list_sorted , error_D_list_sorted / dt_list_sorted, label = '$\mbox{'+str(label_list[i])+' D}$',  linestyle = linestyles[1], linewidth = 2, color = tableau20[i+3], marker ='x', markersize = 10)

    plt.legend(loc=4)
    figure = plt.gcf() # get current figure
    figure.set_size_inches(20, 15)
	#~ plt.savefig(outFileName, format=str(extension[1:]), dpi=100)
    plt.savefig(outFileName, format=str(extension[1:]), dpi=100)		
    plt.show()
