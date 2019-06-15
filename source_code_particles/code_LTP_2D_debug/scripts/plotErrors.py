import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
import pylab





if (len(sys.argv[1:]) == 3) :	

	inFileName=sys.argv[1]
	outFileName=sys.argv[2] # parameter for when we will automatically save figs
	option=sys.argv[3]
	extension = os.path.splitext(outFileName)[1] # either .eps or .png	
	#load data from inFileName (arg 1)
	data=np.genfromtxt(fname=inFileName, names = True)
	if (option == 'N') :	
		N_list = data['N']
		error_list = data['error']
		print "N List = " , N_list
		print "errors List = " , error_list
		
		params = {'legend.fontsize': 80,
		          'legend.linewidth': 2}
		matplotlib.rcParams.update({'font.size': 50})
		pylab.rcParams.update(params)
		plt.rc('text', usetex=True)
		plt.figure(10)
		plt.plot(N_list , error_list )
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel('$N$')
		plt.ylabel('$L^{\infty}$ error')
		figure = plt.gcf() # get current figure
		figure.set_size_inches(20, 15)		
		plt.savefig(outFileName, format='eps', dpi=1000)
		
	
	elif (option == 'dt') :
		dt_list = data['dt']
		error_list = data['error']
		print "dt List = " , dt_list
		print "dt List = " , error_list
		
		params = {'legend.fontsize': 80,
		          'legend.linewidth': 2}
		matplotlib.rcParams.update({'font.size': 50})
		pylab.rcParams.update(params)
		plt.rc('text', usetex=True)
		plt.figure(10)
		plt.plot(dt_list , error_list )
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel('$dt$')
		plt.ylabel('$L^{\infty}$ error')
		figure = plt.gcf() # get current figure
		figure.set_size_inches(20, 15)		
		plt.savefig(outFileName, format='eps', dpi=1000)
	
	elif (option == 'epsilon') :
		epsilon_list = data['epsilon']
		error_list = data['error']
		print "epsilon List = " , epsilon_list
		print "error List = " , error_list
		
		params = {'legend.fontsize': 40,
		          'legend.linewidth': 2}
		matplotlib.rcParams.update({'font.size': 50})
		pylab.rcParams.update(params)
		plt.rc('text', usetex=True)
		plt.figure(10)
		plt.plot(epsilon_list , error_list,label='SP', marker='x',linestyle = 'None')
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel('$\epsilon$')
		plt.ylabel('$L^{\infty}$ error')
		plt.legend(loc=1)
		
		# ERROR FOR LTP, N=2500, dt = 0.001 = 0.0582447947973
		#~ error_LTP_list=[0.0582447947973 for i in range(0,len(epsilon_list))]
		#~ plt.plot(epsilon_list , error_LTP_list,label='LTP' )
		
		
		#~ mng = plt.get_current_fig_manager()
		#~ mng.resize(*mng.window.maxsize())
		m,b = np.polyfit(epsilon_list, error_list, 1)
		print m , b
		plt.plot(epsilon_list, m * epsilon_list + b, color='red')
		#~ plt.savefig(outFileName, format='png', dpi=1000, bbox_inches='tight')
		figure = plt.gcf() # get current figure
		figure.set_size_inches(20, 15)
		plt.savefig(outFileName, format=str(extension[1:]), dpi=100)
		plt.show()
	
	else :
		#raise ValueError("parameter option is to plot error against 'dt' or 'N'")
		print "parameter option is to plot error against 'dt' or 'N'"




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
	
	option=sys.argv[2]
	outFileName = sys.argv[3]
	nb_of_files = len(sys.argv[4:])
	print "\n Number of files from which errors will be computed = " , nb_of_files , "\n"
	print "option choosen : " , option
	print "output file name  : " , outFileName
	
	extension = os.path.splitext(outFileName)[1] # either .eps or .png	
	
	label_list = [os.path.splitext(f)[0] for f in sys.argv[4:]]
	label_list = [label.split('/')[1] for label in label_list]
	label_list = ["$"+label.split('_')[0]+"_{"+label.split('_')[1]+"}$" for label in label_list]
	plt.figure(10)
	if (option == 'epsilon') :					
		for i in range(0,int(nb_of_files)) :			
			inFileName=sys.argv[4+i]
			print "input file name = " , inFileName
			
			data=np.genfromtxt(fname=inFileName, names = True)			
			epsilon_list = data['epsilon']
			error_list = data['error']
			
			#~ delta_list = 
			print "all errors : " , error_list
			
			plt.xscale('log')
			plt.yscale('log')
			plt.xlabel('$\epsilon$')
			plt.ylabel('$L^{\infty}$ error')
			plt.plot(epsilon_list , error_list, label = '$'+str(label_list[i])+'$')
			plt.legend(loc=1)
	if (option == 'N') :			
		plt.legend(loc=2)
		for i in range(0,int(nb_of_files)) :
			inFileName=sys.argv[3+i]
			data=np.genfromtxt(fname=inFileName, names = True)			
			N_list = data['N']
			error_list = data['error']
			plt.xscale('log')
			plt.yscale('log')
			plt.xlabel('$N$')
			plt.ylabel('$L^{\infty}$ error')
			plt.plot(N_list , error_list, label = '$'+str(label_list[i])+'$')

	
	figure = plt.gcf() # get current figure
	figure.set_size_inches(20, 15)
	plt.savefig(outFileName, format=str(extension[1:]), dpi=100)
		
	plt.show()
