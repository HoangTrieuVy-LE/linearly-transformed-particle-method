#!/usr/bin/env python
import os
import sys
#import commands
import time
import glob
import numpy as np


# This script takes 2 arguments :
#   - a directory to find all data files that contain parameters for the code
#   - a directory to write all output files

# REQUIREMENTS : -all data files must have a '.csv' extension and must be in 'inputs_dir' directory
#                -all output files will be written with the same name as the data files but with '.out' as extension
#                 e.g. : if file name in inputs_dir is data_1.csv the corresponding output file will be written in outputs_dir/data_1.out

inputs_dir = sys.argv[1]
outputs_dir = sys.argv[2]



#if the output directory doesn't exist, we create it
if (not(os.path.isdir(str(outputs_dir)))) :
    os.makedirs(str(outputs_dir))


input_files_list = glob.glob(str(inputs_dir)+'/*.csv')


print (input_files_list)

for f_in in input_files_list :
    file_name_no_ext = f_in.split('/')[-1].split('.')[0]
    f_out = str(outputs_dir)+'/'+file_name_no_ext+'.out'
    print(f_out)
#    f_out = str(outputs_dir)+'/'+f_in+'.out'
#    cmd = './code_LTP_2D_diffusion.py '+str(f_in)+' '+str(f_out)    
    cmd = 'python code_LTP_2D_diffusion.py '+str(f_in)+' '+str(f_out)    

    print ("launch : " , cmd )
#    os.system("date")
    
#    launch_date = time.strftime("%Y-%m-%d_%H:%M:%S", time.gmtime())    
#    log_file=open("logs_launchRuns.log", 'a')
#    log_file.write('launch = '+str(cmd)+'\n')
#    log_file.write('date = '+str(launch_date)+'\n')
#    log_file.close()
    os.system(str(cmd))

#inputs_dir = sys.argv[1]
#outputs_dir = sys.argv[2]

#if the output directory doesn't exist, we create it
#if (not(os.path.isdir(str(outputs_dir)))) :
#    os.makedirs(str(outputs_dir))
#
#input_files_list = glob.glob(str(inputs_dir)+'/*.csv')
#
#
#print input_files_list
#
#for f_in in input_files_list :
#	
#    file_name_no_ext = f_in.split('/')[-1].split('.')[0]
#    f_out = str(outputs_dir)+'/'+file_name_no_ext+'.out'
#    cmd = './code_LTP_2D_diffusion.py '+str(f_in)+' '+str(f_out)    
#    print "launch : " , cmd 
#    os.system("date")
#    
#    launch_date = time.strftime("%Y-%m-%d_%H:%M:%S", time.gmtime())    
#    log_file=open("logs_launchRuns.log", 'a')
#    log_file.write('launch = '+str(cmd)+'\n')
#    log_file.write('date = '+str(launch_date)+'\n')
#    log_file.close()
#    os.system(str(cmd))
#
# commands.getstatusoutput returns a status and the output result as a tuple. So we take only second arguemnt as an intger(the result)
# nb_runs = int(commands.getstatusoutput('ls '+str(inputs_dir)+'/*.csv | wc -l')[1])

#~ for i in range(1,int(nb_runs)) :
    #~ 
    #~ cmd = './code_LTP_2D.py ./data/inputs/data_'+str(i)+'.csv ./data/outputs/data_'+str(i)+'.out'
    #~ print cmd
    #~ os.system(str(cmd))
    
