import sys
import os
import csv
import glob

import numpy as np

input_dir_name = sys.argv[1] # format = ../outputs/output.txt
output_error_name = sys.argv[2]

list_errors_file = glob.glob(str(input_dir_name)+'*.out')
list_data_file = []
list_error_L_infty = []
list_error_L_infty_X = []
list_error_L_infty_D = []
list_error_L1 = []
list_error_L2 = []
list_dt = []
list_Nx = []
list_Ny = []
list_epsilon = []

nb_runs = len(list_errors_file)

for fname in list_errors_file :    
    with open(fname) as f :
        lines = f.readlines()   
        data_fname = lines[1].replace("\t",",").split(',')[0]
        print data_fname
        error_L_infty = lines[1].replace("\t",",").split(',')[1]
        error_L1 = lines[1].replace("\t",",").split(',')[2]
        error_L2 = lines[1].replace("\t",",").split(',')[3]
        error_L_infty_X = lines[1].replace("\t",",").split(',')[4]
        error_L_infty_D = lines[1].replace("\t",",").split(',')[5]
        list_error_L_infty.append(error_L_infty)
        list_error_L_infty_X.append(error_L_infty_X)
        list_error_L_infty_D.append(error_L_infty_D)
        list_error_L1.append(error_L1)
        list_error_L2.append(error_L2)
        list_data_file.append(data_fname)
    
    with open("../../"+str(data_fname), 'rb') as fdata :        
        reader = csv.reader(fdata, delimiter='=')
        for row in reader:
            if (row[0] == 'epsilon') : 	epsilon = float(row[1])
            if (row[0] == 'dt') : 		dt = float(row[1])
            if (row[0] == 'Nx') : 		Nx = float(row[1])
            if (row[0] == 'Ny') : 		Ny = float(row[1])
        
    list_dt.append(dt)
    list_epsilon.append(epsilon)
    list_Nx.append(Nx)
    list_Ny.append(Ny)
    

print list_error_L_infty_X
ferror = open(str(output_error_name), 'w')
ferror.write('dt\tepsilon\tNx\tNy\tL_infinity_error\tL1_error\tL2_error\tLinfinity_error_X\tLinfinity_error_D\n')
for i in range(0,nb_runs) :
    ferror.write(str(list_dt[i])+'\t'+str(list_epsilon[i])+'\t'+str(list_Nx[i])+'\t'+str(list_Ny[i])+'\t'+str(list_error_L_infty[i])+'\t'+str(list_error_L1[i])+'\t'+str(list_error_L2[i])+'\t'+str(list_error_L_infty_X[i])+'\t'+str(list_error_L_infty_D[i]))
    
ferror.close()
