import sys
import os
import csv

import numpy as np

name_input_file = sys.argv[1] # format = ../outputs/output.txt
val_dt = sys.argv[2]
val_Nx = sys.argv[3]

errors_file = open("errors_dt="+str(val_dt)+"_Nx="+str(val_Nx)+'.out','a')

data = np.genfromtxt(fname=name_input_file, names = True)

with open(name_input_file) as f:
    content = f.readlines()
#remove whitespace characters like `\n` at the end of each line
content = [x.strip() for x in content]
list_name_data_file = [f.split('.csv')[0]+".csv" for f in content[1:]]

print "list_name_data_file = " , list_name_data_file
#~ print result

list_Linfty_error = np.genfromtxt(fname=name_input_file, names = True)['Linfinity_error']
list_L1_error = np.genfromtxt(fname=name_input_file, names = True)['L1_error']
list_L2_error = np.genfromtxt(fname=name_input_file, names = True)['L2_error']

i=0
errors_file.write('epsilon\tLinfinity_error\tL1_error\tL2_error\n')
for name_data_file in list_name_data_file :
    Linfty_error = list_Linfty_error[i]
    L1_error = list_L1_error[i]
    L2_error = list_L2_error[i]
    print "file location : " , "../../"+str(name_data_file)
    with open("../../"+str(name_data_file), 'rb') as f:
        reader = csv.reader(f, delimiter='=')
        dt = 0.
        Nx= 0.
        epsilon = 0.
        for row in reader:
            if (row[0] == 'dt') : dt = float(row[1])
            elif (row[0] == 'Nx') : Nx = float(row[1])
            elif (row[0] == 'epsilon') : epsilon = float(row[1])
        print dt , Nx , epsilon
        if (abs(float(val_dt) - float(dt)) < 0.00000001) and (abs(float(val_Nx) - float(Nx)) < 0.00000001):
            errors_file.write(str(epsilon)+'\t'+str(Linfty_error)+'\t'+str(L1_error)+'\t'+str(L2_error)+'\n')
    i+=1
errors_file.close()
