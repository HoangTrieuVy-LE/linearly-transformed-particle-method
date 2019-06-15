import numpy as np
import sys

name_file = sys.argv[1]

indices_part = np.genfromtxt(fname=name_file, usecols = 0)
det_Dk = np.genfromtxt(fname=name_file, usecols = 1)
X_0 = np.genfromtxt(fname=name_file, usecols = 2)
X_1 = np.genfromtxt(fname=name_file, usecols = 3)
delta_X = np.genfromtxt(fname=name_file, usecols = 4)
time = np.genfromtxt(fname=name_file, usecols = 5)

sorted_det_Dk = np.sort(det_Dk)
args_sorted_det_Dk = np.argsort(det_Dk)
indices_part = indices_part[args_sorted_det_Dk]
time = time[args_sorted_det_Dk]
X_0 = X_0[args_sorted_det_Dk]
X_1 = X_1[args_sorted_det_Dk]
delta_X = delta_X[args_sorted_det_Dk]



N = len(det_Dk)

#~ arg_sort_labels = np.argsort(labels_list_1)
#~ labels = labels_list_1[arg_sort_labels]

for i in range(0,int(N)) : 	
	print i, indices_part[i] , sorted_det_Dk[i] , X_0[i], X_1[i], delta_X[i], time[i]


