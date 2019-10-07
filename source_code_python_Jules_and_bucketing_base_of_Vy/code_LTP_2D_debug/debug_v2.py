import numpy as np

data_ltp = np.genfromtxt('ltp.debug', usenames=True)
data_sp = np.genfromtxt('sp.debug', usenames=True)

ltp_tmp_val_x= data_ltp['tmp_val_x']
ltp_tmp_val_y= data_ltp['tmp_val_y']
sp_tmp_val_x= data_sp['tmp_val_x']
sp_tmp_val_y= data_sp['tmp_val_y']

list_errors_x = []
list_errors_y = []

N=1024
for j in range(0,2*N) : 
	print (j , abs(ltp_tmp_val_x[j] - sp_tmp_val_x[j]) , abs(ltp_tmp_val_y[j] - sp_tmp_val_y[j]))
