import numpy as np

U_ltp = np.genfromtxt('U_ltp_t=0.051.log')
U_sp = np.genfromtxt('U_sp_t=0.051.log')
N=2500
list_errors_x = []
list_errors_y = []
for k in range(0,N) : 
    list_errors_x.append(abs(U_sp[k,0] - U_ltp[k,0]))
    list_errors_y.append(abs(U_sp[k,1] - U_ltp[k,1]))
    print( k , U_sp[k,0] , U_sp[k,1] , U_ltp[k,0], U_ltp[k,1] , abs(U_sp[k,0] - U_ltp[k,0]) , abs(U_sp[k,1] - U_ltp[k,1]))
print (max(list_errors_x))
print (max(list_errors_y))
