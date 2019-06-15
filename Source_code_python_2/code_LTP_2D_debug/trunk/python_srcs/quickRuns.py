import numpy as np
import os

for i in range(2,10) : 
	print i
	command = "python code_LTP_2D.py data/inputs/data_"+str(i)+".csv data/outputs/output.txt"
	print "launch command : " , command
	os.system(str(command))


#~ os.system()
