To compil fortran modules :

1) cd ./trunk/python_srcs/

2) python setup.py build_ext --inplace



To launch code with one run :
           ./code_LTP_2D.py data.csv data.out



To launch code with several runs :
	./launchRuns ./data/inputs/LTP_dt=0.0001_Nx=50-100-200-400 
					./data/outputs

Comments :

- code_LTP_2D.py is an executable but it can be launched as 'python code_LTP_2D.py data.csv output.txt

- data.csv is a csv file which contains all parameters the code needs


- data.out is a txt file which contains the relative errors in L^\infty, L^1 and L^2 norms
  

it will also save matrices at several time steps during the run
- if we launch several runs : 

we must give as arguments : 
	

1) the directory to read the input files 
	

2) the directory to write the output files 
















