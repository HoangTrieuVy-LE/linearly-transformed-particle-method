import numpy
#  ============= Possible initial data ================================
#Lx1 = -0.5 
#Lx2 = 0.5
#Ly1 = -1.  
#Ly2 = 1.  

### for reverse swirling velocity field test. Use rhoini_cone.
#Lx1 = 0. 
#Lx2 = 1.
#Ly1 = 0.  
#Ly2 = 1.  


### For linear velocity field test. Use rhoini_gaus.
Lx1 = -2. 
Lx2 = 2.
Ly1 = -2.  
Ly2 = 2.  


r = 1.
sigma_gaus = 0.2
x0 = 0.5 # x center of initial density
y0 = 0.25 # y center of initial density

outfile_name="TEST.txt"
myfile=open('./out/'+str(outfile_name),'a')
myfile.write('\n\n#______________________________________________________________________#\n')
myfile.write('#\n')
myfile.write('### NEW RUN\n')
myfile.write('#______________________________________________________________________#\n\n')
myfile.write("Lx1 = "+str(Lx1)+'\n')
myfile.write("Lx2 = "+str(Lx2)+'\n')
myfile.write("Ly1 = "+str(Ly1)+'\n')
myfile.write("Ly2 = "+str(Ly2)+'\n')
myfile.write("sigma_gaus = "+str(sigma_gaus)+'\n')
myfile.write("x0 = "+str(x0)+'\n')
myfile.write("y0 = "+str(y0)+'\n')
myfile.close()
def rhoini_disque(x, y, r=1.):
    n = numpy.sqrt(x*x+y*y)
    if n<=r:
        return 1.
    else:
        return 0.  

def rhoini_rectangle(x, y):
    if Lx1<x<Lx2:
        ax = 1.
    else:
        ax = 0.  

    if Ly1<y<Ly2:
        ay =1.
    else:       
        ay = 0.
    return ax * ay

#def rhoini_gaus(x, y):
    #if Lx1<x<Lx2:
        ##ax =1+0.01*cos(k*x)
        ##ax = math.exp(-30 * (x - 0.5) * (x - 0.5)) + 2 * math.exp(-50 * (x + 0.3) * (x + 0.3))
        #ax = math.exp(- 2.*x  * x  /(sigma**2)) 
    #else:
        #ax = 0.  

    #if Ly1<y<Ly2:
        #ay =1.
        ##ay = math.exp(-y*y/2)/sqrt(2*3.14)
    #else:       
        #ay = 0.
    #return ax * ay
   
#Gaussienne centree en 0 par defaut (changer ? )
def rhoini_gaus(x, y):
	if Lx1<x<Lx2:
        #ax =1+0.01*cos(k*x)
        #ax = math.exp(-30 * (x - 0.5) * (x - 0.5)) + 2 * math.exp(-50 * (x + 0.3) * (x + 0.3))
		ax = numpy.exp(-x  * x  /(sigma_gaus**2)) 
	else:
		ax = 0.  
	if Ly1<y<Ly2:
		ay = numpy.exp(- y  * y  /(sigma_gaus**2))
	else:       
		ay = 0.
	return ax * ay
    
def grad_rohini_gaus(x,y) : 
    grad_gaus=numpy.zeros(2)
    if (Lx1<x<Lx2) & (Ly1<y<Ly2):
        #ax =1+0.01*cos(k*x)
        #ax = math.exp(-30 * (x - 0.5) * (x - 0.5)) + 2 * math.exp(-50 * (x + 0.3) * (x + 0.3))
        ax = numpy.exp(-x  * x  /(sigma_gaus**2)) 
        ax_deriv= (-2*x/sigma_gaus**2) * numpy.exp(-x  * x  /(sigma_gaus**2)) 
        ay = numpy.exp(- y  * y  /(sigma_gaus**2))
        ay_deriv = (-2*y/sigma_gaus**2) * numpy.exp(- y  * y  /(sigma_gaus**2))
        grad_gaus[0] = ax_deriv * ay
        grad_gaus[1] = ax * ay_deriv
    return grad_gaus

def rhoini_cone(x,y) :
	return max(0,1-(20./3.)*numpy.sqrt((x-x0)**2+(y-y0)**2))

