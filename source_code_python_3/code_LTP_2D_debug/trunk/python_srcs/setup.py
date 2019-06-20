import sys
#sys.path.append('../../trunk/fortran_srcs/')

extra_link_args=[]

if sys.platform=='darwin':	extra_link_args=['-framework','veclib', '-fopenmp' , '-lgomp', '/openmp']

from numpy.distutils.core import setup, Extension

extra_link_args.append('-fopenmp') # option for fortran OpenMP

opt_flags = ['-fopenmp' , '-lgomp'  ,'-g','-fbounds-check' ]

#~ opt_flags = ['g']	

files1 = ['../../trunk/fortran_srcs/calculsfortran_rec.f90'\
         ,'../../trunk/fortran_srcs/calculsfortran_ini.f90'\
         ,'../../trunk/fortran_srcs/calculs_2D.f90'\
         ,'../../trunk/fortran_srcs/2D/mod_particle_2D.f90'\
         ,'../../trunk/fortran_srcs/2D/MPI_2D_spatialisation.f90'\
         ,'../../trunk/fortran_srcs/2D/MPI_2D_structures.f90'\
         ,'../../trunk/fortran_srcs/2D/tri_casier_method.f90'\

                ] # fortran modules used by main python code
                        
ext1 = Extension(name='_calculsfor_f90',
                 sources = files1,
                 extra_compile_args=opt_flags,
                 extra_f90_compile_args=opt_flags,
                 extra_link_args=extra_link_args)
                 



setup(name        = "_modules_f90",
      version     = '0.1',
      description = "modules for a python code",
      author      = "Fred a partir du fichier de Martin",
      author_email= '',
      url         = '',
      ext_modules = [ext1])
