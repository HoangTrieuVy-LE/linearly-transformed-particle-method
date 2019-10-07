import os
import platform
import subprocess

libs_gfortran = ['gfortran']
libs_mpifortran = ['mpi_f90', 'mpi_f77']

path_base = os.path.join(
    'build',
    'temp.' + platform.system().lower() + '-'
    + platform.machine() + '-'
    + '.'.join(platform.python_version_tuple()[:2]))

def build_objects_from_fortran(sources):
    objects = []

    for source in sources:
        path_dir, name = source.rsplit(os.path.sep, 1)
        path_dir_objet = os.path.join(path_base, path_dir)
        if not os.path.exists(path_dir_objet):
            os.makedirs(path_dir_objet)
        path_objet = os.path.join(
            path_dir_objet,
            os.path.splitext(name)[0] + '.o')
        objects.append(os.path.relpath(path_objet))
        command_compile_fortran_mod = (
            'mpif90 ' + ' -O3 -fPIC -J ' + path_dir_objet + ' '
            + source + ' -c -o ' + path_objet)
        print(command_compile_fortran_mod)

        code = subprocess.check_output(command_compile_fortran_mod, shell=True)

        print(code)

    return objects

if __name__ == '__main__':
    path_source = 'pack0/source'
    objects = build_objects_from_fortran([path_source + '/modf.f90'])
    print(objects)