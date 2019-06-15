#!/bin/bash
#SBATCH -J PartsBarenblatt
#SBATCH -N 2
#SBATCH -n 2
#SBATCH --cpus-per-task=20
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-core=1
#SBATCH --time=48:00:00

#(--ntasks-per-node) * (-N ) = -n
#--ntasks-per-core=1 # OBLIGATOIRE

### options SBATCH  :
#-n, --ntasks=ntasks         number of tasks to run
#-N, --nodes=N               number of nodes on which to run (N = min[-max])
#--ntasks-per-node=n         number of tasks to invoke on each node
#--ntasks-per-core=n         number of tasks to invoke on each core
#-c, --cpus-per-task=ncpus   number of cpus required per task


module purge
module load intel chdb/intelmpi
#module load fftw/3.3.4_openmp
module load python/2.7.6


#export OMP_NUM_THREADS=#SBATCH --cpus-per-task
export OMP_NUM_THREADS=20
MKL_NUM_THREADS=$OMP_NUM_THREADS 
#export OMP_PROC_BIND=true 

INPUT='./data/inputs/'
OUTPUT='./data/outputs_chdb/' #OUTPUT must not exist when it's launched

#rmdir $OUTPUT 2> /tmp/null

### launch runs
mpirun chdb --on-error error.txt --report report.txt --in-dir $INPUT --in-type csv --out-dir $OUTPUT --out-files %out-dir%/%path% --command "./code_LTP_2D.py %in-dir%/%name% %out-dir%/%name%.out > log/output_${SLURM_JOBID}_%name%.log"

### get infos on job (time, memory, etc...)
#jobinfo ${SLURM_JOBID}
#infoincidentjob
