#!/bin/bash
#  Job name is myjob, this will be displayed in the list of running jobs
#PBS -N xz
#  stdout (console output) and stderr (error messages) are passed together into one output textfile
#PBS -j oe
#  One node and 56 processor cores are requested on this node (maximum 128 processor cores are possible).
#PBS -l nodes=1:ppn=128
#  Maximum computing time of the job in hours:minutes:seconds.
#PBS -l walltime=1-23:00:00
# It is assumed that the program is called myfile and is located in the directory /home/htc/myname/myfolder.
 #SBATCH --mem=32768 
jobdir="/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/emiModelBddcSdc"
pwd
# mkdir -p "${jobdir}"
 
cd "${jobdir}"
pwd
 
# Application programs that employ multiple processor cores are always started with the slurm command srun.
# The option -B *:*:* causes all processor cores requested with the directive "#PBS -l nodes=1:ppn=xx" to be used by the myfile program (e.g. by threaded computing).
# The -n1 specification causes exactly one instance of the program myfile to be started (which is typical for the present Kaskade7 applications - other specifications for the -n value could make sense e.g. for MPI programs).
# Just use this like you would use your console to start your program.
 
srun -B *:*:* -n1 ./emiModel --refine 1

