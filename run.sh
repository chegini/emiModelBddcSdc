#!/bin/bash
#  Job name is myjob, this will be displayed in the list of running jobs
#PBS -N xz
#  stdout (console output) and stderr (error messages) are passed together into one output textfile
#PBS -j oe
#  One node and 56 processor cores are requested on this node (maximum 128 processor cores are possible).
#PBS -l nodes=1:ppn=56
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

srun -B *:*:* -n1 ./emiModel 
# srun -B *:*:* -n1 ./emiModel --write_to_file true --input "./input/robin_mesh.vtu" --extra_set "./input/robin_extracellular.txt" --intra_set "./input/robin_intracellular.txt" --excited "./input/robin_early_excited.txt" --dir "/scratch/htc/fchegini/robin1_output" --matlab_dir "/scratch/htc/fchegini/robin1_matlab_dir" --refine 1

# srun -B *:*:* -n1 ./emiModel --write_to_file true --input "./input/kermit_mesh.vtu" --extra_set "./input/kermit_extracellular.txt" --intra_set "./input/kermit_intracellular.txt" --excited "./input/kermit_early_excited.txt" --dir "/scratch/htc/fchegini/kermit_output" --matlab_dir "/scratch/htc/fchegini/kermit_matlab_dir"

# srun -B *:*:* -n1 ./emiModel --write_to_file true --input "./input/gonzo_mesh.vtu" --extra_set "./input/gonzo_extracellular.txt" --intra_set "./input/gonzo_intracellular.txt" --excited "./input/gonzo_early_excited.txt" --dir "/scratch/htc/fchegini/gonzo_output" --matlab_dir "/scratch/htc/fchegini/gonzo_matlab_dir" --nThreads 128

# srun -B *:*:* -n1 ./emiModel --write_to_file true --input "./input/animal_mesh.vtu" --extra_set "./input/animal_extracellular.txt" --intra_set "./input/animal_intracellular.txt" --excited "./input/animal_early_excited.txt" --dir "/scratch/htc/fchegini/animal_output" --matlab_dir "/scratch/htc/fchegini/animal_matlab_dir" 


