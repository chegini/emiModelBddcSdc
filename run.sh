#!/bin/bash
#  Job name is myjob, this will be displayed in the list of running jobs
#PBS -N r4
#  stdout (console output) and stderr (error messages) are passed together into one output textfile
#PBS -j sol
#  One node and 56 processor cores are requested on this node (maximum 128 processor cores are possible).
#PBS -l nodes=1:ppn=rTB
#  Maximum computing time of the job in hours:minutes:seconds.
##PBS -l walltime=29-23:00:00
#PBS -l walltime=10-23:00:00
# It is assumed that the program is called myfile and is located in the directory /home/htc/myname/myfolder.
#####SBATCH --mem=32768 ls -l
#SBATCH --mem=102400
##SBATCH --mem=1000000
#SBATCH --partition=big 
##SBATCH --partition=high-mem
##SBATCH --nodelist=htc-cmp023
##SBATCH --partition=cno

jobdir="/data/numerik/people/fchegini/project/fatemeh/MicroCard/current_papers/BDDC_petsc_SDC_2023/gatevariables/20Nov2023/BDDC/emiModelBddcSdc"
pwd
# mkdir -p "${jobdir}"
 
cd "${jobdir}"
pwd

# mass + stiffness + next + time step
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny1nz.vtu" --extra_set "./input/extracellular3nx3ny1nz.txt" --intra_set "./input/intracellular3nx3ny1nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny1nz_mass_stiffness_diff_time/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny1nz_mass_stiffness_diff_time/matlab_dir" --interfacetypes 7 --T 0.02
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny2nz.vtu" --extra_set "./input/extracellular3nx3ny2nz.txt" --intra_set "./input/intracellular3nx3ny2nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny2nz_mass_stiffness_diff_time/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny2nz_mass_stiffness_diff_time/matlab_dir" --interfacetypes 7 --T 0.02
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny3nz.vtu" --extra_set "./input/extracellular3nx3ny3nz.txt" --intra_set "./input/intracellular3nx3ny3nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny3nz_mass_stiffness_diff_time/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny3nz_mass_stiffness_diff_time/matlab_dir" --interfacetypes 7 --T 0.02
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid4nx4ny4nz.vtu" --extra_set "./input/extracellular4nx4ny4nz.txt" --intra_set "./input/intracellular4nx4ny4nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube4nx4ny4nz_mass_stiffness_diff_time/output" --matlab_dir "/scratch/htc/fchegini/cube4nx4ny4nz_mass_stiffness_diff_time/matlab_dir" --interfacetypes 7 --T 0.02
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid5nx5ny5nz.vtu" --extra_set "./input/extracellular5nx5ny5nz.txt" --intra_set "./input/intracellular5nx5ny5nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube5nx5ny5nz_mass_stiffness_diff_time/output" --matlab_dir "/scratch/htc/fchegini/cube5nx5ny5nz_mass_stiffness_diff_time/matlab_dir" --interfacetypes 7 --T 0.02
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid6nx6ny6nz.vtu" --extra_set "./input/extracellular6nx6ny6nz.txt" --intra_set "./input/intracellular6nx6ny6nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube6nx6ny6nz_mass_stiffness_diff_time/output" --matlab_dir "/scratch/htc/fchegini/cube6nx6ny6nz_mass_stiffness_diff_time/matlab_dir" --interfacetypes 7 --T 0.02
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid7nx7ny7nz.vtu" --extra_set "./input/extracellular7nx7ny7nz.txt" --intra_set "./input/intracellular7nx7ny7nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube7nx7ny7nz_mass_stiffness_diff_time/output" --matlab_dir "/scratch/htc/fchegini/cube7nx7ny7nz_mass_stiffness_diff_time/matlab_dir" --interfacetypes 7 --T 0.02
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid8nx8ny8nz.vtu" --extra_set "./input/extracellular8nx8ny8nz.txt" --intra_set "./input/intracellular8nx8ny8nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube8nx8ny8nz_mass_stiffness_diff_time/output" --matlab_dir "/scratch/htc/fchegini/cube8nx8ny8nz_mass_stiffness_diff_time/matlab_dir" --interfacetypes 7 --T 0.02




# new mesh weak scability solution for T
#srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny3nz_new.vtu" --extra_set "./input/extracellular3nx3ny3nz.txt" --intra_set "./input/intracellular3nx3ny3nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny1nz_solution/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny1nz_solution/matlab_dir" --interfacetypes 7 --T 1.0
#srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny3nz.vtu" --extra_set "./input/extracellular3nx3ny3nz.txt" --intra_set "./input/intracellular3nx3ny3nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny1nz_solution/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny1nz_solution/matlab_dir" --interfacetypes 7 --T 10.0 --run_implicit_CG_BDDC false
#srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny1nz.vtu" --extra_set "./input/extracellular3nx3ny1nz.txt" --intra_set "./input/intracellular3nx3ny1nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny1nz/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny1nz/matlab_dir" --interfacetypes 7 --T 4.5 --run_implicit_CG_BDDC false

# new mesh weak scability
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny1nz.vtu" --extra_set "./input/extracellular3nx3ny1nz.txt" --intra_set "./input/intracellular3nx3ny1nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny1nz/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny1nz/matlab_dir" --interfacetypes 7
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny2nz.vtu" --extra_set "./input/extracellular3nx3ny2nz.txt" --intra_set "./input/intracellular3nx3ny2nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny2nz/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny2nz/matlab_dir" --interfacetypes 7
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny3nz.vtu" --extra_set "./input/extracellular3nx3ny3nz.txt" --intra_set "./input/intracellular3nx3ny3nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny3nz/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny3nz/matlab_dir" --interfacetypes 7
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid4nx4ny4nz.vtu" --extra_set "./input/extracellular4nx4ny4nz.txt" --intra_set "./input/intracellular4nx4ny4nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube4nx4ny4nz/output" --matlab_dir "/scratch/htc/fchegini/cube4nx4ny4nz/matlab_dir" --interfacetypes 7
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid5nx5ny5nz.vtu" --extra_set "./input/extracellular5nx5ny5nz.txt" --intra_set "./input/intracellular5nx5ny5nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube5nx5ny5nz/output" --matlab_dir "/scratch/htc/fchegini/cube5nx5ny5nz/matlab_dir" --interfacetypes 7
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid6nx6ny6nz.vtu" --extra_set "./input/extracellular6nx6ny6nz.txt" --intra_set "./input/intracellular6nx6ny6nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube6nx6ny6nz/output" --matlab_dir "/scratch/htc/fchegini/cube6nx6ny6nz/matlab_dir" --interfacetypes 7
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid7nx7ny7nz.vtu" --extra_set "./input/extracellular7nx7ny7nz.txt" --intra_set "./input/intracellular7nx7ny7nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube7nx7ny7nz/output" --matlab_dir "/scratch/htc/fchegini/cube7nx7ny7nz/matlab_dir" --interfacetypes 7
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid8nx8ny8nz.vtu" --extra_set "./input/extracellular8nx8ny8nz.txt" --intra_set "./input/intracellular8nx8ny8nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube8nx8ny8nz/output" --matlab_dir "/scratch/htc/fchegini/cube8nx8ny8nz/matlab_dir" --interfacetypes 7
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid9nx9ny9nz.vtu" --extra_set "./input/extracellular9nx9ny9nz.txt" --intra_set "./input/intracellular9nx9ny9nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube9nx9ny9nz/output" --matlab_dir "/scratch/htc/fchegini/cube9nx9ny9nz/matlab_dir" --interfacetypes 7
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid10nx10ny10nz.vtu" --extra_set "./input/extracellular10nx10ny10nz.txt" --intra_set "./input/intracellular10nx10ny10nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube10nx10ny10nz/output" --matlab_dir "/scratch/htc/fchegini/cube10nx10ny10nz/matlab_dir" --interfacetypes 7

# new mesh strong scability
#srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny3nz.vtu" --extra_set "./input/extracellular3nx3ny3nz.txt" --intra_set "./input/intracellular3nx3ny3nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny3nzr1/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny3nzr1/matlab_dir" --interfacetypes 7 --refine 1
#srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny3nz.vtu" --extra_set "./input/extracellular3nx3ny3nz.txt" --intra_set "./input/intracellular3nx3ny3nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny3nzr2/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny3nzr2/matlab_dir" --interfacetypes 7 --refine 2
#srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny3nz.vtu" --extra_set "./input/extracellular3nx3ny3nz.txt" --intra_set "./input/intracellular3nx3ny3nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny3nzr3/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny3nzr3/matlab_dir" --interfacetypes 7 --refine 3
srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny3nz.vtu" --extra_set "./input/extracellular3nx3ny3nz.txt" --intra_set "./input/intracellular3nx3ny3nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny3nzr4/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny3nzr4/matlab_dir" --interfacetypes 7 --refine 4



#srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny3nz.vtu" --extra_set "./input/extracellular3nx3ny3nz.txt" --intra_set "./input/intracellular3nx3ny3nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny3nzr2o2/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny3nzr2o2/matlab_dir" --interfacetypes 7 --refine 2 --order 2

# more refinement
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny3nz.vtu" --extra_set "./input/extracellular3nx3ny3nz.txt" --intra_set "./input/intracellular3nx3ny3nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny3nzO2r1/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny3nzO2r1/matlab_dir" --interfacetypes 7 --order 2 --refine 1
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny3nz.vtu" --extra_set "./input/extracellular3nx3ny3nz.txt" --intra_set "./input/intracellular3nx3ny3nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny3nzO2r2/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny3nzO2r2/matlab_dir" --interfacetypes 7 --order 2 --refine 2
# srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny3nz.vtu" --extra_set "./input/extracellular3nx3ny3nz.txt" --intra_set "./input/intracellular3nx3ny3nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny3nzO2r3/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny3nzO2r3/matlab_dir" --interfacetypes 7 --order 2 --refine 3

#srun -B *:*:* -n1 ./emiModel --input "./input/emiGrid3nx3ny3nz.vtu" --extra_set "./input/extracellular3nx3ny3nz.txt" --intra_set "./input/intracellular3nx3ny3nz.txt" --excited "./input/emiGrid_earlyExited.txt" --dir "/scratch/htc/fchegini/cube3nx3ny3nzO2r1_T/output" --matlab_dir "/scratch/htc/fchegini/cube3nx3ny3nzO2r1_T/matlab_dir" --interfacetypes 7 --order 2 --refine 1 --T 1.0

# Application programs that employ multiple processor cores are always started with the slurm command srun.
# The option -B *:*:* causes all processor cores requested with the directive "#PBS -l nodes=1:ppn=xx" to be used by the myfile program (e.g. by threaded computing).
# The -n1 specification causes exactly one instance of the program myfile to be started (which is typical for the present Kaskade7 applications - other specifications for the -n value could make sense e.g. for MPI programs).
# Just use this like you would use your console to start your program.

#srun -B *:*:* -n1 ./emiModel --refine 0

# tetsing parallization
# ./emiModel --input "./input/10Cells3d_10extra_mesh.vtu" --extra_set "./input/10Cells3d_10extra_list_extracellular.txt" --intra_set "./input/10Cells3d_10extra_list_intracellular.txt" --excited "./input/10Cells3d_10extra_early_excited.txt"

# srun -B *:*:* -n1 ./emiModel
#----------------------------------------------------------------------------------------------------------------------------------------
# BDDC scabality 
#----------------------------------------------------------------------------------------------------------------------------------------
# scabality test
# srun -B *:*:* -n1 ./emiModel --input "./input/10Cells3d_10extra_mesh.vtu" --extra_set "./input/10Cells3d_10extra_list_extracellular.txt" --intra_set "./input/10Cells3d_10extra_list_intracellular.txt" --excited "./input/10Cells3d_10extra_early_excited.txt" --dir "/scratch/htc/fchegini/cells20/output" --matlab_dir "/scratch/htc/fchegini/cells20/matlab_dir"  --refine 2
# srun -B *:*:* -n1 ./emiModel --input "./input/20Cells3d_20extra_mesh.vtu" --extra_set "./input/20Cells3d_20extra_list_extracellular.txt" --intra_set "./input/20Cells3d_20extra_list_intracellular.txt" --excited "./input/20Cells3d_20extra_early_excited.txt" --dir "/scratch/htc/fchegini/cells40/output" --matlab_dir "/scratch/htc/fchegini/cells40/matlab_dir"  --refine 2
# srun -B *:*:* -n1 ./emiModel --input "./input/40Cells3d_40extra_mesh.vtu" --extra_set "./input/40Cells3d_40extra_list_extracellular.txt" --intra_set "./input/40Cells3d_40extra_list_intracellular.txt" --excited "./input/40Cells3d_40extra_early_excited.txt" --dir "/scratch/htc/fchegini/cells80/output" --matlab_dir "/scratch/htc/fchegini/cells80/matlab_dir"  --refine 2

# srun -B *:*:* -n1 ./emiModel --input "./input/pepe_combi_domi.vtu" --extra_set "./input/pepe_combi_domi_extracellular.txt" --intra_set "./input/pepe_combi_domi_intracellular.txt" --excited "./input/pepe_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells20_pepe/output" --matlab_dir "/scratch/htc/fchegini/cells20_pepe/matlab_dir" --interfacetypes 5
# srun -B *:*:* -n1 ./emiModel --input "./input/rizzo_combi_domi.vtu" --extra_set "./input/rizzo_combi_domi_extracellular.txt" --intra_set "./input/rizzo_combi_domi_intracellular.txt" --excited "./input/rizzo_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells40_rizzo/output" --matlab_dir "/scratch/htc/fchegini/cells40_rizzo/matlab_dir" --interfacetypes 5
# srun -B *:*:* -n1 ./emiModel --input "./input/robin_mesh.vtu" --extra_set "./input/robin_extracellular.txt" --intra_set "./input/robin_intracellular.txt" --excited "./input/robin_early_excited.txt" --dir "/scratch/htc/fchegini/cells88/output" --matlab_dir "/scratch/htc/fchegini/cells88/matlab_dir"  --interfacetypes 5
# srun -B *:*:* -n1 ./emiModel --input "./input/kermit_mesh.vtu" --extra_set "./input/kermit_extracellular.txt" --intra_set "./input/kermit_intracellular.txt" --excited "./input/kermit_early_excited.txt" --dir "/scratch/htc/fchegini/cells172/output" --matlab_dir "/scratch/htc/fchegini/cells172/matlab_dir"  --interfacetypes 5
# srun -B *:*:* -n1 ./emiModel --input "./input/gonzo_mesh.vtu" --extra_set "./input/gonzo_extracellular.txt" --intra_set "./input/gonzo_intracellular.txt" --excited "./input/gonzo_early_excited.txt" --dir "/scratch/htc/fchegini/cells416/output" --matlab_dir "/scratch/htc/fchegini/cells416/matlab_dir"  
#srun -B *:*:* -n1 ./emiModel --input "./input/animal_mesh.vtu" --extra_set "./input/animal_extracellular.txt" --intra_set "./input/animal_intracellular.txt" --excited "./input/animal_early_excited.txt" --dir "/scratch/htc/fchegini/cells830/output" --matlab_dir "/scratch/htc/fchegini/cells830/matlab_dir"  

# optimality test
# srun -B *:*:* -n1 ./emiModel --input "./input/robin_mesh.vtu" --extra_set "./input/robin_extracellular.txt" --intra_set "./input/robin_intracellular.txt" --excited "./input/robin_early_excited.txt" --dir "/scratch/htc/fchegini/cells88_refine1/output" --matlab_dir "/scratch/htc/fchegini/cells88_refine1/matlab_dir"  --refine 1
# srun -B *:*:* -n1 ./emiModel --input "./input/robin_mesh.vtu" --extra_set "./input/robin_extracellular.txt" --intra_set "./input/robin_intracellular.txt" --excited "./input/robin_early_excited.txt" --dir "/scratch/htc/fchegini/cells88_refine2/output" --matlab_dir "/scratch/htc/fchegini/cells88_refine2/matlab_dir"  --refine 2
# srun -B *:*:* -n1 ./emiModel --input "./input/robin_mesh.vtu" --extra_set "./input/robin_extracellular.txt" --intra_set "./input/robin_intracellular.txt" --excited "./input/robin_early_excited.txt" --dir "/scratch/htc/fchegini/cells88_order2/output" --matlab_dir "/scratch/htc/fchegini/cells88_order2/matlab_dir"  --order 2

# srun -B *:*:* -n1 ./emiModel --input "./input/robin_mesh.vtu" --extra_set "./input/robin_extracellular.txt" --intra_set "./input/robin_intracellular.txt" --excited "./input/robin_early_excited.txt" --dir "/scratch/htc/fchegini/cells88_refine1_order2/output" --matlab_dir "/scratch/htc/fchegini/cells88_refine1_order2/matlab_dir"  --refine 1 --order 2


# srun -B *:*:* -n1 ./emiModel --input "./input/pepe_combi_domi.vtu" --extra_set "./input/pepe_combi_domi_extracellular.txt" --intra_set "./input/pepe_combi_domi_intracellular.txt" --excited "./input/pepe_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells20_order2_new/output" --matlab_dir "/scratch/htc/fchegini/cells20_order2_new/matlab_dir"  --order 2
# srun -B *:*:* -n1 ./emiModel --input "./input/pepe_combi_domi.vtu" --extra_set "./input/pepe_combi_domi_extracellular.txt" --intra_set "./input/pepe_combi_domi_intracellular.txt" --excited "./input/pepe_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells20_refine1_new/output" --matlab_dir "/scratch/htc/fchegini/cells20_refine1_new/matlab_dir"  --refine 1
# srun -B *:*:* -n1 ./emiModel --input "./input/pepe_combi_domi.vtu" --extra_set "./input/pepe_combi_domi_extracellular.txt" --intra_set "./input/pepe_combi_domi_intracellular.txt" --excited "./input/pepe_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells20_refine2_new/output" --matlab_dir "/scratch/htc/fchegini/cells20_refine2_new/matlab_dir"  --refine 2
#srun -B *:*:* -n1 ./emiModel --input "./input/pepe_combi_domi.vtu" --extra_set "./input/pepe_combi_domi_extracellular.txt" --intra_set "./input/pepe_combi_domi_intracellular.txt" --excited "./input/pepe_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells20_refine1_order2_new/output" --matlab_dir "/scratch/htc/fchegini/cells20_refine1_order2_new/matlab_dir"  --refine 1 --order 2




# corners
# srun -B *:*:* -n1 ./emiModel --input "./input/pepe_combi_domi.vtu" --extra_set "./input/pepe_combi_domi_extracellular.txt" --intra_set "./input/pepe_combi_domi_intracellular.txt" --excited "./input/pepe_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells20_pepe_corner/output" --matlab_dir "/scratch/htc/fchegini/cells20_pepe_corner/matlab_dir" --interfacetypes 1
# srun -B *:*:* -n1 ./emiModel --input "./input/rizzo_combi_domi.vtu" --extra_set "./input/rizzo_combi_domi_extracellular.txt" --intra_set "./input/rizzo_combi_domi_intracellular.txt" --excited "./input/rizzo_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells40_rizzo_corner/output" --matlab_dir "/scratch/htc/fchegini/cells40_rizzo_corner/matlab_dir" --interfacetypes 1
# srun -B *:*:* -n1 ./emiModel --input "./input/robin_mesh.vtu" --extra_set "./input/robin_extracellular.txt" --intra_set "./input/robin_intracellular.txt" --excited "./input/robin_early_excited.txt" --dir "/scratch/htc/fchegini/cells88_corner/output" --matlab_dir "/scratch/htc/fchegini/cells88_corner/matlab_dir"  --interfacetypes 1
# srun -B *:*:* -n1 ./emiModel --input "./input/kermit_mesh.vtu" --extra_set "./input/kermit_extracellular.txt" --intra_set "./input/kermit_intracellular.txt" --excited "./input/kermit_early_excited.txt" --dir "/scratch/htc/fchegini/cells172_corner/output" --matlab_dir "/scratch/htc/fchegini/cells172_corner/matlab_dir"  --interfacetypes 1
# srun -B *:*:* -n1 ./emiModel --input "./input/gonzo_mesh.vtu" --extra_set "./input/gonzo_extracellular.txt" --intra_set "./input/gonzo_intracellular.txt" --excited "./input/gonzo_early_excited.txt" --dir "/scratch/htc/fchegini/cells416_corner/output" --matlab_dir "/scratch/htc/fchegini/cells416_corner/matlab_dir"  --interfacetypes 1
# srun -B *:*:* -n1 ./emiModel --input "./input/animal_mesh.vtu" --extra_set "./input/animal_extracellular.txt" --intra_set "./input/animal_intracellular.txt" --excited "./input/animal_early_excited.txt" --dir "/scratch/htc/fchegini/cells830_corner/output" --matlab_dir "/scratch/htc/fchegini/cells830_corner/matlab_dir"  --interfacetypes 1

# corners + edge
# srun -B *:*:* -n1 ./emiModel --input "./input/pepe_combi_domi.vtu" --extra_set "./input/pepe_combi_domi_extracellular.txt" --intra_set "./input/pepe_combi_domi_intracellular.txt" --excited "./input/pepe_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells20_pepe_corner_edge/output" --matlab_dir "/scratch/htc/fchegini/cells20_pepe_corner_edge/matlab_dir" --interfacetypes 3
# srun -B *:*:* -n1 ./emiModel --input "./input/rizzo_combi_domi.vtu" --extra_set "./input/rizzo_combi_domi_extracellular.txt" --intra_set "./input/rizzo_combi_domi_intracellular.txt" --excited "./input/rizzo_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells40_rizzo_corner_edge/output" --matlab_dir "/scratch/htc/fchegini/cells40_rizzo_corner_edge/matlab_dir" --interfacetypes 3
# srun -B *:*:* -n1 ./emiModel --input "./input/robin_mesh.vtu" --extra_set "./input/robin_extracellular.txt" --intra_set "./input/robin_intracellular.txt" --excited "./input/robin_early_excited.txt" --dir "/scratch/htc/fchegini/cells88_corner_edge/output" --matlab_dir "/scratch/htc/fchegini/cells88_corner_edge/matlab_dir"  --interfacetypes 3
# srun -B *:*:* -n1 ./emiModel --input "./input/kermit_mesh.vtu" --extra_set "./input/kermit_extracellular.txt" --intra_set "./input/kermit_intracellular.txt" --excited "./input/kermit_early_excited.txt" --dir "/scratch/htc/fchegini/cells172_corner_edge/output" --matlab_dir "/scratch/htc/fchegini/cells172_corner_edge/matlab_dir"  --interfacetypes 3
# srun -B *:*:* -n1 ./emiModel --input "./input/gonzo_mesh.vtu" --extra_set "./input/gonzo_extracellular.txt" --intra_set "./input/gonzo_intracellular.txt" --excited "./input/gonzo_early_excited.txt" --dir "/scratch/htc/fchegini/cells416_corner_edge/output" --matlab_dir "/scratch/htc/fchegini/cells416_corner_edge/matlab_dir"  --interfacetypes 3
# srun -B *:*:* -n1 ./emiModel --input "./input/animal_mesh.vtu" --extra_set "./input/animal_extracellular.txt" --intra_set "./input/animal_intracellular.txt" --excited "./input/animal_early_excited.txt" --dir "/scratch/htc/fchegini/cells830_corner_edge/output" --matlab_dir "/scratch/htc/fchegini/cells830_corner_edge/matlab_dir"  --interfacetypes 3

# edge + face
# srun -B *:*:* -n1 ./emiModel --input "./input/pepe_combi_domi.vtu" --extra_set "./input/pepe_combi_domi_extracellular.txt" --intra_set "./input/pepe_combi_domi_intracellular.txt" --excited "./input/pepe_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells20_pepe_edge_face/output" --matlab_dir "/scratch/htc/fchegini/cells20_pepe_edge_face/matlab_dir" --interfacetypes 6
# srun -B *:*:* -n1 ./emiModel --input "./input/rizzo_combi_domi.vtu" --extra_set "./input/rizzo_combi_domi_extracellular.txt" --intra_set "./input/rizzo_combi_domi_intracellular.txt" --excited "./input/rizzo_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells40_rizzo_edge_face/output" --matlab_dir "/scratch/htc/fchegini/cells40_rizzo_edge_face/matlab_dir" --interfacetypes 6
# srun -B *:*:* -n1 ./emiModel --input "./input/robin_mesh.vtu" --extra_set "./input/robin_extracellular.txt" --intra_set "./input/robin_intracellular.txt" --excited "./input/robin_early_excited.txt" --dir "/scratch/htc/fchegini/cells88_edge_face/output" --matlab_dir "/scratch/htc/fchegini/cells88_edge_face/matlab_dir"  --interfacetypes 6
# srun -B *:*:* -n1 ./emiModel --input "./input/kermit_mesh.vtu" --extra_set "./input/kermit_extracellular.txt" --intra_set "./input/kermit_intracellular.txt" --excited "./input/kermit_early_excited.txt" --dir "/scratch/htc/fchegini/cells172_edge_face/output" --matlab_dir "/scratch/htc/fchegini/cells172_edge_face/matlab_dir"  --interfacetypes 6
# srun -B *:*:* -n1 ./emiModel --input "./input/gonzo_mesh.vtu" --extra_set "./input/gonzo_extracellular.txt" --intra_set "./input/gonzo_intracellular.txt" --excited "./input/gonzo_early_excited.txt" --dir "/scratch/htc/fchegini/cells416_edge_face/output" --matlab_dir "/scratch/htc/fchegini/cells416_edge_face/matlab_dir"  --interfacetypes 6
# srun -B *:*:* -n1 ./emiModel --input "./input/animal_mesh.vtu" --extra_set "./input/animal_extracellular.txt" --intra_set "./input/animal_intracellular.txt" --excited "./input/animal_early_excited.txt" --dir "/scratch/htc/fchegini/cells830_edge_face/output" --matlab_dir "/scratch/htc/fchegini/cells830_edge_face/matlab_dir"  --interfacetypes 6


#corner + face

#srun -B *:*:* -n1 ./emiModel --input "./input/pepe_combi_domi.vtu" --extra_set "./input/pepe_combi_domi_extracellular.txt" --intra_set "./input/pepe_combi_domi_intracellular.txt" --excited "./input/pepe_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells20_pepe_corner_face/output" --matlab_dir "/scratch/htc/fchegini/cells20_pepe_corner_face/matlab_dir" --interfacetypes 5
#srun -B *:*:* -n1 ./emiModel --input "./input/rizzo_combi_domi.vtu" --extra_set "./input/rizzo_combi_domi_extracellular.txt" --intra_set "./input/rizzo_combi_domi_intracellular.txt" --excited "./input/rizzo_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells40_rizzo_corner_face/output" --matlab_dir "/scratch/htc/fchegini/cells40_rizzo_corner_face/matlab_dir" --interfacetypes 5
# srun -B *:*:* -n1 ./emiModel --input "./input/robin_mesh.vtu" --extra_set "./input/robin_extracellular.txt" --intra_set "./input/robin_intracellular.txt" --excited "./input/robin_early_excited.txt" --dir "/scratch/htc/fchegini/cells88_corner_face/output" --matlab_dir "/scratch/htc/fchegini/cells88_corner_face/matlab_dir"  --interfacetypes 5
# srun -B *:*:* -n1 ./emiModel --input "./input/kermit_mesh.vtu" --extra_set "./input/kermit_extracellular.txt" --intra_set "./input/kermit_intracellular.txt" --excited "./input/kermit_early_excited.txt" --dir "/scratch/htc/fchegini/cells172_corner_face/output" --matlab_dir "/scratch/htc/fchegini/cells172_corner_face/matlab_dir"  --interfacetypes 5
#srun -B *:*:* -n1 ./emiModel --input "./input/gonzo_mesh.vtu" --extra_set "./input/gonzo_extracellular.txt" --intra_set "./input/gonzo_intracellular.txt" --excited "./input/gonzo_early_excited.txt" --dir "/scratch/htc/fchegini/cells416_corner_face/output" --matlab_dir "/scratch/htc/fchegini/cells416_corner_face/matlab_dir"  --interfacetypes 5
# srun -B *:*:* -n1 ./emiModel --input "./input/animal_mesh.vtu" --extra_set "./input/animal_extracellular.txt" --intra_set "./input/animal_intracellular.txt" --excited "./input/animal_early_excited.txt" --dir "/scratch/htc/fchegini/cells830_corner_face/output" --matlab_dir "/scratch/htc/fchegini/cells830_corner_face/matlab_dir"  --interfacetypes 5
#

# srun -n1 ./emiModel --input "./input/pepe_combi_domi.vtu" --extra_set "./input/pepe_combi_domi_extracellular.txt" --intra_set "./input/pepe_combi_domi_intracellular.txt" --excited "./input/pepe_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells20_pepe_corner_face/output" --matlab_dir "/scratch/htc/fchegini/cells20_pepe_corner_face/matlab_dir" --interfacetypes 6
# srun -B *:*:* -n1 ./emiModel --input "./input/rizzo_combi_domi.vtu" --extra_set "./input/rizzo_combi_domi_extracellular.txt" --intra_set "./input/rizzo_combi_domi_intracellular.txt" --excited "./input/rizzo_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/cells40_rizzo_corner_face/output" --matlab_dir "/scratch/htc/fchegini/cells40_rizzo_corner_face/matlab_dir" --interfacetypes 6
# srun -B *:*:* -n1 ./emiModel --input "./input/robin_mesh.vtu" --extra_set "./input/robin_extracellular.txt" --intra_set "./input/robin_intracellular.txt" --excited "./input/robin_early_excited.txt" --dir "/scratch/htc/fchegini/cells88_corner_face/output" --matlab_dir "/scratch/htc/fchegini/cells88_corner_face/matlab_dir"  --interfacetypes 6
# srun -B *:*:* -n1 ./emiModel --input "./input/kermit_mesh.vtu" --extra_set "./input/kermit_extracellular.txt" --intra_set "./input/kermit_intracellular.txt" --excited "./input/kermit_early_excited.txt" --dir "/scratch/htc/fchegini/cells172_corner_face/output" --matlab_dir "/scratch/htc/fchegini/cells172_corner_face/matlab_dir"  --interfacetypes 6
# srun -B *:*:* -n1 ./emiModel --input "./input/gonzo_mesh.vtu" --extra_set "./input/gonzo_extracellular.txt" --intra_set "./input/gonzo_intracellular.txt" --excited "./input/gonzo_early_excited.txt" --dir "/scratch/htc/fchegini/cells416_corner_face/output" --matlab_dir "/scratch/htc/fchegini/cells416_corner_face/matlab_dir"  --interfacetypes 6
# srun -B *:*:* -n1 ./emiModel --input "./input/animal_mesh.vtu" --extra_set "./input/animal_extracellular.txt" --intra_set "./input/animal_intracellular.txt" --excited "./input/animal_early_excited.txt" --dir "/scratch/htc/fchegini/cells830_corner_face/output" --matlab_dir "/scratch/htc/fchegini/cells830_corner_face/matlab_dir"  --interfacetypes 6
#

# LNS 
#srun -B *:*:* -n1 ./emiModel --input "./input/animal_mesh.vtu" --extra_set "./input/animal_extracellular.txt" --intra_set "./input/animal_intracellular.txt" --excited "./input/animal_early_excited.txt" --dir "/scratch/htc/fchegini/LNS/output" --matlab_dir "/scratch/htc/fchegini/LNS/matlab_dir"  --T_ 20000
#srun -B *:*:* -n1 ./emiModel --input "./input/kermit_mesh.vtu" --extra_set "./input/kermit_extracellular.txt" --intra_set "./input/kermit_intracellular.txt" --excited "./input/kermit_early_excited.txt" --dir "/scratch/htc/fchegini/LNS_kermit/output" --matlab_dir "/scratch/htc/fchegini/LNS_kermit/matlab_dir"  --T_ 20000
# srun -B *:*:* -n1 ./emiModel --input "./input/robin_mesh.vtu" --extra_set "./input/robin_extracellular.txt" --intra_set "./input/robin_intracellular.txt" --excited "./input/robin_early_excited.txt" --dir "/scratch/htc/fchegini/LNS_robin/output" --matlab_dir "/scratch/htc/fchegini/LNS_robin/matlab_dir"  --T_ 20000
#----------------------------------------------------------------------------------------------------------------------------------------
# srun -B *:*:* -n1 ./emiModel --input "./input/pepe_combi_domi.vtu" --dir "./output_rescaled0" --matlab_dir "./matlab_dir_rescaled0"
# srun -B *:*:* -n1 ./emiModel --input "./input/pepe_combi_domi_smaller.vtu" --dir "./output_rescaled1" --matlab_dir "./matlab_dir_rescaled1"
# srun -B *:*:* -n1 ./emiModel --input "./input/pepe_combi_domi_smaller_more.vtu" --dir "./output_rescaled2" --matlab_dir "./matlab_dir_rescaled2"

#srun -B *:*:* -n1 ./emiModel --write_to_file false --input "./input/robin_combi_domi.vtu" --extra_set "./input/robin_combi_domi_extracellular.txt" --intra_set "./input/robin_combi_domi_intracellular.txt" --excited "./input/robin_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/robin/robin_combi_output" --matlab_dir "/scratch/htc/fchegini/robin/robin_combi_matlab_dir"  
#srun -B *:*:* -n1 ./emiModel --write_to_file false --input "./input/robin_combi_domi_smaller.vtu" --extra_set "./input/robin_combi_domi_extracellular.txt" --intra_set "./input/robin_combi_domi_intracellular.txt" --excited "./input/robin_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/robin/robin_combi_smaller_output" --matlab_dir "/scratch/htc/fchegini/robin/robin_combi_smaller_matlab_dir"  
#srun -B *:*:* -n1 ./emiModel --write_to_file false --input "./input/robin_sep_domi.vtu" --extra_set "./input/robin_sep_domi_extracellular.txt" --intra_set "./input/robin_sep_domi_intracellular.txt" --excited "./input/robin_sep_domi_excited.txt" --dir "/scratch/htc/fchegini/robin/robin_sep_output" --matlab_dir "/scratch/htc/fchegini/robin/robin_sep_matlab_dir"  
#srun -B *:*:* -n1 ./emiModel --write_to_file false --input "./input/robin_sep_domi_smaller.vtu" --extra_set "./input/robin_sep_domi_extracellular.txt" --intra_set "./input/robin_sep_domi_intracellular.txt" --excited "./input/robin_sep_domi_excited.txt" --dir "/scratch/htc/fchegini/robin/robin_sep_smaller_output" --matlab_dir "/scratch/htc/fchegini/robin/robin_sep_smaller_matlab_dir"  


#srun -B *:*:* -n1 ./emiModel --write_to_file false --input "./input/robin_mesh.vtu" --extra_set "./input/robin_extracellular.txt" --intra_set "./input/robin_intracellular.txt" --excited "./input/robin_early_excited.txt" --dir "/scratch/htc/fchegini/robin/robin_output" --matlab_dir "/scratch/htc/fchegini/robin/robin_matlab_dir"  

#srun -B *:*:* -n1 ./emiModel --write_to_file false --input "./input/kermit_mesh.vtu" --extra_set "./input/kermit_extracellular.txt" --intra_set "./input/kermit_intracellular.txt" --excited "./input/kermit_early_excited.txt" --dir "/scratch/htc/fchegini/kermit/kermit_output" --matlab_dir "/scratch/htc/fchegini/kermit/kermit_matlab_dir" 

#srun -B *:*:* -n1 ./emiModel --write_to_file false --input "./input/gonzo_mesh.vtu" --extra_set "./input/gonzo_extracellular.txt" --intra_set "./input/gonzo_intracellular.txt" --excited "./input/gonzo_early_excited.txt" --dir "/scratch/htc/fchegini/gonzo/gonzo_output" --matlab_dir "/scratch/htc/fchegini/gonzo/gonzo_matlab_dir"

#srun -B *:*:* -n1 ./emiModel --write_to_file false --input "./input/animal_mesh.vtu" --extra_set "./input/animal_extracellular.txt" --intra_set "./input/animal_intracellular.txt" --excited "./input/animal_early_excited.txt" --dir "/scratch/htc/fchegini/animal/animal_output" --matlab_dir "/scratch/htc/fchegini/animal/animal_matlab_dir"

# srun -B *:*:* -n1 ./emiModel --write_to_file false --input "./input/rizzo_combi_domi.vtu" --extra_set "./input/rizzo_combi_domi_extracellular.txt" --intra_set "./input/rizzo_combi_domi_intracellular.txt" --excited "./input/rizzo_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/rizzo/rizzo_output" --matlab_dir "/scratch/htc/fchegini/rizzo/rizzo_matlab_dir"
# srun -B *:*:* -n1 ./emiModel --write_to_file false --input "./input/rizzo_combi_domi.vtu" --extra_set "./input/rizzo_combi_domi_extracellular.txt" --intra_set "./input/rizzo_combi_domi_intracellular.txt" --excited "./input/rizzo_combi_domi_excited.txt" --dir "/scratch/htc/fchegini/rizzo/rizzo_output" --matlab_dir "/scratch/htc/fchegini/rizzo/rizzo_matlab_dir"







