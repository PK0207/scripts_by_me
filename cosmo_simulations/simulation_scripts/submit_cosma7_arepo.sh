#!/bin/bash -l

#SBATCH --ntasks 112
#SBATCH -J volume1_arepo #job name
#SBATCH -o std_output.%J.out
#SBATCH -e std_error.%J.err
#SBATCH -p cosma7
#SBATCH -A dp004
#SBATCH --exclusive
#SBATCH -t 72:00:00
#SBATCH --mail-type=END 
#SBATCH --mail-user=prerana.kottapalli@durham.ac.uk


module purge
module load intel_comp/2018
module load intel_mpi/2018
module load parallel_hdf5/1.8.20
module load gsl/2.4
module load fftw/3.3.7

#cd /cosma7/data/dp004/azadehf/APOSTLE_AURIGA/apostle_MR/S5_DMO/run

# Run the program
mpirun  -np $SLURM_NTASKS  ./Arepo param_c7.txt 

 

