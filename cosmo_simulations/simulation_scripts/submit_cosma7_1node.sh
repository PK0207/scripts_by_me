#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=28
#SBATCH -J volume1_swift_singlenode   #job name 
#SBATCH -o std_output.%J.out
#SBATCH -e std_error.%J.err
#SBATCH -p cosma7
#SBATCH -A dp004
#SBATCH --exclusive
#SBATCH -t 72:00:00
#SBATCH --mail-user=prerana.kottapalli@durham.ac.uk

module purge
module load intel_comp/2021.1.0 compiler
module load intel_mpi/2018
module load ucx/1.10.1
module load fftw/3.3.9
module load parallel_hdf5/1.10.6
module load parmetis/4.0.3-64bit
module load gsl/2.5

echo $SLURM_CPUS_PER_TASK

# Run the program
./swift -y 100 -Y 100 --threads=$SLURM_CPUS_PER_TASK --cosmology --self-gravity --pin parameter.yml

 

