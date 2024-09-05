The executab;es used to run the simulations are swift_mpi for 4 node run of SWIFT, swift for the single node run of SWIFT, and Arepo for the 4 node MPI run of Arepo.
walclock.py is the python program used to get the wallclock times for the timings plot in analyse_tasks.ipynb and wallclock.txt is that output.
parameters.yml are the parameters used for both of the SWIFT runs
param_c7.txt is the parameter file used for the Arepo run
The sbatch scripts used to submit these simulation runs to the compute nodes on cosma are the ones starting with submit_
