#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --output=slurm_out/new/%A.out

srun julia /scratch/work/condeil1/EnergySystemModeling.jl/.triton/exe/opt/run_FTR_nodal.jl