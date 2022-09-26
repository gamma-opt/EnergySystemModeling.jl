#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --array=1-3360
#SBATCH --output=slurm_array_out/new/%A_%a.out

n=$SLURM_ARRAY_TASK_ID
instance=`sed -n "${n} p" instances1.txt`      # Get n-th line (1-indexed) of the file

srun julia /scratch/work/condeil1/EnergySystemModeling.jl/.triton/exe/opt/run_nodal_clust.jl ${instance} "" "false"
# srun julia /scratch/work/condeil1/EnergySystemModeling.jl/.triton/exe/opt/run_nodal_clust.jl ${instance} "dc" "false"
# srun julia /scratch/work/condeil1/EnergySystemModeling.jl/.triton/exe/opt/run_nodal_clust.jl ${instance} "dc_new" "false"
# srun julia /scratch/work/condeil1/EnergySystemModeling.jl/.triton/exe/opt/run_nodal_clust.jl ${instance} "dc_new_nosun" "true"
# srun julia /scratch/work/condeil1/EnergySystemModeling.jl/.triton/exe/opt/run_nodal_clust.jl ${instance} "" "true"

# srun julia /scratch/work/condeil1/EnergySystemModeling.jl/.triton/exe/opt/run_fix.jl ${instance}