#!/bin/bash
#SBATCH --time=14:00:00
#SBATCH --mem=2G

srun julia /scratch/work/condeil1/EnergySystemModeling.jl/aggregation/form_aggreg_instance.jl