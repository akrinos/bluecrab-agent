#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --qos=unlim
#SBATCH --mem=600000
#SBATCH --time=48:00:00
#SBATCH -n 8

julia code/bcmodel_runall.jl