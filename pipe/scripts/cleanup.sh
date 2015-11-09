#!/bin/bash
#SBATCH -e logs/cleanup.err
#SBATCH -o logs/cleanup.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute
#SBATCH --mem=10
# Remove extra files
rm $1