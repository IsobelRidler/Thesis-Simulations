#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-30
#SBATCH --job-name=lassoarrayT5e
#SBATCH --mem-per-cpu=12G
#SBATCH --time=12:00:00
#SBATCH --mail-user=isobel.ridler@kcl.ac.uk
#SBATCH --mail-type=FAIL
#SBATCH --requeue
#SBATCH --output=/scratch/users/k1620435/output/array_%A_%a.out

#Loading R
module load apps/R/3.6.0
#Calling the R script you want to run
Rscript /scratch/users/k1620435/scripts/lassoT5e.R
