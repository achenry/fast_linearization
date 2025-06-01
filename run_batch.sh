#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --time=04:00:00
#SBATCH --partition=shas-testing
#SBATCH --output=run_batch.out

module purge
module load python

source /curc/sw/anaconda3/latest
conda activate weis-env

python main.py

mv results /projects/aohe7145/OpenFAST-results/linearization