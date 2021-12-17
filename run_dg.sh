#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH -J dg
#SBATCH -o dg.%A.out
#SBATCH -e dg.%A.err

source /nas/longleaf/home/wwbooker/miniconda3/etc/profile.d/conda.sh
conda deactivate
conda activate pysam
python full_genotyper.py