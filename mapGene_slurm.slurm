#!/bin/bash
BATCH --time=04:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=1024       # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=mapGene
#SBATCH --error=/work/unobioinformatics/luppala/job.%J.err
#SBATCH --output=/work/unobioinformatics/luppala/job.%J.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=luppala@unomaha.edu

module load biopython

python mapProbes.py

hostname
sleep 60
