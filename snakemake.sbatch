#!/bin/bash

#SBATCH -J snakeflow
#SBATCH -p batch
#SBATCH -n 8

# activate conda env

source ~/.bashrc; conda activate variant

# run jobs

SLURM_ARGS="-p {cluster.partition} -N {cluster.nodes} -n {cluster.ntasks} -c {cluster.ncpus} -J {cluster.job-name} -o {cluster.output} -e {cluster.error}

snakemake -j 999 --use-conda --cluster-config cluster.yaml --cluster "sbatch $SLURM_ARGS"

