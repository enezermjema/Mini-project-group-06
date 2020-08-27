#!/bin/bash

#SBATCH -J snakeflow
#SBATCH -p batch
#SBATCH -n 4

# activate conda env

source ~/.bashrc; conda activate variant

# run jobs

#SLURM_ARGS="-p {cluster.partition} -N {cluster.nodes} -c {cluster.ncpus} -J {cluster.job-name} -o {cluster.output} -e {cluster.error}"

snakemake -j 666 --use-conda --cluster-config cluster.yaml --cluster 'sbatch -p {cluster.partition} -c {cluster.ncpus} -o {cluster.output} -e {cluster.error}'

