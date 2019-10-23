#!/bin/bash
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=8G                          # Memory total in MB (for all cores)

source activate chemfp

set -eux

Rscript 06_calculating_tas_similarity_processing.R "$1" "$2" "$3"
