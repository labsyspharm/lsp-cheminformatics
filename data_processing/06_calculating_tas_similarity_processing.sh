#!/bin/bash
#SBATCH -t 0-2:00                        # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=2G                          # Memory total in MB (for all cores)

set -eux

Rscript 06_calculating_tas_similarity_processing.R "$1" "$2" "$3"
