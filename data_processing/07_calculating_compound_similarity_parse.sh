#!/bin/bash
#SBATCH -c 2
#SBATCH -t 0-8:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=15G                          # Memory total in MB (for all cores)

source activate base

set -eux

python3 07_calculating_compound_similarity_parse.py "$1" - | \
  sort -k1,1 -k2,2 -S 10G --parallel=2 > "$2"
