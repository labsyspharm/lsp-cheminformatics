#!/bin/bash
#SBATCH -t 0-4:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=8G                          # Memory total in MB (for all cores)

source activate chemfp

set -eux

if [ "$1" == "$2" ]; then
  simsearch -t 0.2 --NxN --memory -k all -o "$3" "$1"
else
  simsearch -t 0.2 --memory -k all -o "$3" -q "$2" "$1"
fi
