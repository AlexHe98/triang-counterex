#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24  # need extra cpu to run manager
#SBATCH --mem=96G
#SBATCH --time=1-0:00  # time (D-HH:MM)
#SBATCH -o ../output/remTunTrefoilP23M96G.%j.%N.out  # STDOUT
#SBATCH -e ../output/remTunTrefoilP23M96G.%j.%N.out  # STDERR

set -e

# Use 23 processes, and wait 1 hour to print updates.
~/software/bin/regina-python removeTunnel.py cPcbbbadu 23 3600

