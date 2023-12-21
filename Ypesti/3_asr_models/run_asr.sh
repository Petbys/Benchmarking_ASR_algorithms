#!/bin/bash -l
#SBATCH -A uppmax2023-2-8
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 15:00:00
#SBATCH -J test_asr
#SBATCH --mail-type=All
#SBATCH --mail-user petter.bystrom.8041@student.uu.se

module load conda
module load paml/4.10.7
pip install biopython
pip install Pandas
pip install PhyloPandas
pip install DendroPy
pip install pyasr
pip install toytree

python pyasr_test.py --fasta --tree --dir