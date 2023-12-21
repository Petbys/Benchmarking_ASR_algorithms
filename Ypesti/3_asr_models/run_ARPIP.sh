#!/bin/bash -l
#SBATCH -A uppmax2023-2-8
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 15:00:00
#SBATCH -J test_asr
#SBATCH --mail-type=All
#SBATCH --mail-user neh-eyong.langsi.0046@student.uu.se

module load bioinfo-tools
module load ARPIP

ARPIP params=conf_ypesti.txt
