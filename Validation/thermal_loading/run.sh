#!/bin/bash -l
#$ -m ea
#$ -N thermal
#$ -e error
#$ -l h_rt=80:00:00
#$ -pe omp 16
module purge
module load fenics/2019.1.0
run_fenics.sh python3 thermal.py
