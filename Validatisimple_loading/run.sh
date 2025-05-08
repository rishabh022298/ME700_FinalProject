#!/bin/bash -l
#$ -m ea
#$ -N Jint
#$ -e error
#$ -l h_rt=80:00:00
#$ -pe omp 4
module purge
module load fenics/2019.1.0
run_fenics.sh python3 JInt.py
