#!/bin/bash
#SBATCH --job-name=anotaciones_create
#SBATCH --chdir=/mnt/data/bioinfo3/ysosa/proyecto-equipo-A/outs
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=anotacion-%j.out
#SBATCH --error=anotacion-%j.err

module load r/4.6.0

Rscript ../R/anotacion.R
