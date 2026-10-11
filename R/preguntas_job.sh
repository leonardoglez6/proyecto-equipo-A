#!/bin/bash
#SBATCH --job-name=preguntas_4
#SBATCH --chdir=/mnt/data/bioinfo3/ysosa/proyecto-equipo-A/outs
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=preguntas-%j.out
#SBATCH --error=preguntas-%j.err

module load r/4.6.0

Rscript ../R/preguntas.R
