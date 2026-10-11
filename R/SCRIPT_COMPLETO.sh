#!/bin/bash
#SBATCH --job-name=anotaciones_create
#SBATCH --chdir=/mnt/data/bioinfo3/ysosa/proyecto-equipo-A/outs
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=anotacion-%j.out
#SBATCH --error=anotacion-%j.err

# Crear carpeta para almacenar los dator
mkdir ../datos

# 1er script
sbatch preguntas_job.sh

# 2do script
sbatch eD_SX_manual.sh

# 3er script
sbatch eD_SX_auto.sh

# 4to script
sbatch anotacion_job.sh

# 5to script
sbatch MT_job.sh
