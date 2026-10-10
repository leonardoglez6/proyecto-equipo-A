#!/bin/bash
#SBATCH --job-name=MT_RP_dbl_rscript
#SBATCH --chdir=/mnt/data/bioinfo3/ysosa/proyecto-equipo-A/outs
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=QC_rscript-%j.out
#SBATCH --error=QC_rscript-%j.err

module load r/4.6.0

Rscript ../R/MT_RP_Dbl.R
