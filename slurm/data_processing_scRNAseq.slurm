#!/bin/bash
#SBATCH --job-name="data_processing_scRNAseq"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=5-00:00:00
#SBATCH --mem=16gb
#SBATCH --gres=gpu:1
#SBATCH --output="data_processing_scRNAseq-%j.o"
#SBATCH --error="data_processing_scRNAseq-%j.e"
#SBATCH --mail-user=zhaok220@wfu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --account="generalGrp"
#SBATCH --partition=gpu

source /deac/csc/khuriGrp/software/zhaok220/bin/activate

cd /deac/csc/khuriGrp/zhaok220/data_processing_scRNA/scripts/
python3 data_processing_scRNAseq.py

slurm_mem_report

exit 0

