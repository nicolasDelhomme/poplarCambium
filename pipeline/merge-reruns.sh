#!/bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nathaniel.street@umu.se

## stop on error
set -ex



module load bioinfo-tools samtools

cd /mnt/picea/projects/aspseq/24_Bellini_cambium/BelliniCambium/star/

samtools merge T89-2_merged.bam FCC4LJ7ACXX-TRDPOPuljTAMRAAPEI-213*_STAR.bam T89-2_FCC4CL0ACXX_L2_TRDPOPuljTAMRAAPEI-213*_STAR.bam
samtools merge T89-3_merged.bam FCC4LJ7ACXX-TRDPOPuljTAJRAAPEI-210*_STAR.bam T89-3_FCC4CL0ACXX_L2_TRDPOPuljTANRAAPEI-214*_STAR.bam
samtools merge T89-4_merged.bam FCC4LJ7ACXX-TRDPOPuljTAORAAPEI-215*_STAR.bam T89-4_FCC4CL0ACXX_L2_TRDPOPuljTAORAAPEI-215*.bam

samtools merge OP42-1_merged.bam FCC4LJ7ACXX-TRDPOPuljTAIRAAPEI-209*_STAR.bam OP1_FCC4CL0ACXX_L2_TRDPOPuljTAIRAAPEI-209*_STAR.bam
samtools merge OP42-5_merged.bam FCC4LJ7ACXX-TRDPOPuljTALRAAPEI-212*_STAR.bam OP5_FCC4CL0ACXX_L2_TRDPOPuljTALRAAPEI-212*_STAR.bam

