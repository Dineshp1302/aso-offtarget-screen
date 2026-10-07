#!/bin/bash
#SBATCH --job-name=aso_prod
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/scratch/%u/aso/logs/prod_%j.out

module load bowtie2 samtools
cd /scratch/$USER/aso

bowtie2 -f -x ref/GRCh38_noalt_as -U results/candidates.fa \
  --end-to-end -a -N 1 -L 10 -i C,1,0 \
  --mp 6,6 --rdg 100,100 --rfg 100,100 \
  --score-min L,-18,0 --no-unal -p 16 \
| samtools sort -@ 4 -o results/hits_genome.bam -

samtools index results/hits_genome.bam
samtools view -c results/hits_genome.bam
echo "DONE"
