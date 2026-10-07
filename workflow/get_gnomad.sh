#!/bin/bash
#SBATCH --job-name=gnomad_dl
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/scratch/%u/aso/logs/gnomad_dl_%j.out

cd /scratch/$USER/aso/ref
BASE=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes
wget -c $BASE/gnomad.genomes.v4.1.sites.chr18.vcf.bgz
wget -c $BASE/gnomad.genomes.v4.1.sites.chr18.vcf.bgz.tbi
echo "DONE"
ls -lh gnomad.genomes.v4.1.sites.chr18.vcf.bgz
