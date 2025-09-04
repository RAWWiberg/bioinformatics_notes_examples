#!/bin/bash -l
 
#SBATCH -A snic2022-5-83
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 2-00:00:00
#SBATCH -J fastqc
#SBATCH -C usage_mail
#SBATCH --mail-user=axel.wiberg@ebc.uu.se
#SBATCH --mail-type=ALL

# TMP Directory
#$SNIC_TMP

#JOB ID
echo "JOB: ${SLURM_JOB_ID}"
echo ""

ml load bioinfo-tools
ml load cutadapt/4.0
ml load FastQC/0.11.9
ml load fastq_screen/0.11.1
ml load MultiQC/1.12

# move to temp dir
cd /proj/snic2021-6-30/selection_lines_gene_expression/reads/rrna_removed/

# Run fastqc
echo ""
echo "Running fastqc"

fastqc -t 20 --noextract *.fq.gz

multiqc ./

echo "Done"
