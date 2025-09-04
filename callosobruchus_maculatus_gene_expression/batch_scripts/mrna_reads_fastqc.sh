#!/bin/bash -l
 
#SBATCH -A snic2022-5-83
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 4-00:00:00
#SBATCH -J adapter_trimming
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
ml load

sample=${1}

# move to temp dir
cd ${SNIC_TMP}

# Copy over new raw reads

cp /proj/snic2021-6-30/selection_lines_gene_expression/reads/rrna_removed/${sample}/*R1.fastq.gz ${SNIC_TMP}/${sample}_R1.fq.gz 
cp /proj/snic2021-6-30/selection_lines_gene_expression/reads/rrna_removed/${sample}/*R2.fastq.gz ${SNIC_TMP}/${sample}_R2.fq.gz 

# Run fastqc
echo ""
echo "Running fastqc"
fastqc -t 10 --noextract ${SNIC_TMP}/adtr/${sample}_merge_adTr_R1.fastq.gz ${SNIC_TMP}/adtr/${sample}_merge_adTr_R2.fastq.gz

# copy the output back 
cp ${SNIC_TMP}/adtr/${sample}*.zip /proj/snic2021-6-30/selection_lines_gene_expression/reads/adtr/
cp ${SNIC_TMP}/adtr/${sample}*.html /proj/snic2021-6-30/selection_lines_gene_expression/reads/adtr/

echo "Done"
