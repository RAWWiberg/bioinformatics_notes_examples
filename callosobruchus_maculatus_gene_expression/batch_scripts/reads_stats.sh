#!/bin/bash -l
 
#SBATCH -A snic2022-5-83
#SBATCH -p core
#SBATCH -n 15
#SBATCH -t 2-00:00:00
#SBATCH -J cmac_reads_stats
#SBATCH -C usage_mail
#SBATCH --mail-user=axel.wiberg@ebc.uu.se
#SBATCH --mail-type=ALL

# TMP Directory
#$SNIC_TMP

#JOB ID
echo "JOB: ${SLURM_JOB_ID}"
echo ""

ml load bioinfo-tools
ml load SeqKit/0.15.0

reads=${1}

echo "READS: ${reads}"

# move to reads dir
cd /proj/snic2021-6-30/selection_lines_gene_expression/reads/${reads}

ls -lh ./*fastq.gz

# Copy over reference ribosomal sequences


echo "Running seqkit stats..."

seqkit stats ./*.fastq.gz -j 10 -b -a -T > ./${reads}_fastq_stats.tab

echo "Done"
