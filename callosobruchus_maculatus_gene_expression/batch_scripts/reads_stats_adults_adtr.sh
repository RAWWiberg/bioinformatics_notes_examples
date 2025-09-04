#!/bin/bash -l
 
#SBATCH -A snic2022-5-83
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 2-00:00:00
#SBATCH -J cmac_reads_stats_adults
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

# move to reads dir
cd /proj/snic2021-6-30/selection_lines_gene_expression/adults/reads/adtr

ls -lh ./*fastq.gz

echo "Running seqkit stats..."

seqkit stats ./*.fastq.gz -j 10 -b -a -T > ./adtr_fastq_stats.tab

echo "Done"
