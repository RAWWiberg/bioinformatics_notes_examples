#!/bin/bash -l
 
#SBATCH -A snic2022-5-83
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 3-00:00:00
#SBATCH -J cmac_adults_rrna_removal
#SBATCH -C usage_mail
#SBATCH --mail-user=axel.wiberg@ebc.uu.se
#SBATCH --mail-type=ALL

# TMP Directory
#$SNIC_TMP

#JOB ID
echo "JOB: ${SLURM_JOB_ID}"
echo ""

ml load bioinfo-tools
ml load SortMeRNA/4.3.3

sample=${1}

echo "SAMPLE: ${sample}"

# move to temp dir
cd ${SNIC_TMP}

# create output directories
mkdir ${SNIC_TMP}/rrna_reads
mkdir ${SNIC_TMP}/reads	

# copy over raw reads

cp /proj/snic2021-6-30/selection_lines_gene_expression/adults/reads/adtr/${sample}_adTr_R1.fastq.gz ${SNIC_TMP}/${sample}_adTr_R1.fastq.gz
cp /proj/snic2021-6-30/selection_lines_gene_expression/adults/reads/adtr/${sample}_adTr_R2.fastq.gz ${SNIC_TMP}/${sample}_adTr_R2.fastq.gz

ls -lh ${SNIC_TMP}/

# Copy over reference ribosomal sequences

cp ${HOME}/projects/Cmac_selection_lines_gene_expression/rrna/ncbi_cdhit95.fa ${SNIC_TMP}/

echo "Running sortmerna..."


gunzip *fastq.gz

# run sortmerna
R1=${SNIC_TMP}/${sample}_adTr_R1.fastq
R2=${SNIC_TMP}/${sample}_adTr_R2.fastq

ls -lh ${SNIC_TMP}/

sortmerna -out2 -paired_in -v -threads 20 \
--workdir ${SNIC_TMP}/ \
--ref ${SNIC_TMP}/ncbi_cdhit95.fa \
-reads ${R1} -reads ${R2} \
-fastx \
--no-best \
--num_alignments 1 \
-aligned ${SNIC_TMP}/rrna_reads/${sample}_rrna \
-other ${SNIC_TMP}/reads/${sample}_mrna

echo "Done"
# copy the output back
# keep all the reads and the log file

cp ${SNIC_TMP}/rrna_reads/${sample}*.log /proj/snic2021-6-30/selection_lines_gene_expression/adults/reads/rrna_removed/

cp ${SNIC_TMP}/rrna_reads/${sample}_rrna* /proj/snic2021-6-30/selection_lines_gene_expression/adults/reads/rrna_removed/

cp ${SNIC_TMP}/reads/${sample}_mrna* /proj/snic2021-6-30/selection_lines_gene_expression/adults/reads/rrna_removed/
