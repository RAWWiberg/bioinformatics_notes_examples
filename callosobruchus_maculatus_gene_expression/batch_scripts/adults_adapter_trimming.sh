#!/bin/bash -l
 
#SBATCH -A snic2022-5-83
#SBATCH -p core
#SBATCH -n 5
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

sample=${1}

# move to temp dir
cd ${SNIC_TMP}

# create output directories
mkdir ${SNIC_TMP}/adtr

# copy over raw reads

#cp /proj/snic2021-6-30/selection_lines_gene_expression/reads/rrna_removed/${sample}* ${SNIC_TMP}/

cp /proj/snic2021-6-30/delivery06076/INBOX/VC-3247/220617_A00181_0525_AHMGWYDSX3/Sample_VC-3247-${sample}/*R1_001.fastq.gz ${SNIC_TMP}/${sample}_R1.fq.gz 
cp /proj/snic2021-6-30/delivery06076/INBOX/VC-3247/220617_A00181_0525_AHMGWYDSX3/Sample_VC-3247-${sample}/*R2_001.fastq.gz ${SNIC_TMP}/${sample}_R2.fq.gz

# run cutadapt
echo "Running cutadapt..."
R1=${SNIC_TMP}/${sample}_R1.fq.gz
R2=${SNIC_TMP}/${sample}_R2.fq.gz

cutadapt \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
--cores 5 \
-m 75 \
-o ${SNIC_TMP}/adtr/${sample}_adTr_R1.fastq.gz \
-p ${SNIC_TMP}/adtr/${sample}_adTr_R2.fastq.gz \
${R1} ${R2}

# copy the output back
cp ${SNIC_TMP}/adtr/${sample}_adTr* /proj/snic2021-6-30/selection_lines_gene_expression/adults/reads/adtr/

echo "Done"

# Run fastqc
echo ""
echo "Running fastqc"
fastqc -t 5 --noextract ${SNIC_TMP}/adtr/${sample}_adTr_R1.fastq.gz ${SNIC_TMP}/adtr/${sample}_adTr_R2.fastq.gz

# copy the output back 
cp ${SNIC_TMP}/adtr/${sample}*.zip /proj/snic2021-6-30/selection_lines_gene_expression/adults/reads/adtr/
cp ${SNIC_TMP}/adtr/${sample}*.html /proj/snic2021-6-30/selection_lines_gene_expression/adults/reads/adtr/

echo "Done"

