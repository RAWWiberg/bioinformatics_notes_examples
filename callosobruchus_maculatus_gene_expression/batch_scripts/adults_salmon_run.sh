#!/bin/bash -l
 
#SBATCH -A snic2022-5-83
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 2-00:00:00
#SBATCH -J cmac_adults_salmon_run
#SBATCH -C usage_mail
#SBATCH --mail-user=axel.wiberg@ebc.uu.se
#SBATCH --mail-type=ALL

# TMP Directory
#$SNIC_TMP

#JOB ID
echo "JOB: ${SLURM_JOB_ID}"
echo ""

ml load bioinfo-tools
ml load Salmon/1.6.0

sample=${1}

# Move to temp dir

cd ${SNIC_TMP}

# Copy over reference and index
cp -r ${HOME}/projects/Cmac_selection_lines_gene_expression/p1.GS.annotation.v1.transcripts* ${SNIC_TMP}/
cp ${HOME}/projects/Cmac_selection_lines_gene_expression/p1.GS.annotation.v1_trans2gene.tab ${SNIC_TMP}/

# Copy over sample read files
cp /proj/snic2021-6-30/selection_lines_gene_expression/adults/reads/adtr/${sample}_adTr_R1.fastq.gz ${SNIC_TMP}/${sample}_R1.fq.gz
cp /proj/snic2021-6-30/selection_lines_gene_expression/adults/reads/adtr/${sample}_adTr_R2.fastq.gz ${SNIC_TMP}/${sample}_R2.fq.gz

echo "Run salmon..."

salmon quant \
-l A \
-i ${SNIC_TMP}/p1.GS.annotation.v1.transcripts \
-g ${SNIC_TMP}/p1.GS.annotation.v1_trans2gene.tab \
-o ${SNIC_TMP}/${sample}_salmon_quant \
-1 ${SNIC_TMP}/${sample}_R1.fq.gz \
-2 ${SNIC_TMP}/${sample}_R2.fq.gz  \
--seqBias --gcBias --posBias -p 10

echo "Done"

ls -lh ${SNIC_TMP}/

# Copy the output back
cp -r ${SNIC_TMP}/${sample}_salmon_quant ${HOME}/projects/Cmac_selection_lines_gene_expression/adults/salmon/
