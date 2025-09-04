#!/bin/bash -l
 
#SBATCH -A snic2022-5-83
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 2-00:00:00
#SBATCH -J salmon_index
#SBATCH -C usage_mail
#SBATCH --mail-user=axel.wiberg@ebc.uu.se
#SBATCH --mail-type=ALL

# TMP Directory
#$SNIC_TMP

#JOB ID
echo "JOB: ${SLURM_JOB_ID}"
echo ""

ml load bioinfo-tools
ml load cufflinks/2.2.1
ml load Salmon/1.6.0

# Move to temp dir

cd ${SNIC_TMP}

echo "Copying files over..."

cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/pt_036_001_hifiasm_20201223.primary.fasta ${SNIC_TMP}/

cp ${HOME}/projects/Cmac_selection_lines_gene_expression/rrna/ncbi_silva_coleoptera_cdhit95.fa ${SNIC_TMP}/

cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/p1.GS.annotation.v1.transcripts.fa ${SNIC_TMP}/

# Add decoy sequences to transcriptome

cat ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta ${SNIC_TMP}/ncbi_silva_coleoptera_cdhit95.fa | grep ">" | sed 's; .*;;g' > ${SNIC_TMP}/salmon_decoys.txt
sed -i 's;>;;' ${SNIC_TMP}/salmon_decoys.txt 

cat ${SNIC_TMP}/p1.GS.annotation.v1.transcripts.fa ${SNIC_TMP}/pt_036_001_hifiasm_20201223.primary.fasta ${SNIC_TMP}/ncbi_silva_coleoptera_cdhit95.fa > ${SNIC_TMP}/p1.GS.annotation.v1.transcripts_decome.fa

echo "Done"

echo "Creating salmon index..."

salmon index -p 10 -t p1.GS.annotation.v1.transcripts_decome.fa -i p1.GS.annotation.v1.transcripts_decome --decoys salmon_decoys.txt -k 31

echo "Done"

# Copy the files back
cp -r p1.GS.annotation.v1.transcripts_decome* /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/
cp ${SNIC_TMP}/salmon_decoys.txt /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/
