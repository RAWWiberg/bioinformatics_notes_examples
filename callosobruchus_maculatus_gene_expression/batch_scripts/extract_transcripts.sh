#!/bin/bash -l
 
#SBATCH -A snic2022-5-83
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 2-00:00:00
#SBATCH -J cmac_extract_transcripts
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

# Copy over annotation and reference genome

cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/pt_036_001_hifiasm_20201223.primary.fasta ./

cp /proj/snic2021-6-30/douglas/p1.GS.annotation.v1+Tor/p1.GS.annotation.v1+Tor.gtf ./

echo "Extracting transcripts..."

gffread ./p1.GS.annotation.v1+Tor.gtf -w ./p1.GS.annotation.v1.transcripts.fa -g ./pt_036_001_hifiasm_20201223.primary.fasta

echo "Done"

# Copy the files back
cp -r p1.GS.annotation.v1.transcripts.fa /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/
