#!/bin/bash -l
 
#SBATCH -A snic2022-5-83
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 2-00:00:00
#SBATCH -J cmac_star_rsem_run_adults
#SBATCH -C usage_mail
#SBATCH --mail-user=axel.wiberg@ebc.uu.se
#SBATCH --mail-type=ALL

# TMP Directory
#$SNIC_TMP

#JOB ID
echo "JOB: ${SLURM_JOB_ID}"
echo ""

ml load bioinfo-tools
ml load samtools/1.14
ml load star/2.7.9a
ml load rsem/1.3.3

sample=${1}

# Move to temp dir

cd ${SNIC_TMP}

# Copy over reference and index
cp -r /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/pt_036_001_hifiasm_20201223.primary.fasta ${SNIC_TMP}/
cp -r /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/star_index ${SNIC_TMP}/
#cp -r ${HOME}/projects/Cmac_selection_lines_gene_expression/star_index ${SNIC_TMP}/
cp /proj/snic2021-6-30/delivery04381/INBOX/pt_036/analysis/pt_036_001/rsem_genome_dir/* ${SNIC_TMP}/
#cp ${HOME}/projects/Cmac_selection_lines_gene_expression/rsem_genome_dir/* ${SNIC_TMP}/

# Copy over sample read files
cp /proj/snic2021-6-30/selection_lines_gene_expression/adults/reads/rrna_removed/${sample}*mrna_fwd.fq.gz ${SNIC_TMP}/${sample}_R1.fq.gz
cp /proj/snic2021-6-30/selection_lines_gene_expression/adults/reads/rrna_removed/${sample}*mrna_rev.fq.gz ${SNIC_TMP}/${sample}_R2.fq.gz

# Map with STAR
echo "Run STAR..."

star  --runThreadN 5 \
--runMode alignReads \
--quantMode TranscriptomeSAM \
--genomeDir ${SNIC_TMP}/star_index \
--readFilesIn ${SNIC_TMP}/${sample}_R1.fq.gz ${SNIC_TMP}/${sample}_R2.fq.gz \
--readFilesCommand zcat \
--outFilterScoreMinOverLread 0 \
--outFilterMatchNminOverLread 0.5 \
--outFilterMatchNmin 0 \
--alignEndsType EndToEnd \
--outFileNamePrefix ${SNIC_TMP}/${sample} \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx

echo "Done"
ls -lh ${SNIC_TMP}/

echo "${sample}Aligned.sortedByCoord.out.bam"
samtools index ${SNIC_TMP}/${sample}Aligned.sortedByCoord.out.bam
samtools flagstat ${SNIC_TMP}/${sample}Aligned.sortedByCoord.out.bam > ${SNIC_TMP}/${sample}Aligned.sortedByCoord.out.flagstats
samtools stat ${SNIC_TMP}/${sample}Aligned.sortedByCoord.out.bam > ${SNIC_TMP}/${sample}Aligned.sortedByCoord.out.stats
echo "${sample}Aligned.toTranscriptome.out.bam"
samtools flagstat ${SNIC_TMP}/${sample}Aligned.toTranscriptome.out.bam > ${SNIC_TMP}/${sample}Aligned.toTranscriptome.out.flagstats
samtools stats ${SNIC_TMP}/${sample}Aligned.toTranscriptome.out.bam > ${SNIC_TMP}/${sample}Aligned.toTranscriptome.out.stats

#echo "Run rsem on output..."

rsem-calculate-expression --paired-end --alignments -p 5 ${SNIC_TMP}/${sample}Aligned.toTranscriptome.out.bam ${SNIC_TMP}/rsem_genome ${sample}_star_rsem

#echo "Done"

ls -lh ${SNIC_TMP}/
ls -lh ${SNIC_TMP}/${sample}_star_rsem*

# Copy the output back
cp ${SNIC_TMP}/${sample}Aligned.sortedByCoord.out.bam ${HOME}/projects/Cmac_selection_lines_gene_expression/adults/star_rsem/
cp ${SNIC_TMP}/${sample}Aligned.sortedByCoord.out.bam.bai ${HOME}/projects/Cmac_selection_lines_gene_expression/adults/star_rsem/
cp ${SNIC_TMP}/${sample}Aligned*flagstats ${HOME}/projects/Cmac_selection_lines_gene_expression/adults/star_rsem/
cp ${SNIC_TMP}/${sample}Aligned*.stats ${HOME}/projects/Cmac_selection_lines_gene_expression/adults/star_rsem/

cp ${SNIC_TMP}/${sample}_star_rsem.genes.results ${HOME}/projects/Cmac_selection_lines_gene_expression/adults/star_rsem/
cp ${SNIC_TMP}/${sample}_star_rsem.isoforms.results ${HOME}/projects/Cmac_selection_lines_gene_expression/adults/star_rsem/
cp ${SNIC_TMP}/${sample}_star_rsem.stat ${HOME}/projects/Cmac_selection_lines_gene_expression/adults/star_rsem/
cp ${SNIC_TMP}/${sample}Log.final.out ${HOME}/projects/Cmac_selection_lines_gene_expression/adults/star_rsem/

#cp ${SNIC_TMP}/*Unmapped.out.mate1 ${HOME}/projects/Cmac_selection_lines_gene_expression/adults/star_rsem/


