# *M. cliftonense* Adult vs. Hatchling contrasts.

### Table of contents
0. [Install packages](#0packages)
1. [Trim and filter reads](#1reads_prep)
2. [Map/Quantify](#2readmapquant)

This document details example code used for the processing of data from the Positional RNA-Seq experiments.
For clarity, I have ommitted, in most cases, full paths to input and output files.
Data were processed in the same way for all three species. I therefore only show one example.

### 0. Install packages, make environments, obtain docker containers <a name="0packages"></a>

	conda create -n salmon
	conda activate salmon
	conda install salmon
	conda deactivate
	conda list --export -n salmon > salmon_env.txt


### 1. Trim and filter reads in prep for mapping <a name="1reads_prep"><\a>

Trim TruSeq adapters from all reads.  
use the Illumina TruSeq adapters: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`

This was done for many individual read files. 
Here I only give example code for a single instance.

	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  --minimum-length 50 --nextseq-trim=25 --cores 20 \
	-o ${READ}_R1.adTr.fastq -p ${READ}_R2.adTr.fastq \
	${READ}_R1_001_MM_1.fastq ${READ}_R2_001_MM_1.fastq	


### 2. Mapping + quantification <a name="2readmapquant"></a>

	conda activate salmon

For each species, the reference annotation and a gene to transcript map needs to be provided to salmon.
In this document I use the placeholders `${reference_annotation` and `${ref_gene2trans_map}` for these files.
I use the placeholder `${spp}` for the species name

	trimDir=${HOME}/data/mac_genomes/genome_Maccli/canu_assembly/final_annotation/AvH_contrast/trimmed_reads

Make the salmon index

	salmon index -t ${reference_annotation} -i ./${spp}_index.cds

Then we simply run `salmon quant` with the input reads to perform pseudo-mapping and quantification at the same time.

	for i in ${samples}; do
	    read1=$(ls ${trimDir}/*${i}*_R1.adTr.fastq)
	    read2=$(ls ${trimDir}/*${i}*_R2.adTr.fastq)
	    echo "Running salmon for: ${i} - $(basename ${read1}) $(basename ${read2})"
	    salmon quant -l A -p 30 --seqBias --gcBias --posBias -g ${reference_dir}/Maccli_v2_rnm_trans2gn_slm.txt -i ./Maccli_v2.cds -1 ${read1} -2 ${read2} --output ${i}; 
	done

	conda deactivate