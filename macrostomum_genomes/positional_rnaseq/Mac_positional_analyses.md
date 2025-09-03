# Positional RNA-seq data analyses

### Table of contents
0. [Install packages](#0packages)
1. [Reads trimming](#1trimming)
	1. [Adapter trimming](#1-1adapter)
2. [Read mapping + quantification](#2readmapquant)

This document details example code used for the processing of data from the Positional RNA-Seq experiments.
For clarity, I have ommitted, in most cases, full paths to input and output files.
Data were processed in the same way for all three species. I therefore only show one example.

### 0. Install packages, make environments, obtain docker containers <a name="0packages"></a>

	conda create -n salmon
	conda activate salmon
	conda install salmon
	conda deactivate
	conda list --export -n salmon > salmon_env.txt

### 1. Read trimming <a name="1trimming"></a>
#### 1.1 Adapter trimming <a name="1-1adapter"></a>

Trim TruSeq adapters from all reads.  
use the Illumina TruSeq adapters: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`

This was done for many individual read files. 
Here I only give example code for a single instance.

	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --minimum-length 50 --nextseq-trim=25 --cores 20 -o ${READ}_adTr.fastq.gz ${READ}_R1_001_MM_1.fastq.gz

Since fragment samples were sequenced on two lanes, I concatenate reads from the same fragment samples.

### 2. Read mapping + quantification <a name="2readmapquant"></a>

	conda activate salmon

For each species, the reference annotation and a gene to transcript map needs to be provided to salmon.
In this document I use the placeholders `${reference_annotation` and `${ref_gene2trans_map}` for these files.
I use the placeholder `${spp}` for the species name

	trimDir=${HOME}/data/positional_rnaseq/trimmed_reads

Make the salmon index

	salmon index -t ${reference_annotation} -i ./${spp}_index.cds

Map and quantify with salmon

	fragments=$(ls ${trimDir}/*fq.gz | sed 's;.fq.gz;;g')
	for fragment in $(echo ${fragments}); do \
		fragment=$(basename ${fragment}) && \
		salmon quant -p 20 --index ${spp}_index --libType A -r ${trimDir}/${fragment}.fq.gz --geneMap ${ref_gene2trans_map} --output ${fragment}_readsmapped.out; \
	done

	conda deactivate