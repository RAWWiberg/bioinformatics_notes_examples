1. [Trim and filter reads](#1reads_prep)
2. [Map/Quantify](#2map_quant)
3. [Adult v. Hatchling contrast](#3AvH)


### 1. Trim and filter reads in prep for mapping <a name="1reads_prep"><\a>

	rawReadsDir=/home/scharer_group/raw_reads/sequencing_2020/RNA_Illumina

	maccli_A_reads="20201013075730868-60718056 20201013075928557-60718057 20201013080128925-60718059 20201013080334400-60718060"
	maccli_H_reads="20201013080537244-60718061 20201013080738726-60718062 20201013080941370-60718063 20201013081145613-60718064"
	maccli_reads=$(echo -n "${maccli_A_reads} " && echo "${maccli_H_reads}

Copy files over and unzip

	mkdir ${HOME}/data/mac_genomes/genome_Maccli/canu_assembly/final_annotation/AvH_contrast/raw_reads
	cd ${HOME}/data/mac_genomes/genome_Maccli/canu_assembly/final_annotation/AvH_contrast/raw_reads
	
	for i in $(echo ${maccli_reads}); do \
	       echo "gunzip -c ${rawReadsDir}/${i}/*R1_001_MM_1.fastq.gz > ./$(basename ${rawReadsDir}/${i}/*R1_001_MM_1.fastq.gz | sed 's;.gz;;g')" && \
	       echo "gunzip -c ${rawReadsDir}/${i}/*R2_001_MM_1.fastq.gz > ./$(basename ${rawReadsDir}/${i}/*R2_001_MM_1.fastq.gz | sed 's;.gz;;g')"; \
	done > copy_unzip.cmds
	
	parallel -a copy_unzip.cmds -j 20


#### 1.1 Trimming reads with cutadapt

	rawReadsDir=${HOME}/data/mac_genomes/genome_Maccli/canu_assembly/final_annotation/AvH_contrast/raw_reads
	trimDir=${HOME}/data/mac_genomes/genome_Maccli/canu_assembly/final_annotation/AvH_contrast/trimmed_reads
	mkdir ${trimDir}
	cd ${trimDir}
	
	# Trim the raw reads
	for read in $(ls -1 ${rawReadsDir}/*_R1_001_MM_1.fastq | sed 's;_R1_001_MM_1.fastq;;g'); do \
	    echo "" && \
	    newRead=$(basename ${read}) && \
	    echo "Running: cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  --minimum-length 50 --nextseq-trim=25 --cores 20 -o ${trimDir}/${newRead}_R1.adTr.fastq -p ${trimDir}/${newRead}_R2.adTr.fastq ${read}_R1_001_MM_1.fastq ${read}_R2_001_MM_1.fastq" &&\
	    echo "" && \
	    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  --minimum-length 50 --nextseq-trim=25 --cores 20 -o ${trimDir}/${newRead}_R1.adTr.fastq -p ${trimDir}/${newRead}_R2.adTr.fastq ${read}_R1_001_MM_1.fastq ${read}_R2_001_MM_1.fastq;
	done


### 2. Map reads and quantify <a name="2map_quant"><\a>

	trimDir=${HOME}/data/mac_genomes/genome_Maccli/canu_assembly/final_annotation/AvH_contrast/trimmed_reads
	reference_dir=${HOME}/data/mac_genomes/genome_Maccli/canu_assembly/final_annotation
	mapping_dir=${HOME}/data/mac_genomes/genome_Maccli/canu_assembly/final_annotation/AvH_contrast/salmon_out
	cd ${mapping_dir}

	conda activate salmon

Salmon requires a `transcript-to-gene` map file.

	cd ${reference_dir}
	grep "mRNA" Maccli_v2_rnm.gff3 | awk 'BEGIN{OFS=FS="\t"}{print $9}' | sed 's:;:\t:g' | sed 's;[a-zA-Z]*=;;g' > ./Maccli_v2_rnm_trans2gn_slm.txt

Salmon also requires an index of the reference sequences.

	cd ${mapping_dir}
	salmon index -t ${reference_dir}/cds/Maccli_v2.cds.fa -i ./Maccli_v2.cds

Then we simply run `salmon quant` with the input reads to perform pseudo-mapping and quantification at the same time.

	samples="Maccli_A_1 Maccli_A_2 Maccli_A_3 Maccli_A_4 Maccli_H_1 Maccli_H_2 Maccli_H_3 Maccli_H_4"

Map the trimmed reads

	for i in ${samples}; do
	    read1=$(ls ${trimDir}/*${i}*_R1.adTr.fastq)
	    read2=$(ls ${trimDir}/*${i}*_R2.adTr.fastq)
	    echo "Running salmon for: ${i} - $(basename ${read1}) $(basename ${read2})"
	    salmon quant -l A -p 30 --seqBias --gcBias --posBias -g ${reference_dir}/Maccli_v2_rnm_trans2gn_slm.txt -i ./Maccli_v2.cds -1 ${read1} -2 ${read2} --output ${i}; 
	done

Exit the conda environment

	conda deactivate


### 3. Perform Adult v. Hatchling gene expression contrast <a name="3AvH"><\a>

	samples="Maccli_A_1 Maccli_A_2 Maccli_A_3 Maccli_A_4 Maccli_H_1 Maccli_H_2 Maccli_H_3 Maccli_H_4"
	
	for samp in ${samples}; do \
	    scp $alcedo:/home/axel/data/mac_genomes/genome_Maccli/canu_assembly/final_annotation/AvH_contrast/salmon_out/${samp}/quant.sf ./${samp}_quant.sf;
	done

See the R notebook `AvH_contrasts.Rmd`















	
