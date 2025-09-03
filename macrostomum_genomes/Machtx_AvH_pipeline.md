1. [Trim and filter reads](#1reads_prep)
2. [Map/Quantify](#2map_quant)
3. [Adult v. Hatchling contrast](#3AvH)


### 1. Trim and filter reads in prep for mapping <a name="1reads_prep"><\a>

	rawReadsDir=/home/scharer_group/raw_reads/2015.12/renamed/PE

	machtx_A_reads="Machtx_A_1 Machtx_A_2 Machtx_A_3"
	machtx_H_reads="Machtx_H_1 Machtx_H_2 Machtx_H_3"
	machtx_reads=$(echo -n "$(echo ${machtx_A_reads} | sed 's;\(_[0-9]\);\1_L001;g') " && echo -n "$(echo ${machtx_A_reads} | sed 's;\(_[0-9]\);\1_L002;g') " && echo -n "$(echo ${machtx_H_reads} | sed 's;\(_[0-9]\);\1_L001;g') " && echo -n "$(echo ${machtx_H_reads} | sed 's;\(_[0-9]\);\1_L002;g')")

Copy files over and unzip

	mkdir ${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/final_annotation/AvH_contrast/raw_reads
	cd ${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/final_annotation/AvH_contrast/raw_reads
	
	for i in $(echo ${machtx_reads}); do \
	       echo "gunzip -c ${rawReadsDir}/${i}_R1.fq.gz > ./$(basename ${rawReadsDir}/${i}_R1.fq.gz | sed 's;.gz;;g')" && \
	       echo "gunzip -c ${rawReadsDir}/${i}_R2.fq.gz > ./$(basename ${rawReadsDir}/${i}_R2.fq.gz | sed 's;.gz;;g')"; \
	done > copy_unzip.cmds
	
	parallel -a copy_unzip.cmds -j 20

There are two files per sample. Combine all reads for a sample into one file

	for i in $(echo ${machtx_A_reads}); do \
	       echo "cat ${i}_L00*_R1.fq > ${i}_all_R1.fq" && \
	       echo "cat ${i}_L00*_R2.fq > ${i}_all_R2.fq"; \
	done > concatenate_reads.cmds
	
	for i in $(echo ${machtx_H_reads}); do \
	       echo "cat ${i}_L00*_R1.fq > ${i}_all_R1.fq" && \
	       echo "cat ${i}_L00*_R2.fq > ${i}_all_R2.fq"; \
	done >> concatenate_reads.cmds
	
	parallel -a concatenate_reads.cmds -j 20

Remove the raw-reads files

	for i in $(ls *L00*.fq | grep -v "_all_"); do \
	    rm ${i};
	done


#### 1.1 Trimming reads with cutadapt

	rawReadsDir=${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/final_annotation/AvH_contrast/raw_reads
	trimDir=${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/final_annotation/AvH_contrast/trimmed_reads
	mkdir ${trimDir}
	cd ${trimDir}
	
	# Trim the raw reads
	for read in $(ls -1 ${rawReadsDir}/*_R1.fq | sed 's;_R1.fq;;g'); do \
	    echo "" && \
	    newRead=$(basename ${read}) && \
	    echo "Running: cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  --minimum-length 50 --nextseq-trim=25 --cores 20 -o ${trimDir}/${newRead}_R1.adTr.fq -p ${trimDir}/${newRead}_R2.adTr.fq ${read}_R1.fq ${read}_R2.fq" &&\
	    echo "" && \
	    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  --minimum-length 50 --nextseq-trim=25 --cores 20 -o ${trimDir}/${newRead}_R1.adTr.fq -p ${trimDir}/${newRead}_R2.adTr.fq ${read}_R1.fq ${read}_R2.fq;
	done


### 2. Map reads and quantify <a name="2map_quant"><\a>

	trimDir=${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/final_annotation/AvH_contrast/trimmed_reads
	reference_dir=${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/final_annotation
	mapping_dir=${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/final_annotation/AvH_contrast/salmon_out
	cd ${mapping_dir}

	conda activate salmon

Salmon requires a `transcript-to-gene` map file.

	cd ${reference_dir}
	grep "mRNA" Machtx_v2_rnm.gff3 | awk 'BEGIN{OFS=FS="\t"}{print $9}' | sed 's:;:\t:g' | sed 's;[a-zA-Z]*=;;g' > ./Machtx_v2_rnm_trans2gn_slm.txt

Salmon also requires an index of the reference sequences.

	cd ${mapping_dir}
	salmon index -t ${reference_dir}/cds/Machtx_v2.cds.fa -i ./Machtx_v2.cds

Then we simply run `salmon quant` with the input reads to perform pseudo-mapping and quantification at the same time.

	samples="Machtx_A_1 Machtx_A_2 Machtx_A_3 Machtx_H_1 Machtx_H_2 Machtx_H_3"

Map the trimmed reads

	for i in ${samples}; do
	    read1=$(ls ${trimDir}/*${i}*_R1.adTr.fq)
	    read2=$(ls ${trimDir}/*${i}*_R2.adTr.fq)
	    echo "Running salmon for: ${i} - $(basename ${read1}) $(basename ${read2})"
	    salmon quant -l A -p 30 --seqBias --gcBias --posBias -g ${reference_dir}/Machtx_v2_rnm_trans2gn_slm.txt -i ./Machtx_v2.cds -1 ${read1} -2 ${read2} --output ${i}; 
	done

Exit the conda environment

	conda deactivate


### 3. Perform Adult v. Hatchling gene expression contrast <a name="3AvH"><\a>

	samples="Machtx_A_1 Machtx_A_2 Machtx_A_3 Machtx_H_1 Machtx_H_2 Machtx_H_3"
	
	for samp in ${samples}; do \
	    scp $alcedo:/home/axel/data/mac_genomes/genome_Machtx/canu_assembly/final_annotation/AvH_contrast/salmon_out/${samp}/quant.sf ./${samp}_quant.sf;
	done

See the R notebook `AvH_contrasts.Rmd`















	
