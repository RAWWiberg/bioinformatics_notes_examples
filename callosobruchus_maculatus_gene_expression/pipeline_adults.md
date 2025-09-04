### Table of contents
0. [Packages and software](#0packages)
1. [Read pre-processing](#preprocess)
	1. [Adapter trimming](#adapter_trimming)
	2. [rRNA contamination](#rrna)
2. [Expression quantification](#expression)
	1. [Simple expression](#simple_expr)
	2. [Allele-specific expression](#ase)
3. [Expression analyses](#expr_analyses)

X. [References](#references)



### 0. Packages and software <a name="0packages"></a>
Here are some details on software used.

sortmerna, v. 4.3.5
cutadapt, v. 4.0
fastqc, v. 0.11.9
cdhit, v. 4.8.1
seqkit, v. 0.15.0

cufflinks, v. 2.2.1
salmon, v. 1.6.0
star, v. 2.7.9a

	workDir=${HOME}/projects/callosobruchus/

The commands below are example commands used in the analysis. 
For clarity, full paths are not shown and placeholder filenames are used.
Raw reads for the analysis can be found in NCBI under the BioProject: XXXX

### 0. Get the raw reads


Collect statistics for the read files. 

	seqkit stats ./*.fq.gz -j 10 -b -a -T > ./fastq_stats.tab

See the batch script ```reads_stats_adults.sh```


### 1.Read pre-processing <a name="preprocess"></a>


#### 1.1 Adapter trimming <a name="adapter_trimming"></a>  

	ls -1 *fastq.gz > adults_samples.list	

Trim Illumina TruSeq adapters with ```cutadapt```.  

	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --cores 10 -o ${sample}_adTr_R1.fastq.gz -p ${sample}_adTr_R2.fastq.gz ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz

This is performed in the batch script ```adults_adapter_trimming.sh```.

	for sample in $(cat adults_samples.list); do echo "sbatch -o ${sample}_adults_adapter_trimming.out adults_adapter_trimming.sh ${sample}"; done > submit_adults_adapter_trimming.sh
	. ./sumbit_adults_adapter_trimming.sh

Collect statistics for the final read files. 

	seqkit stats ./*.fq.gz -j 10 -b -a -T > ./fastq_stats.tab

See the batch script ```reads_stats_adults_adtr.sh```


#### 1.2 rRNA contamination <a name="rrna"></a>

I manually collected and downloaded all coleoptera rRRNA sequences from NCBI and SILVA databases.
See the file `data/rrna/ncbi_rrna.fa`  

SILVA: n = 13,768 sequences
NCBI: n = 2,646 sequences (search term: ("Coleoptera"[Organism] OR Coleoptera[All Fields]) AND rRNA[All Fields] AND biomol_rrna[PROP])

Cluster the sequences with cdhit to 95% similarity

	cd-hit-est -c 0.95 -i cmac_rrna.fa -o cmac_rrna_cdhit95.fa

The software ```sortmerna``` maps reads to an rRNA database. 
I label all reads that do not map to this database as "good" mRNA reads and discard the rest.
The .log file also gives the percent of reads mapping to the rRNA database.
	
	sortmerna -out2 -paired_in -v -threads 20 \
	--workdir ${SNIC_TMP} \
	--ref cmac_rrna.fa \
	-reads ${R1} -reads ${R2} \
	-fastx \
	--no-best \
	--num_alignments 1 \
	-aligned ${sample}_rrna \
	-other ${sample}_mrna
	

See the batch script ```rrna_removal.sh```.

Then run for the remaining samples that were only sequenced once.

	for sample in $(cat adults_samples.list); do echo "sbatch -o ${sample}_rrna_removal.out rrna_removal.sh ${sample}"; done > submit_rrna_removal_adults.sh
	. ./submit_rrna_removal_adults.sh

Collect some statistics from the logs.

	grep -A 1 "Coverage by database:" *log | grep "ncbi_silva_coleoptera" > sortmerna_mapping_rates.txt

Collect statistics for the final read files. See the batch script ```reads_stats.sh```

	sbatch -o rrna_removed_stats.out reads_stats.sh rrna_removed

Run FastQC for the final read files. See the batch script ```reads_fastqc.sh```

	sbatch -o rrna_reads_fastqc.out reads_fastqc.sh rrna_removed

Then summarise FastQC and SortMeRNA output with MultiQC


The generated file ```rrna_removed_fastq_stats.tab ``` needs to be manually mofidied a bit to make sample names more easily parseable.

See the script ```rnotebooks/reads_assessment.Rmd``` for more details of the assessment of trimming and read cleaning steps.


### 2. Expression quantification <a name="expression"></a>
#### 2.1 Simple expression <a name="simple_expr"></a>

Quantify expression simply at each annotated transcript.

	annotation=p1.GS.annotation.v1.gtf
	assembly=pt_036_001_hifiasm_20201223.primary.fasta

First extract all transcripts from the annotation `.gtf` file to make a gene to transcript mapping file for salmon.

	grep -P "\ttranscript\t" p1.GS.annotation.v1.gtf | awk	'BEGIN{OFS=FS="\t"}{print $9}' | sed 's;\(cmac_.*\)\(.t[0-9]\);\1\2\t\1;' > p1.GS.annotation.v1_trans2gene.tab
	
	trans2gene=p1.GS.annotation.v1_trans2gene.tab


##### 2.1.1 Salmon expression quantification
##### 2.1.1.1 Extract transcript sequences from genome and create salmon index

See the batch script ```extract_transcripts.sh```	

	gffread ${annotation} -w p1.GS.annotation.v1.transcripts.fa -g ${assembly}

	salmon index -p 3 -i ./p1.GS.annotation.v1.transcripts -t ./p1.GS.annotation.v1.transcripts.fa


##### 2.1.1.2 Map trimmed and cleaned reads to transcripts

Index the transcripts

	salmon quant -1 ${R1} -2 ${R2} -i p1.GS.annotation.v1.transcripts --seqBias --gcBias --posBias -p 10 -g ${trans2gene}

See the batch script ```salmon_run.sh```

	sbatch -o ${sample}_adults_salmon_run.out adults_salmon_run.sh ${sample}

	for sample in $(cat adults_samples.list); do echo "sbatch -o ${sample}_adults_salmon_run.out adults_salmon_run.sh ${sample}"; done > submit_adults_salmon_run.sh
	. ./submit_adults_salmon_run.sh


Collect some statistics from the log files.

	grep "Mapping rate" *_salmon_quant/logs/salmon_quant.log | sed 's;\[.*\];;g' | sed 's;_salmon_quant.*: Mapping rate =;;g' > salmon_mapping_rates.txt

	tar -czvf adults_salmon_logs.tar.gz *adults_salmon_sun.out


##### 2.1.2 STAR alignment of reads to genome


See the batch scripts `adults_start_resem_run.sh`


### 3. Expression analyses <a name="expr_analyses"></a>

For the full analysis of expression data as computed by `salmon` above, see the R scripts in the folder `rnotebooks`
