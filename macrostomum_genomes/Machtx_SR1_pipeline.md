# Assembly, and annotation for *Macrostomum hystrix*

### Table of contents
0. [Packages and software](#0packages)
1. [Genome assembly [part 1]](#genome_assembly1)
	1. [Primary assembly](#prime_assembly)
	2. [Purge haplotigs](#purge_hapl)
	3. [Arrow polishing](#arrow)
	4. [Pilon polishing](#pilon)
	5. [Final assembly](#final_assembly)
2. [Mitochondiral genome identification and assembly](#mt_genome)
3. [Genome assembly [part 2]](#genome_assembly2)
	1. [Repeatmasking](#repeatmasking)
	2. [Re-mapping reads to final genome assembly](#re_mapping)
4. [Annotation](#annotation)
	1. [SL sequence identification](#SL_sequence_id)
	2. [RNA mapping to final assembly](#rna_mapping)
		1. [Trimming: Adapter sequences](#adapter_trimming)
		2. [Trimming: SL sequence](#SL_trimming)
		3. [Mapping: All RNA-seq reads](#RNAseq_reads_mapping)
	3. [BRAKER](#braker)
	4. [IsoSeq](#isoseq)
	5. [PASA](#pasa)
	6. [Final annotation](#final_annotation)
		1. [Annotation cleaning and collecting stats](#final_annotation_stats)
		2. [Mapping: SL RNA-seq reads](#SL_reads_mapping)
5. [References](#references)

This document details example code used for assembly and annotation of the *M. cliftonense* genome.
For clarity, I have ommitted, in most cases, full paths to input and output files.

### 0. Install packages, make environments, obtain docker containers <a name="0packages"></a>

Using conda v. 4.8.5
**BUSCO**

	conda create -n busco
	conda activate busco
	conda install busco
	conda deactivate
	conda list --export -n busco > busco_env.txt

**PB-assembly**

	conda create -n pb-assembly
	conda activate pb-assembly
	conda install pb-assembly
	conda deactivate
	conda list --export -n pb-assembly  > pb-assembly_env.txt

**proteinortho**

	conda create -n proteinortho
	conda activate proteinortho
	conda install proteinortho
	conda deactivate
	conda list --export -n proteinortho  > proteinortho_env.txt


To re-create these conda environments use the files given in the `conda_env_files` folder. e.g. 

	conda create -n pb-assembly --file pb-assembly_env.txt

For some tools I used available docker containers:

**PASA**

	docker pull pasapipeline/pasapipeline:latest

**MITObim**

This docker image contains MITObim 1.6 and 1.9.1, make sure you know which one you are using, they are labelled.

	docker pull chrishah/mitobim

My own helper scripts that are used throughout can be found in the folder `scripts`
Where these were used I prefix with a variable as a placeholder for the path to these scripts: `${scriptsDir}`

### 1. Genome assembly [part 1] <a name="genome_assembly1"></a>
#### 1.1. Primary assembly <a name="prime_assembly"></a>
Here I run canu on the raw PacBio reads to make an initial assembly.

First make .fastq files from the PacBio .bam files.

	PBreadsDir="${HOME}/data/mac_genomes/genome_Machtx/raw_reads/gDNA/PacBio"
	samtools fastq -@ 10 -0 ${PBreadsDir}/r54273_20200623_142718/1_A01/m54273_200623_144107.subreads.fastq ${PBreadsDir}/r54273_20200623_142718/1_A01/m54273_200623_144107.subreads.bam
	samtools fastq -@ 10 -0 ${PBreadsDir}/r54273_20200625_142010/1_A03/m54273_200625_143339.subreads.fastq ${PBreadsDir}/r54273_20200625_142010/1_A03/m54273_200625_143339.subreads.bam

Now run canu, including the read-correction steps.

	canu -d ~/data/genome_Machtx/canu_assembly -p Machtx_SR1 genomeSize=230m maxThreads=50 stopOnLowCoverage=10 \
	-pacbio-raw \
	${PBreadsDir}/r54273_20200623_142718/1_A01/m54273_200623_144107.subreads.fastq \
	${PBreadsDir}/r54273_20200625_142010/1_A03/m54273_200625_143339.subreads.fastq


Collect N50, assembly size, etc.

	assembly-stats *.fasta

Run BUSCO

	conda activate busco

	busco -i Machtx_SR1.contigs.fasta  --cpu 10 --out Machtx_SR1.contigs.busco -m genome -l metazoa_odb10

	conda deactivate


#### 1.2 Purge haplotigs <a name="purge_hapl"></a>
Because the initial assemblies were much larger than expect (~2x), a likely scenario is that we have assembled both haplotypes as separate contigs. The following steps identify these redundant contigs and filters out the shorter of the two haplotypes.

Map original PacBio reads	

	PBreadsDir="${HOME}/data/mac_genomes/genome_Machtx/raw_reads/gDNA/PacBio"
	minimap2 -t 10 --secondary=no -a -x map-pb Machtx_SR1.contigs.fasta ${PBreadsDir}/r54273_20200623_142718/1_A01/m54273_200623_144107.subreads.fastq ${PBreadsDir}/r54273_20200625_142010/1_A03/m54273_200625_143339.subreads.fastq | samtools sort -m 1G -@ 5 -o Machtx_SR1.contigs_readsmapped_nomult.bam -T tmp.ali
	samtools sort Machtx_SR1.contigs_readsmapped_nomult.bam > Machtx_SR1.contigs_readsmapped_nomult_srt.bam

Get the coverage histogram

	purge_haplotigs hist -g ../Machtx_SR1.contigs.fasta -b ../Machtx_SR1.contigs_readsmapped_nomult.bam -t 10
	purge_haplotigs cov -i Machtx_SR1.contigs_readsmapped_nomult.bam.gencov -l 2 -m 22 -h 190
	purge_haplotigs purge -g ../Machtx_SR1.contigs.fasta -c coverage_stats.csv -b ../Machtx_SR1.contigs_readsmapped_nomult.bam -t 10 -a 60

N50, assembly size, etc.

	assembly-stats *.fasta 

Run BUSCO for primary purged assembly

	conda activate busco
	busco -i curated_l2_m22_h190_a60.fasta --cpu 10 --out Machtx_SR1.purge_haplotigs_curated_l2_m22_h190_a60.busco -m genome -l metazoa_odb10
	conda deactivate


**[BLOBTOOLS?]**

#### 1.3 Arrow polishing <a name="arrow"></a>
Here I conduct several rounds of "polishing" with the raw PacBio reads.

Start polishing with `arrow`

	conda activate pb-assembly

	PBreadsDir="${HOME}/data/mac_genomes/genome_Machtx/raw_reads/gDNA/PacBio"
	ls -1 ${PBreadsDir}/*/*/*subreads.bam > bams.fofn
	pbmm2 align --sort -J 8 -j 20 ./curated_l2_m22_h190_a60.fasta bams.fofn ./Machtx_SR1_curated_l2_m22_h190_a60_pbmm2readsmapped_srt.bam
	pbindex Machtx_SR1_curated_l2_m22_h190_a60_pbmm2readsmapped_srt.bam


Polish w/ Arrow: (1st iteration)

	gcpp --algorithm arrow -j 10 -o ./Machtx_SR1_curated_l2_m22_h190_a60_arrow1.fasta --reference ./curated_l2_m22_h190_a60.fasta ./Machtx_SR1_curated_l2_m22_h190_a60_pbmm2readsmapped_srt.bam

	pbmm2 align --sort -J 8 -j 30 ./Machtx_SR1_curated_l2_m22_h190_a60_arrow1.fasta bams.fofn ./Machtx_SR1_curated_l2_m22_h190_a60_arrow1_pbmm2readsmapped_srt.bam

	pbindex ./Machtx_SR1_curated_l2_m22_h190_a60_arrow1_pbmm2readsmapped_srt.bam

Polish w/ Arrow (x iteration)

	iter=2 # Current iteration
	x=$((${iter}-1)) # Which genome to use
	ls -1 ${reads_Dir}/*/*/*subreads.bam > bams.fofn

	gcpp --algorithm arrow -j 40 -o ./Machtx_SR1_curated_l2_m22_h190_a60_arrow${iter}.fasta --reference ./Machtx_SR1_curated_l2_m22_h190_a60_arrow${x}.fasta ./Machtx_SR1_curated_l2_m22_h190_a60_arrow${x}_pbmm2readsmapped_srt.bam

	pbmm2 align --sort -J 8 -j 40 ./Machtx_SR1_curated_l2_m22_h190_a60_arrow${iter}.fasta bams.fofn ./Machtx_SR1_curated_l2_m22_h190_a60_arrow${iter}_pbmm2readsmapped_srt.bam

	pbindex ./Machtx_SR1_curated_l2_m22_h190_a60_arrow${iter}_pbmm2readsmapped_srt.bam

	conda deactivate

Run BUSCO

	conda activate busco
	busco -f -i Machtx_SR1_curated_l2_m22_h190_a60_arrow${iter}.fasta --cpu 20 --out Machtx_SR1_curated_l2_m22_h190_a60_arrow${iter}.busco -m genome -l metazoa_odb10
	conda deactivate


Evaluate all assemblies with quast.

	quast.py -t 30 -o ../quast_assessment --eukaryote --gene-finding \
	--bam ./Machtx_SR1_curated_l2_m22_h190_a60_pbmm2readsmapped_srt.bam,./Machtx_SR1_curated_l2_m22_h190_a60_arrow1_pbmm2readsmapped_srt.bam,./Machtx_SR1_curated_l2_m22_h190_a60_arrow2_pbmm2readsmapped_srt.bam,./Machtx_SR1_curated_l2_m22_h190_a60_arrow3_pbmm2readsmapped_srt.bam,./Machtx_SR1_curated_l2_m22_h190_a60_arrow4_pbmm2readsmapped_srt.bam \
	./curated_l2_m22_h190_a60.fasta ./Machtx_SR1_curated_l2_m22_h190_a60_arrow1.fasta ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2.fasta ./Machtx_SR1_curated_l2_m22_h190_a60_arrow3.fasta ./Machtx_SR1_curated_l2_m22_h190_a60_arrow4.fasta


#### 1.4 Pilon polishing <a name="pilon"></a>
Here I conduct several iterations of "polishing" with Illumina short read data.

	conda activate pb-assembly

Trim adapters from illumina reads

	IlluminareadsDir="/home/scharer_group/raw_reads/sequencing_2020/gDNA_Illumina/20201013083400233-60718076"
	
	adapters=${HOME}/packages/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa

	trimmomatic PE -threads 10 -phred33 -trimlog ./BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1.log -summary ./BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1.summary ${IlluminareadsDir}/BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R1_001_MM_1.fastq.gz ${IlluminareadsDir}/BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R2_001_MM_1.fastq.gz ./BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R1_001_MM_1_adTr.fq.gz ./BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R1_001_MM_1_unpaired_adTr.fq.gz ./BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R2_001_MM_1_adTr.fq.gz ./BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R2_001_MM_1_unpaired_adTr.fq.gz ILLUMINACLIP:${adapters}:2:30:10

Polish with pilon
use the arrow2 assembly

	cp Machtx_SR1_curated_l2_m22_h190_a60_arrow2.fasta Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon0.fasta

	bwa mem -t 20 ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2.fasta $trimReadsDir/BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R1_001_MM_1_adTr.fq.gz $trimReadsDir/BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R2_001_MM_1_adTr.fq.gz | samtools view -Sb | samtools sort -@ 20 > ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2_D126_127-bwaIlluminaReadsMapped.srt.bam

	samtools index Machtx_SR1_curated_l2_m22_h190_a60_arrow2_D126_127-bwaIlluminaReadsMapped.srt.bam

Map reads with bwa mem

Polish with pilon (x iteration)

	iter=1 # Current iteration
	x=$((${iter}-1)) # Which genome to use

	bwa index ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon${x}.fasta

	bwa mem -t 30 ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon${x}.fasta BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R1_001_MM_1_adTr.fq.gz BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R2_001_MM_1_adTr.fq.gz | samtools view -Sb | samtools sort -@ 20 > ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon${x}_D126_127-bwaIlluminaReadsMapped.srt.bam

	pilon --genome ./Machtx_SR1_curated_l2_m22_h190_a60_arrow${x}.fasta --bam ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2_D126_127-bwaIlluminaReadsMapped.srt.bam --output ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon${iter} --minmq 10 --outdir ./ --tracks --changes --threads 50

	. ${scriptsDir}/count_changes.sh Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon${iter}.changes | awk 'BEGIN{FS=":"}{print $2}'	

	conda deactivate
	

Run BUSCO for the assembly

	conda activate busco
	busco -f -i ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon${iter}.fasta --cpu 50 --out Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon${iter}.busco -m genome -l metazoa_odb10
	conda deactivate

Run QUAST for all assemblies

	bams=$(echo *.bam | sed 's; ;,;g')
	bams=$(echo $bams | sed 's;Machtx;./Machtx;g')
	fastas=$(echo *pilon*fasta)
	fastas=$(echo $fastas | sed 's;Machtx;./Machtx;g')
	quast.py -t 20 -o ../quast_assessment_pilon --eukaryote --gene-finding --bam ${bams} ${fastas}


#### 1.5 Final assembly <a name="final_assembly"></a>

I get rid of a bunch of superfluous labels in the fasta headers

	sed -i 's;|.*;;g' Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5.fasta

Remove lower-case letters left over from arrow polishing runs.

	seqkit seq --upper-case ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5.fasta > ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_cln.fasta
	mv ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_cln.fasta ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5.fasta
	

**mt Genome**
Here I identify any contigs that look like they might be the mitochondrial genome. 
I decided that I would assemble and annotate these separately (see **6. mitochondrial genome**). 
So here I first remove contigs from the final genome assembly that get hits for available COI sequences. 
I then carry on with annotation and comparative analyses while also assembling the mt genome in parallel

Blast Machtx COI sequences against the final genome assembly

	makeblastdb -in Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5.fasta -dbtype nucl -title Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5 -out Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5.fasta -parse_seqids -hash_index
	blastn -db ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5.fasta -query ${mitoDir}/Machtx_COI.fasta -outfmt 6

I get no hits, whaaaaat!?!?!

This is probably because the gDNA fragments were size-selected to be as long or longer than the mtGenome.
Additionally, coverage was lower for the M. hystrix sequencing runs than for the M. cliftonense runs
So very few reads (if any) will be produced from this.

There is no contig in the Machtx assembly that is likely to be the mt genome. 
I therefore label the current assembly as `nomito` and carry on with annotation and comparative analyses of the genome assembly.
In [Section 2](#mt_genome) I *de novo* assemble the mt genome from Illumina short reads.

	cp ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5.fasta ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito.fasta

Lets also re-name the file.

	cp ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito.fasta ./Machtx_SR1_v2.fasta

Compute basic stats for each contig

	fastaStats.py -f ./Machtx_SR1_v2.fasta -t p -s n > ./Machtx_SR1_v2.fasta.stats



### 2. Mitochondrial genome identification and assembly <a name="mt_genome"></a>

**MITObim 2**
I use the Maclig mt_genome as an initial bait, I use the de novo assembly option of MITObim2

Copy over the read files

	trimReadsDir=${HOME}/data/genome_Machtx/raw_reads/gDNA/Illumina
	cp $trimReadsDir/BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_*_001_MM_1.fastq.gz ./

Convert to interleaved format and clean to make them compatible with MITObim

	rename.sh in=BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R1_001_MM_1_adTr.fq.gz in2=BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R2_001_MM_1_adTr.fq.gz out=Machtx_SR1_illumina_gdna_adTr_interleaved.fastq prefix=read_
	
	rm ./BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_*_001_MM_1_adTr.fq.gz
	sed -i 's; ;\/;g' Machtx_SR1_illumina_gdna_adTr_interleaved.fastq
	sed -i 's;\(\/[[:digit:]]\):;\1;g' Machtx_SR1_illumina_gdna_adTr_interleaved.fastq

Create MIRA manifest file. 
Should look like:

	___________________________________________________________
	project = initial-mapping-testpool-to-Maclig-mtgenome

	job=genome,mapping,accurate

	parameters = -NW:mrnl=0 -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no

	readgroup
	is_reference
	data = Maclig_mtgenome.fna
	strain = Maclig-mtgenome
	
	readgroup = reads
	autopairing
	data = Machtx_SR1_illumina_gdna_adTr_interleaved.fastq
	technology = solexa
	strain = Machtx_SR1  
	___________________________________________________________

The next steps are all done within the mitobim docker container. 
Spin up docker container.
e.g.

	docker run --rm -it -v ${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/mitochondrial_shenanigans/mitobim2:/home/data chrishah/mitobim /bin/bash

Run MIRA

	cd ./data
	mira manifest.conf

Run MITObim

	date=$(date | sed 's; ;_;g' | sed 's;:;;g')
	/home/src/scripts/MITObim.pl -start 1 -end 50 --pair -sample Machtx_SR1 -ref Maclig_mtgenome -readpool ./Machtx_SR1_illumina_gdna_adTr_interleaved.fastq -maf ./initial-mapping-testpool-to-Maclig-mtgenome_assembly/initial-mapping-testpool-to-Maclig-mtgenome_d_results/initial-mapping-testpool-to-Maclig-mtgenome_out.maf &> log_$date

Make sure I'm allowed to use the files that have been created 

	chown --recursive 1010:1004 ./*
	
	exit

Copy final assembly .fasta file and clean it up.

Remove IUPAC characters
Rename the contig so its clear that its the mtGenome
Replace any 'X' characters with 'N'

	seqkit seq -u ./iteration37/Machtx_SR1-Maclig_mtgenome-it37_noIUPAC.fasta | sed 's;>scaf0.*;>scaf0_mito;g' > ./iteration37_Machtx_SR1-Machtx_SR1-Maclig_mtgenome-it37_noIUPAC.fasta
	seqkit seq -w 60 -u ./iteration37_Machtx_SR1-Machtx_SR1-Maclig_mtgenome-it37_noIUPAC.fasta > Machtx_SR1_mitobim_mtgenome_noIUPAC.fasta
	sed '/^[^>]/s/X/N/g' Machtx_SR1_mitobim_mtgenome.fasta > Machtx_SR1_mitobim_mtgenome_noX.fasta

Get the reads that were used for the final iteration

	cp ./iteration37/Machtx_SR1-readpool-it37.fastq ./iteration37_Machtx_SR1-readpool-it37.fastq

Sanity check for the presence of COI genes

	makeblastdb -in ./Machtx_SR1_mitobim_mtgenome_noIUPAC.fasta -dbtype nucl -title Machtx_SR1_mitobim_mtgenome_noIUPAC -out Machtx_SR1_mitobim_mtgenome_noIUPAC.fasta -parse_seqids -hash_index
	blastn -db ./Machtx_SR1_mitobim_mtgenome_noIUPAC.fasta -query ../Machtx_COI.fasta -outfmt 6
	blat ./Machtx_SR1_mitobim_mtgenome_noIUPAC.fasta ../Machtx_COI.fasta Machtx_COI_Machtx_SR1_mitobim_mtgenome.psl

This produces hits for all sequences.

Map the caught reads from iteration7 to the unpadded assembly
	
	conda activate pb-assembly

	bwa index ./Machtx_SR1_mitobim_mtgenome_noIUPAC.fasta

	bwa mem -p -t 50 ./Machtx_SR1_mitobim_mtgenome_noIUPAC.fasta ./iteration37_Machtx_SR1-readpool-it37.fastq | samtools view -Sb | samtools sort -@ 50 > ./Machtx_SR1_mitobim_mtgenome_noIUPAC_iteration37IlluminareadsMapped.srt.bam
	samtools index Machtx_SR1_mitobim_mtgenome_noIUPAC_iteration37IlluminareadsMapped.srt.bam
	samtools flagstat Machtx_SR1_mitobim_mtgenome_noIUPAC_iteration37IlluminareadsMapped.srt.bam
	
	bamCoverage -p 50 --bam Machtx_SR1_mitobim_mtgenome_noIUPAC_iteration37IlluminareadsMapped.srt.bam --outFileFormat bedgraph --outFileName Machtx_SR1_mitobim_mtgenome_noIUPAC_iteration37IlluminareadsMapped.srt.cov.bedgraph

The Mitochondrial annotation server MITOS also independently identifies cox-1 in the same regions and warns that they are split.

Add the mitobim mitochondrial genome to the assembly

	cat ./Machtx_SR1_v2.fasta ${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/mitochondrial_shenanigans/mitobim2/Machtx_SR1_mitobim_mtgenome_noIUPAC.fasta > ./Machtx_SR1_v2_wMito.fasta 



### 3. Genome assembly [part 2] <a name="genome_assembly2"></a>
#### 3.1 Repeatmasking <a name="repeatmasking"></a>

	BuildDatabase -name machtx Machtx_SR1_v2.fasta
	RepeatModeler -database machtx -pa 30 -LTRStruct

	RepeatClassifier -consensi machtx-families.fa

Softmasking

	RepeatMasker -pa 20 -xsmall -engine rmblast -lib machtx-families.fa.classified Machtx_SR1_v2.fasta

Hardmasking

	RepeatMasker -pa 20 -gff -engine rmblast -lib machtx-families.fa.classified Machtx_SR1_v2.fasta


#### 3.2 Re-mapping reads to final genome assembly <a name="re_mapping"></a>
Use the `nomito` or `nomito_rptsftmsk` genomes for this.

Re-map the Illumina reads

	trimReadsDir=${HOME}/data/mac_genomes/genome_Machtx/raw_reads/gDNA/Illumina
	bwa index ./Machtx_SR1_v2.fasta

	bwa mem -t 50 ./Machtx_SR1_v2.fasta $trimReadsDir/BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R1_001_MM_1_adTr.fq.gz $trimReadsDir/BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_R2_001_MM_1_adTr.fq.gz | samtools view -Sb | samtools sort -@ 50 > ./Machtx_SR1_v2_D126_127-bwaIlluminaReadsMapped.srt.bam

	samtools index Machtx_SR1_v2_D126_127-bwaIlluminaReadsMapped.srt.bam

Re-map the PacBio reads to this final assembly

	PBreadsDir="${HOME}/data/mac_genomes/genome_Machtx/raw_reads/gDNA/PacBio"
	ls -1 ${PBreadsDir}/*/*/*subreads.bam > bams.fofn

	pbmm2 align --sort -J 20 -j 50 ./Machtx_SR1_v2.fasta bams.fofn ./Machtx_SR1_v2_pbmm2readsmapped_srt.bam

	pbindex ./Machtx_SR1_v2_pbmm2readsmapped_srt.bam

#### 3.3 QC the final and initial assemblies using merqury <a name="merqury"></a>

First get k-mers from the reads

	sh $MERQURY/best_k.sh 217000000

Here I use the reads from the strain that was sequenced (also used for polishing)
I use the raw canu assembly as well as the final assembly (+ mtGenome)

	meryl k=19 count ${trimReadsDir}/BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002_*_001_MM_1_adTr.fq.gz output ${trimReadsDir}/BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002.meryl
	
	./merqury.sh ${trimReadsDir}/BSSE_QGF_152428_HTWN3DRXX_2_D126_127_AGTTCAGG_CCAACAGA_S1_L002.meryl ./Machtx_SR1.contigs.fasta ././Machtx_SR1_v2_wMito.fasta merq_D126_127_v_Machtx_final_wmito_and_raw.out



### 4. Annotation <a name="annotation"></a>
#### 4.1 SL sequence identification <a name="SL_sequence_id"></a>
Here I use a simple approach to identify spliced leader (SL) sequences.
This approach takes very common k-mers from previously published *de novo* transcriptome assemblies. 
These k-mers are then assembled with velvet. 
The resulting sequences look very similar to previously identified *M. lignano* SL sequences.

Collect the available transcriptome assemblies for Machtx from: 
[https://doi.org/10.5281/zenodo.4543289](https://doi.org/10.5281/zenodo.4543289)
I use the `noPrimer` data

	python3 ${scriptsDir}/100bpFasta.py -f ~/data/assemblies/assemblies_2018/Machtx_20180703_CroCo.Trinity.fasta -n 100 > Machtx_20180703_CroCo.Trinity.100bp.fasta
	jellyfish count -s 100M -m 19 -t 10 ./Machtx_20180703_CroCo.Trinity.100bp.fasta
	jellyfish dump -c mer_counts.jf > mer_counts.fas
	awk 'BEGIN{FS=OFS=" "}{if($2 > 1000) print ">"$2"\n"$1}' mer_counts.fas > common_mers.fasta
	velveth ./ 10 common_mers.fasta -short
	velvetg ./ -cov_cutoff 1 -min_contig_lgth 10 -coverage_mask 0 -exp_cov 3
	cat contigs.fa > Machtx_SL.fasta


#### 4.2 RNA mapping to final assembly <a name="rna_mapping"></a>
Trimming and mapping of RNA-seq reads to the genome assembly.

	cd ${HOME}/data/genome_Machtx/canu_assembly/rna_mapping_v2
	mkdir final_assembly
	mkdir final_assembly_softmasked

##### 4.2.1 Trimming: adapters <a name="adapter_trimming"></a>

This was done for many individual read files. Here I only give example code for a single instance.

	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	--cores 20 --info-file ${outDir}/${newRead}_adTr_cutadapt.info \
	-o ${READ}_adTr_R1_001_MM_1.fastq.gz -p ${READ}_adTr_R2_001_MM_1.fastq.gz \
	${READ}_R1_001_MM_1.fastq.gz ${READ}_R2_001_MM_1.fastq.gz

##### 4.2.2 Trimming: SL sequence <a name="SL_trimming"></a>
In stranded TruSeq libraries, the SL sequence should only be present on R2 reads. 
Nevertheless, to trim, we run cutadapt as for Illumina adapters. 
This ensures that reads which end up too short are thrown out as pairs.

This was done for many individual read files. Here I only give example code for a single instance.

	cutadapt -g CGTAAAGACGGTCTCTTACTGCGAAGACTCAATTTATTGCATG -O 6 --minimum-length 50 \
	--cores 20 -o ${READ}_adTrSLTr_R1_001_MM_1.fastq.gz -p ${READ}_adTrSLTr_R2_001_MM_1.fastq.gz \
	${READ}_adTr_R1_001_MM_1.fastq.gz ${READ}_adTr_R2_001_MM_1.fastq.gz
		

However, to identify which reads contain the SL sequence, we can run cutadapt for only the R2 reads in single-end mode and create an “info file”. Parsing the info file allows us to then extract which reads contain the SL sequence and to map only these to the genome.

	cutadapt -g GCCGTAAAGACGGTCTCTTACTGCGAAGACTCAATTTATTGCATG -O 6 --minimum-length 50 --cores 20 \
	--info-file ${READ}_adTrSLTr_cutadapt-temp.info \
	-o ${READ}_adTrSLTr_R2_001_MM_1-temp.fastq.gz ${READ}_adTr_R2_001_MM_1.fastq.gz
	awk 'BEGIN{FS=OFS="\t"}{if($2 != -1) print $0}' ${READ}_adTrSLTr_cutadapt-temp.info > ${READ}_adTrSLTr_cutadapt_SLreads.info
	awk 'BEGIN{FS=OFS="\t"}{print $1}' ${READ}_adTrSLTr_cutadapt_SLreads.info > ${READ}_adTrSLTr_cutadapt_SLreads.list && \

	rm ${READ}_adTrSLTr_cutadapt-temp.info

	seqkit grep -n -f ${READ}_adTrSLTr_cutadapt_SLreads.list ${READ}_adTrSLTr_R2_001_MM_1-temp.fastq.gz > ${READ}_adTrSLTr_R2_001_MM_1-SLreads.fastq.gz

	rm ${READ}_adTrSLTr_R2_001_MM_1-temp.fastq.gz

##### 4.2.3 Mapping: All RNA-seq reads <a name="RNAseq_reads_mapping"></a>
I map RNA-seq reads to the softmasked (repeats) final assembly.

Make an index for the genome

	STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ${genomeDir}/Machtx_SR1_v2.rptsftmsk.fasta --runThreadN 10 -sjdbOverhang ReadLength-1 --genomeSAindexNbases 12

Map the trimmed reads

	rnaDir=${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/rna_mapping_v2/SL_trimmed_reads/cutadapt

	readFiles_R1=$(echo $(ls ${rnaDir}/*R1*fq.gz) | grep -v "SLreads" | sed 's; ;,;g')
	readFiles_R2=$(echo ${readFiles_R1} | sed 's;_R1.fq.gz;_R2.fq.gz;g')

	ulimit -n 2048 # Need to set a higher limit to the maximum number of open files (originally at 1024)

	STAR --runThreadN 30 --genomeDir ./ --readFilesIn ${readFiles_R1} ${readFiles_R2} \
	--readFilesCommand zcat --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0.5 --outFilterMatchNmin 0 --alignEndsType EndToEnd \
	--outFileNamePrefix ./Machtx_SR1_v2_rptsftmsk_RNA_filt --outSAMstrandField intronMotif --outSAMtype BAM Unsorted

	samtools sort -@ 30 -o Machtx_SR1_v2_rptsftmsk_RNA_filtAligned.out_srt.bam Machtx_SR1_v2_rptsftmsk_RNA_filtAligned.out.bam

	samtools index Machtx_SR1_v2_rptsftmsk_RNA_filtAligned.out_srt.bam

	samtools flagstat Machtx_SR1_v2_rptsftmsk_RNA_filtAligned.out_srt.bam

	bamCoverage -p 30 --bam Machtx_SR1_v2_rptsftmsk_RNA_filtAligned.out_srt.bam --outFileFormat bedgraph --outFileName Machtx_SR1_v2_rptsftmsk_RNA_filtAligned.out_srt.cov.bedgraph

	samtools stats -@ 30 Machtx_SR1_v2_rptsftmsk_RNA_filtAligned.out_srt.bam > Machtx_SR1_v2_rptsftmsk_RNA_filtAligned.out_srt.stats


#### 4.3 BRAKER <a name="braker"><\a>

	braker.pl --verbosity=4 --cores 20 --gff3 --UTR=on --softmasking --species=Machtx --genome=Machtx_v2_rptsftmsk.rptsftmsk.fasta --bam=Machtx_v2_rptsftmsk_filtAligned.out_srt.bam

Run BUSCO for the BRAKER annotation

	conda activate busco
	busco -f -i augustus.hints_utr.codingseq -c 30 -m tran --out augustus.hints_utr.codingseq-busco -l ~/miniconda3/envs/busco/busco_downloads/lineages/metazoa_odb10
	conda deactivate


#### 4.4 IsoSeq <a name="isoseq"><\a>
These steps run the PacBio Isoseq3 pipeline as well as pos-processing steps for the IsoSeq data to produce full length expressed transcripts.

	conda activate pb-assembly

Data come from two “sources”. The same library was sequenced twice.

	sample1 # SEQUEL I, D-BSSE
	sample2 # SEQUEL II, Lausanne

**IsoSeq3 Pipeline**

	readsDir=${HOME}/data/genome_Machtx/canu_assembly/isoseq/raw_reads

Merge sample1 and sample2 bams to process them together.

	ls -1 ${readsDir}/sample*.ccs.bam > ${readsDir}/ccs_bamlist.fofn
	samtools merge -@ 10 -b ${readsDir}/ccs_bamlist.fofn ${readsDir}/sample1_2.ccs.bam
	ls -1 ${readsDir}/sample*.subreads.bam > ${readsDir}/subreads_bamlist.fofn
	samtools merge -@ 10 -b ${readsDir}/subreads_bamlist.fofn ${readsDir}/sample1_2.subreads.bam

Remove primers from ccs reads.

	lima ${readsDir}/sample1_2.ccs.bam ./primers_empirical.fasta ./sample1+2/sample1_2.noprimer_empirical.ccs.bam --isoseq

Refine reads by removing poly A tails and solving chimaeric reads (sample1+2)

	isoseq3 refine -j 20 ./sample1+2/sample1_2.noprimer_empirical.ccs.IsoSeq_Express_Primer_3p--NEB_5p.bam ./primers_empirical.fasta ./sample1+2/sample1_2.noprimer_empirical.ccs.flnc.bam --require-polya

Cluster similar isoforms

	isoseq3 cluster -j 20 ./sample1+2/sample1_2.noprimer_empirical.ccs.flnc.bam ./sample1+2/sample1_2.noprimer_empirical.ccs.clustered.bam --verbose --use-qvs

Polish isoforms

	isoseq3 polish -j 20 ./sample1+2/sample1_2.noprimer_empirical.ccs.clustered.bam ${readsDir}/sample1_2.subreads.bam ./sample1+2/sample1_2.noprimer_empirical.ccs.polished.bam

Map the polished transcripts to the genome assembly

	minimap2 -ax splice -t 10 -uf --secondary=no -C5 ./Machtx_SR1_v2.fasta ./sample1_2.noprimer_empirical.ccs.polished.hq.fastq > ./Machtx_SR1_v2_IsoSeq-sample1_2PolishedHQFLTrans.sam

Use this .sam file for cDNA_Cupcake

**cDNA_Cupcake ToFU Pipeline**

Collapse the transcripts by their mapping locations

	sort -k 3,3 -k 4,4n ./Machtx_SR1_v2_IsoSeq-sample1_2PolishedHQFLTrans.sam > ./Machtx_SR1_v2_IsoSeq-sample1_2PolishedHQFLTrans.srt.sam

	gunzip ./sample1_2.noprimer_empirical.ccs.polished.hq.fastq.gz

	collapse_isoforms_by_sam.py -c 0.98 -i 0.95 --input ./sample1_2.noprimer_empirical.ccs.polished.hq.fastq --fq -s ./Machtx_SR1_v2_IsoSeq-sample1_2PolishedHQFLTrans.srt.sam --dun-merge-5-shorter -o ./sample1_2.noprimer_empirical.ccs.polished.hq

	seqtk seq -A ./sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fq | sed 's;/;_;g' > ./sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa

	gffread sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.gff > sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.gff3

These are the final "non-redundant" transcripts/isoforms.


#### 4.5 PASA <a name="pasa"><\a>
These steps first aligns IsoSeq collapsed transcripts to genome with PASA and uses them as "known" expressed sequence tags (ESTs). These alignments are then updated using the BRAKER annotation.

Copy over the sqlite config files and modify the database name

	mkdir sqlite.confs
	cp /home/axel/packages/PASApipeline/sample_data/sqlite.confs/alignAssembly.config /home/axel/data/genome_Machtx/canu_assembly/pasa/sqlite.confs/
	cp /home/axel/packages/PASApipeline/sample_data/sqlite.confs/annotCompare.config /home/axel/data/genome_Machtx/canu_assembly/pasa/sqlite.confs/

In these files change database name and path: 

	"DATABASE=/tmp/sample_mydb_pasa.sqlite" > "DATABASE=/home/axel/data/genome_Machtx/canu_assembly/pasa/Machtx_SR1_v2_mydb_pasa.sqlite"

	cp /home/axel/data/genome_Machtx/canu_assembly/final_assembly/Machtx_SR1_v2.fasta ./

The PASA step usese the output from BRAKER and isoseq3 files and prepare them

**ISOSEQ3**

	cp /home/axel/data/genome_Machtx/canu_assembly/isoseq/isoseq3/sample1+2/sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa ./
	
	sed -i 's;|.*;;'g ./sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa
	grep ">" ./sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa | sed 's/>//g' > ./sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.names

**BRAKER**	

	cp /home/axel/data/genome_Machtx/canu_assembly/braker_annotation/augustus.hints_utr.gff3 ./
	

Start the `pasa` container

	docker run --rm -it -v ${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/pasa:${HOME}/data/genome_Machtx/canu_assembly/pasa pasapipeline/pasapipeline:latest bash

The next few lines are run within the docker container:

	cd ${HOME}/data/genome_Machtx/canu_assembly/pasa
	PASAHOME="/usr/local/src/PASApipeline"
	
Align and assemble the IsoSeq data using both gmap and blat, this also creates a transcripts database that we will update later on.

	${PASAHOME}/Launch_PASA_pipeline.pl --replace -c sqlite.confs/alignAssembly.config -C -R --ALIGNER gmap,blat --CPU 20 -g Machtx_SR1_v2.fasta -t sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa -f sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.names

	${PASAHOME}/scripts/build_comprehensive_transcriptome.dbi -c sqlite.confs/alignAssembly.config -t sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa --min_per_ID 95 --min_per_aligned 30

**PASA Update: Round 1**  
Add the BRAKER transcript file to the database and compare/update the annotation using the IsoSeq data
	
	${PASAHOME}/misc_utilities/pasa_gff3_validator.pl augustus.hints_utr.gff3
	
	${PASAHOME}/scripts/Load_Current_Gene_Annotations.dbi -c sqlite.confs/alignAssembly.config -g Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito.fasta -P augustus.hints_utr.gff3

Compare the annotations

	${PASAHOME}/Launch_PASA_pipeline.pl -c sqlite.confs/annotCompare.config -A -g Machtx_SR1_v2.fasta -t sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa

**PASA Update: Round 2**  
Add the PASA output from the previous round to the database and compare/update the annotation using the IsoSeq data
	
	${PASAHOME}/misc_utilities/pasa_gff3_validator.pl Machtx_SR1_v2_mydb_pasa.sqlite.gene_structures_post_PASA_updates.23.gff3
	
	${PASAHOME}/scripts/Load_Current_Gene_Annotations.dbi -c sqlite.confs/alignAssembly.config -g Machtx_SR1_v2.fasta -P Machtx_SR1_v2_mydb_pasa.sqlite.gene_structures_post_PASA_updates.23.gff3

Compare the annotations

	${PASAHOME}/Launch_PASA_pipeline.pl -c sqlite.confs/annotCompare.config -A -g Machtx_SR1_v2.fasta -t sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa

Make sure I'm allowed to use the files that have been created 

	chown --recursive 1010:1004 ./*
	
	exit


#### 4.6 Final annotation <a name="final_annotation"></a>
##### 4.6.1 Annotation cleaning and collecting stats <a name="final_annotation_stats"></a>

	cp ${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/pasa/Machtx_SR1_v2_mydb_pasa.sqlite.gene_structures_post_PASA_updates.18.gff3 ./

Remove superfluous entries and empty lines in the gff3 file

	grep -v "^#" Machtx_v2_mydb_pasa.sqlite.gene_structures_post_PASA_updates.18.gff3 | sed '/^$/d' > Machtx_SR1_v2_mydb_pasa.sqlite.gene_structures_post_PASA_updates.18_reduced.gff3
	
Make some changes to the final annotation file.  
Assures that "ID=" attributes field has a unique value for each CDS feature

	python3 ${scriptsDir}/gff_mod.py -gff Machtx_SR1_v2_mydb_pasa.sqlite.gene_structures_post_PASA_updates.18_reduced.gff3 > Machtx_v2.gff3

Remove all "Name=" attributes for gene and mRNA features

	sed -i 's/;Name=.*//g' Machtx_SR1_v2.gff3

Then re-name all the genes so that they follow the current *M. lignano* convention.

	python3 ${scriptsDir}/gff3_rename.py -gff Machtx_SR1_v2.gff3 -prefix Mhtx > Machtx_SR1_v2_rnm.gff3

Extract cds and aa sequences

	gffread ./Machtx_v2_rnm.gff3 -x ./Machtx_SR1_v2.cds.fa -g Machtx_SR1_v2.rptsftmsk.fasta
	gffread ./Machtx_v2_rnm.gff3 -y ./Machtx_SR1_v2.pep.fa -g Machtx_SR1_v2.rptsftmsk.fasta

Run BUSCO for the final annotated transcripts

	conda activate busco
	cd cds
	busco -i Machtx_SR1_v2.cds.fa --cpu 20 --out Machtx_SR1_v2.cds.busco -m trans -l metazoa_odb10
	conda deactivate

Get some stats for the annotations.

	python3 ${scriptsDir}/gff3_stats.py -gff Machtx_SR1_v2_rnm.gff3 -mode g > Machtx_SR1_v2_rnm.gff3.genes.stats
	python3 ${scriptsDir}/gff3_stats.py -gff Machtx_SR1_v2_rnm.gff3 -mode t > Machtx_SR1_v2_rnm.gff3.transcripts.stats

Identify the genes with "invalid" coding sequences (length not multiple of 3)

	awk 'BEGIN{FS=OFS="\t"}{if ($6 == 0) print $2}' Machtx_SR1_v2_rnm.gff3.transcripts.stats | sort -k1,1 | uniq > Machtx_SR1_v2_rnm_invalid.list

I also produce a filtered annotation file where I have remove genes with coding-sequences whose length is not a multiple of 3.

	python3 ${scriptsDir}/gff3_remove.py -gff Machtx_SR1_v2_rnm.gff3 -list Machtx_SR1_v2_rnm_invalid.list > Machtx_SR1_v2_rnm_rmnc.gff3

Finally, another filtered annotation file is produced where I identify and extract only the longest isoform for each gene

	python3 ${scriptsDir}/gff3_longest.py -gff Machtx_SR1_v2_rnm_rmnc.gff3 > Machtx_SR1_v2_rmnc_longest.gff3

##### 4.6.2 Mapping: SL RNA-seq reads <a name="SL_reads_mapping"></a>
Here I identify and map just the reads that contain SL sequence (SLreads)  
	
	genomeDir=${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/final_assembly
	cd ${genomeDir}/final_assembly_softmasked
	
	rnaDir=${HOME}/data/mac_genomes/genome_Machtx/canu_assembly/rna_mapping_v2/SL_trimmed_reads/cutadapt
	readFiles=$(echo $(ls ${rnaDir}/*_adTrSLTr_R2-SLreads.fq) | sed 's; ;,;g')

	ulimit -n 2048 # Need to set a higher limit to the maximum number of open files (originally at 1024)

	STAR --runThreadN 25 --genomeDir ./ --readFilesIn ${readFiles} --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0.5 --outFilterMatchNmin 0 --alignEndsType EndToEnd --outFileNamePrefix ./Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito_rptsftmsk_RNA_filt_SLreads --outSAMstrandField intronMotif --outSAMtype BAM Unsorted

	samtools sort -@ 25 -o Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito_rptsftmsk_RNA_filt_SLreadsAligned.out_srt.bam Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito_rptsftmsk_RNA_filt_SLreadsAligned.out.bam

	samtools index Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito_rptsftmsk_RNA_filt_SLreadsAligned.out_srt.bam

	bamCoverage -p 25 --bam Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito_rptsftmsk_RNA_filt_SLreadsAligned.out_srt.bam --outFileFormat bigwig --outFileName Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito_rptsftmsk_RNA_filt_SLreadsAligned.out_srt.cov.bigwig
	
	bamCoverage -p 25 --bam Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito_rptsftmsk_RNA_filt_SLreadsAligned.out_srt.bam --outFileFormat bedgraph --outFileName Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito_rptsftmsk_RNA_filt_SLreadsAligned.out_srt.cov.bedgraph

Extract the regions with coverage > 100x

	awk 'BEGIN{FS=OFS="\t"}{if ($4 > 100) print $0 }' Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito_rptsftmsk_RNA_filt_SLreadsAligned.out_srt.cov.bedgraph > Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito_rptsftmsk_RNA_filt_SLreadsAligned.out_srt.cov_min100.bedgraph
	
	bedtools merge -d 10 -i Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito_rptsftmsk_RNA_filt_SLreadsAligned.out_srt.cov_min100.bedgraph > Machtx_SR1_curated_l2_m22_h190_a60_arrow2-pilon5_nomito_rptsftmsk_RNA_filt_SLreadsAligned.out_srt.cov_min100.merged.bedgraph

	ulimit -n 1024 # Re-set the limit to the maximum number of open files to the default

Find the overlapping between the bedgraph of SL reads and the annotation .gff

	python3 ${scriptsDir}/gff2coding_region_bed.py -gff Machtx_SR1_v2_rnm.gff3 > coding_annotation.bed
	bedtools intersect -wo -a coding_annotation.bed -b Machtx_SR1_v2_rptsftmsk_RNA_filt_SLreadsAligned.out_srt.cov_min100.merged.bedgraph > overlaps.tab