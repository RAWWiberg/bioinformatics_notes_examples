# Assembly, and annotation for *Macrostomum cliftonensis*

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

First make .fastq files from the PacBio .bam files

	PBreadsDir="${HOME}/data/mac_genomes/genome_Maccli/raw_reads/gDNA/PacBio/r54273_20200330_131313"
	samtools fastq -@ 10 -0 ${PBreadsDir}/1_A03/m54273_200330_132306.subreads.fastq ${PBreadsDir}/1_A03/m54273_200330_132306.subreads.bam
	samtools fastq -@ 10 -0 ${PBreadsDir}/2_C03/m54273_200331_093736.subreads.fastq ${PBreadsDir}/2_C03/m54273_200331_093736.subreads.bam
	
	samtools fasta -@ 10 -0 ${PBreadsDir}/1_A03/m54273_200330_132306.subreads.fasta ${PBreadsDir}/1_A03/m54273_200330_132306.subreads.bam
	samtools fasta -@ 10 -0 ${PBreadsDir}/2_C03/m54273_200331_093736.subreads.fasta ${PBreadsDir}/2_C03/m54273_200331_093736.subreads.bam

Now run canu, including read-correction steps.

	canu -d ./ -p Maccli_23b-7d genomeSize=230m maxThreads=50 stopOnLowCoverage=10 \
	-pacbio-raw \
	${PBreadsDir}/1_A03/m54273_200330_132306.subreads.fastq \
	${PBreadsDir}/2_C03/m54273_200331_093736.subreads.fastq


Collect N50, assembly size, etc.

	assembly-stats *.fasta # this software calculates N50 based on the total assembly size.

Run BUSCO

	conda activate busco
	busco -i Maccli_23b-7d.contigs.fasta  --cpu 10 --out Maccli_23b-7d.contigs.busco -m genome -l metazoa_odb10
	conda deactivate


#### 1.2 Purge haplotigs <a name="purge_hapl"></a>

Because the initial assemblies were much larger than expect (~2x), a likely scenario is that we have assembled both haplotypes as separate contigs. 
The following steps identify these redundant contigs and filters out the shorter of the two haplotypes.

Map original PacBio reads to the assembly

	PBreadsDir="${HOME}/data/mac_genomes/genome_Maccli/raw_reads/gDNA/PacBio/r54273_20200330_131313"
	minimap2 -t 10 --secondary=no -a -x map-pb Maccli_23b-7d.contigs.fasta ${PBreadsDir}/1_A03/m54273_200330_132306.subreads.fastq ${PBreadsDir}/2_C03/m54273_200331_093736.subreads.fastq | samtools sort -m 1G -@ 5 -o Maccli_23b-7d.contigs_readsmapped_nomult.bam -T tmp.ali
	samtools sort Maccli_23b-7d.contigs_readsmapped_nomult.bam > Maccli_23b-7d.contigs_readsmapped_nomult_srt.bam

Create the coverage histogram

	purge_haplotigs hist -g Maccli_23b-7d.contigs.fasta -b Maccli_23b-7d.contigs_readsmapped_nomult.bam -t 4
	purge_haplotigs cov -i Maccli_23b-7d.contigs_readsmapped_nomult.bam.gencov -l 5 -m 55 -h 190
	purge_haplotigs purge -g Maccli_23b-7d.contigs.fasta -c coverage_stats.csv -b Maccli_23b-7d.contigs_readsmapped_nomult.bam -t 4 -a 60

	mv tmp_purge_haplotigs/MISC/Maccli_23b-7d.contigs_readsmapped_nomult.bam.histogram.csv ./

N50, assembly size, etc.from canu report

	assembly-stats *.fasta # this software calculates N50 based on the total assembly size.


Run BUSCO

	conda activate busco
	busco -i curated_l5_m55_h190_a60.fasta --cpu 10 --out Maccli_23b-7d.purge_haplotigs_curated_l5_m55_h190_a60.busco -m genome -l  metazoa_odb10
	conda deactivate

Run blobtools

First blast contigs agains the ncbi nt and uniprot databases

	conda activate blobtools

	blastn -task megablast -query Maccli_23b-7d.contigs.fasta -db /home/scharer_group/dbs/ncbi/nt/nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 > Maccli_23b-7d.contigs.blastn-nt.out
	diamond blastx --threads 10 --query Maccli_23b-7d.contigs.fasta --db /home/scharer_group/dbs/uniprot/uniprot_ref_proteomes.diamond.dmnd --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 > Maccli_23b-7d.contigs.diamond-blastx-uniprot.out

Summarise results with blobtools

	blobtools create -i Maccli_23b-7d.contigs.fasta -o Maccli_23b-7d.contigs.blob --title Maccli_23b-7d.contigs.blob -b ../Maccli_23b-7d.contigs_readsmapped_srt.bam -t Maccli_23b-7d.contigs.blastn-nt.out -t Maccli_23b-7d.contigs.diamond-blastx-uniprot.out -x bestsumorder --db ~/packages/blobtools/data/nodesDB.txt

	blobtools view -i Maccli_23b-7d.contigs.blob.blobDB.json -x bestsumorder

	blobtools plot -i Maccli_23b-7d.contigs.blob.blobDB.json

	conda deactivate


#### 1.3 Arrow polishing <a name="arrow"></a>
Here I conduct several rounds of "polishing" with the raw PacBio reads.

Evaluate initial assembly with quast

	quast.py -o ./quast_assessment --eukaryote --gene-finding --bam Maccli_23b-7d_curated_l5_m55_h190_a60_pbmm2readsmapped_srt.bam curated_l5_m55_h190_a60.fasta

Start polishing with `arrow`

	conda activate pb-assembly

	PBreadsDir="${HOME}/data/mac_genomes/genome_Maccli/raw_reads/gDNA/PacBio/r54273_20200330_131313"
	ls -1 ${PBreadsDir}/*/*subreads.bam > bams.fofn
	pbmm2 align --sort -J 8 -j 20 curated_l5_m55_h190_a60.fasta bams.fofn Maccli_23b-7d_curated_l5_m55_h190_a60_pbmm2readsmapped_srt.bam
	pbindex Maccli_23b-7d_curated_l5_m55_h190_a60_pbmm2readsmapped_srt.bam

Polish w/ Arrow (1st iteration)

	gcpp --algorithm arrow -j 20 -o Maccli_23b-7d_curated_l5_m55_h190_a60_arrow1.fasta --reference curated_l5_m55_h190_a60.fasta ./Maccli_23b-7d_curated_l5_m55_h190_a60_pbmm2readsmapped_srt.bam

	pbmm2 align --sort -J 8 -j 15 ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow1.fasta bams.fofn Maccli_23b-7d_curated_l5_m55_h190_a60_arrow1_pbmm2readsmapped_srt.bam

	pbindex Maccli_23b-7d_curated_l5_m55_h190_a60_arrow1_pbmm2readsmapped_srt.bam

Polish w/ Arrow (x iteration)

	iter=2
	x=$((${iter}-1))

	gcpp --algorithm arrow -j 40 -o Maccli_23b-7d_curated_l5_m55_h190_a60_arrow${iter}.fasta --reference Maccli_23b-7d_curated_l5_m55_h190_a60_arrow${iter}.fasta Maccli_23b-7d_curated_l5_m55_h190_a60_arrow${x}_pbmm2readsmapped_srt.bam

	pbmm2 align --sort -J 8 -j 20 ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow${iter}.fasta bams.fofn Maccli_23b-7d_curated_l5_m55_h190_a60_arrow${iter}_pbmm2readsmapped_srt.bam
	pbindex Maccli_23b-7d_curated_l5_m55_h190_a60_arrow${iter}_pbmm2readsmapped_srt.bam

	conda deactivate

	conda activate busco

	busco -f -i Maccli_23b-7d_curated_l5_m55_h190_a60_arrow${iter}.fasta --cpu 20 --out Maccli_23b-7d_curated_l5_m55_h190_a60_arrow${iter}.busco -m genome -l metazoa_odb10

	conda deactivate

Evaluate all assemblies with quast

	quast.py -t 20 -o ../quast_assessment --eukaryote --gene-finding \
	--bam ./Maccli_23b-7d_curated_l5_m55_h190_a60_pbmm2readsmapped_srt.bam,./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow1_pbmm2readsmapped_srt.bam,./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2_pbmm2readsmapped_srt.bam,./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow3_pbmm2readsmapped_srt.bam,./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow4_pbmm2readsmapped_srt.bam \
	../purge_haplotigs/curated_l5_m55_h190_a60.fasta ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow1.fasta ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2.fasta ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow3.fasta ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow4.fasta


#### 1.4 Pilon polishing <a name="pilon"></a>
Here I conduct several rounds of "polishing" with Illumina short read data.

Trim adapters from illumina reads

	IlluminareadsDir="/home/axel/data/mac_genomes/genome_Maccli/raw_reads/gDNA/Illumina"

	adapters=${HOME}/packages/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa

	trimmomatic PE -threads 10 -phred33 -trimlog ./BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002.log -summary ./BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002.summary ${IlluminareadsDir}/BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R1_001_MM_1.fastq.gz ${IlluminareadsDir}/BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R2_001_MM_1.fastq.gz ./BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R1_001_MM_1_adTr.fq.gz ./BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R1_001_MM_1_unpaired_adTr.fq.gz ./BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R2_001_MM_1_adTr.fq.gz ./BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R2_001_MM_1_unpaired_adTr.fq.gz ILLUMINACLIP:${adapters}:2:30:10

	trimmomatic PE -threads 10 -phred33 -trimlog ./BSSE_QGF_152430_HTWN3DRXX_2_D130_131_TCTCTACT_CGCGGTTC_S3_L002.log -summary ./BSSE_QGF_152430_HTWN3DRXX_2_D130_131_TCTCTACT_CGCGGTTC_S3_L002.summary ${IlluminareadsDir}/BSSE_QGF_152430_HTWN3DRXX_2_D130_131_TCTCTACT_CGCGGTTC_S3_L002_R1_001_MM_1.fastq.gz ${IlluminareadsDir}/BSSE_QGF_152430_HTWN3DRXX_2_D130_131_TCTCTACT_CGCGGTTC_S3_L002_R2_001_MM_1.fastq.gz ./BSSE_QGF_152430_HTWN3DRXX_2_D130_131_TCTCTACT_CGCGGTTC_S3_L002_R1_001_MM_1_adTr.fq.gz ./BSSE_QGF_152430_HTWN3DRXX_2_D130_131_TCTCTACT_CGCGGTTC_S3_L002_R1_001_MM_1_unpaired_adTr.fq.gz ./BSSE_QGF_152430_HTWN3DRXX_2_D130_131_TCTCTACT_CGCGGTTC_S3_L002_R2_001_MM_1_adTr.fq.gz ./BSSE_QGF_152430_HTWN3DRXX_2_D130_131_TCTCTACT_CGCGGTTC_S3_L002_R2_001_MM_1_unpaired_adTr.fq.gz ILLUMINACLIP:${adapters}:2:30:10

Do a quick map to check which sample (128+129 or 130+131) is most similar to the genome

	bwa index ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2.fasta

	bwa mem -t 20 ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2.fasta $trimReadsDir/BSSE_QGF_152430_HTWN3DRXX_2_D130_131_TCTCTACT_CGCGGTTC_S3_L002_R1_001_MM_1_adTr.fq.gz $trimReadsDir/BSSE_QGF_152430_HTWN3DRXX_2_D130_131_TCTCTACT_CGCGGTTC_S3_L002_R2_001_MM_1_adTr.fq.gz | samtools view -Sb | samtools sort -@ 20 > ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2_D130_131-bwaIlluminaReadsMapped-CHECK.srt.bam

	samtools flagstat ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2_D130_131-bwaIlluminaReadsMapped-CHECK.srt.bam
	samtools stats ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2_D130_131-bwaIlluminaReadsMapped-CHECK.srt.bam > ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2_D130_131-bwaIlluminaReadsMapped-CHECK.srt.stats
	grep "^SN" ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2_D130_131-bwaIlluminaReadsMapped-CHECK.srt.stats

	bwa mem -t 20 ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2.fasta $trimReadsDir/BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R1_001_MM_1_adTr.fq.gz $trimReadsDir/BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R2_001_MM_1_adTr.fq.gz | samtools view -Sb | samtools sort -@ 20 > ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2_D128_129-bwaIlluminaReadsMapped-CHECK.srt.bam

	samtools flagstat ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2_D128_129-bwaIlluminaReadsMapped-CHECK.srt.bam
	samtools stats ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2_D128_129-bwaIlluminaReadsMapped-CHECK.srt.bam > ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2_D128_129-bwaIlluminaReadsMapped-CHECK.srt.stats 
	grep "^SN" ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2_D128_129-bwaIlluminaReadsMapped-CHECK.srt.stats

Based on this, sample 128+129 (i.e. 23b-7d or GV23d) has lower error rate and lower number of mismatches.
I therefore perform `pilon` polishing with this sample.

Polishing with pilon
Use the arrow2 assembly

	cp Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2.fasta Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon0.fasta

Map reads with bwa mem

Polish with pilon (x iteration)

	iter=1 # Which iteration is this
	x=$((${iter}-1)) # Which genome to use
	bwa index ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon${x}.fasta

	bwa mem -t 20 ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon${x}.fasta BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R1_001_MM_1_adTr.fq.gz BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R2_001_MM_1_adTr.fq.gz | samtools view -Sb | samtools sort -@ 20 > Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon${x}_D128_129-bwaIlluminaReadsMapped.srt.bam

	pilon --genome ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon${x}.fasta --bam ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon${x}_D128_129-bwaIlluminaReadsMapped.srt.bam --output ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon${iter} --minmq 10 --outdir ./ --tracks --changes --threads 100

	. ${scriptsDir}/count_changes.sh Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon${iter}.changes | awk 'BEGIN{FS=":"}{print $2}'

Run busco

	conda activate busco

	for iter in {0..9}; do busco -f -i ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon${iter}.fasta --cpu 40 --out Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon${iter}.busco -m genome -l metazoa_odb10; done

	conda deactivate

Evaluate all assemblies with quast

	bams=$(echo *.bam | sed 's; ;,;g')
	bams=$(echo $bams | sed 's;Maccli;./Maccli;g')
	fastas=$(echo *pilon*fasta)
	fastas=$(echo $fastas | sed 's;Maccli;./Maccli;g')
	quast.py -t 20 -o ../quast_assessment_pilon --eukaryote --gene-finding --bam ${bams} ${fastas}


#### 1.5 Final assembly <a name="final_assembly"></a>

I get rid of a bunch of superfluous labels in the fasta headers

	sed -i 's;|.*;;g' Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5.fasta

Remove lower-case letters left over from arrow polishing runs.

	seqkit seq --upper-case ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5.fasta > Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5_cln.fasta
	mv Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5_cln.fasta Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5.fasta

**mt Genome**  
Here I identify any contigs that look like they might be the mitochondrial genome. 
I decided that I would assemble and annotate these separately (see **6. mitochondrial genome**). 
So here I first remove contigs from the final genome assembly that get hits for available COI sequences.
I then carry on with annotation and comparative analyses while also assembling the mt genome in parallel

Blast the available COI sequences against assembly, I should hit 1 contig...
	
	maccli_COI="MK690022 MK690023 MK690031 MK690029 MK690032 MK690024 MK690045 MK690030 MK690036 MK690021 MK690046"

Blast Maccli COI sequences against the final genome assembly

	makeblastdb -in Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5.fasta -dbtype nucl -title Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5 -out Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5.fasta -parse_seqids -hash_index

	blastn -db ./Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5.fasta -query Maccli_COI.fasta -outfmt 6

I get a single good hit (tig00002844), this is probably the assembled mt_genome

	seqkit grep -n -w 0 -p tig00002844 Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5.fasta > Maccli_mtgenome.fasta

Remove the contig from the assembly, carry on with annotation and comparative analyses.

	seqkit grep -v -n -w 0 -p tig00002844 Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5.fasta > Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5_nomito.fasta

Now that I have removed the mtGenome we have our final assembly. 
Lets re-name the file.

	cp Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5_nomito.fasta Maccli_GV23d_v2.fasta


Compute basic stats for each contig

	fastaStats.py -f Maccli_GV23d_v2.fasta -t p -s n > Maccli_GV23d_v2.fasta.stats



### 2. Mitochondrial genome identification and assembly <a name="mt_genome"></a>

**MITObim 2** 

I use the Maclig mt_genome as an initial bait 
I use the de novo assembly option of MITObim2

Copy over the read files

	trimReadsDir=${HOME}/data/genome_Maccli/raw_reads/gDNA/Illumina	
	cp ${trimReadsDir}/BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_*_001_MM_1_adTr.fq.gz ./	

Convert Illumina reads to interleaved format and clean to make them compatible with MITObim

	rename.sh in=BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R1_001_MM_1_adTr.fq.gz in2=BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R2_001_MM_1_adTr.fq.gz out=Maccli_23b-7d_illumina_gdna_adTr_interleaved.fq prefix=read_
	
	rm ./BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_*_001_MM_1_adTr.fq.gz
	sed -i 's; ;\/;g' Maccli_23b-7d_illumina_gdna_adTr_interleaved.fq
	sed -i 's;\(\/[[:digit:]]\):;\1;g' Maccli_23b-7d_illumina_gdna_adTr_interleaved.fq

Copy over the Maclig mt genome file. 
This is simply the contig labelled as the mtGenome extracted from the Mlig_3_7 genome assembly.

	cp Maclig_mtgenome.fasta ./

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
	data = Maccli_23b-7d_illumina_gdna_adTr_interleaved.fastq
	technology = solexa
	strain = Maccli_23b-7d	
	___________________________________________________________

The next steps are all done within the mitobim docker container
Spin up docker container. 
e.g.

	docker run --rm -it -v ${HOME}/data/mac_genomes/genome_Maccli/canu_assembly/mitochondrial_shenanigans/mitobim2:/home/data chrishah/mitobim /bin/bash

Run MIRA

	cd ./data
	mira manifest.conf

Run MITObim

	date=$(date | sed 's; ;_;g' | sed 's;:;;g')
	/home/src/scripts/MITObim.pl -start 1 -end 50 --pair -sample Maccli_23b-7d -ref Maclig_mtgenome -readpool ./Maccli_23b-7d_illumina_gdna_adTr_interleaved.fastq -maf ./initial-mapping-testpool-to-Maclig-mtgenome_assembly/initial-mapping-testpool-to-Maclig-mtgenome_d_results/initial-mapping-testpool-to-Maclig-mtgenome_out.maf &> log_$date

Make sure I'm allowed to use the files that have been created 

	chown --recursive 1010:1004 ./*
	
	exit

Copy final assembly .fasta file and clean it up. 

Remove IUPAC characters, 
Rename the contig so it is clear that it is the mtGenome
Replace 'X' characters with 'N'.

	seqkit seq -u ./iteration7/Maccli_23b-7d-Maclig_mtgenome-it7_noIUPAC.fasta | sed 's;>scaf0.*;>scaf0_mito;g' > ./iteration7_Maccli_23b-7d-Maclig_mtgenome-it7_noIUPAC.fasta
	seqkit seq -w 60 -u ./iteration7_Maccli_23b-7d-Maclig_mtgenome-it7_noIUPAC.fasta > Maccli_23b-7d_mitobim_mtgenome_noIUPAC.fasta
	sed '/^[^>]/s/X/N/g' Maccli_23b-7d_mitobim_mtgenome.fasta > Maccli_23b-7d_mitobim_mtgenome_noX.fasta

Get the reads that were used for the final iteration.

		cp ./iteration7/Maccli_23b-7d-readpool-it7.fastq ./iteration7_Maccli_23b-7d-readpool-it7.fastq


Sanity check for presence of COI genes
	
	makeblastdb -in Maccli_23b-7d_mitobim_mtgenome_noIUPAC.fasta -dbtype nucl -title Maccli_23b-7d_mitobim_mtgenome_noIUPAC -out Maccli_23b-7d_mitobim_mtgenome_noIUPAC.fasta -parse_seqids -hash_index
	blastn -db Maccli_23b-7d_mitobim_mtgenome_noIUPAC.fasta -query Maccli_COI.fasta -outfmt 6
	blat Maccli_23b-7d_mitobim_mtgenome_noIUPAC.fasta Maccli_COI.fasta Maccli_COI_Maccli_23b-7d_mitobim_mtgenome.psl

This produces hits for all sequences.

Map the caught reads from iteration7 to the unpadded assembly
	
	conda activate pb-assembly

	bwa index ./Maccli_23b-7d_mitobim_mtgenome_noIUPAC.fasta

	bwa mem -p -t 50 ./Maccli_23b-7d_mitobim_mtgenome_noIUPAC.fasta ./iteration7_Maccli_23b-7d-readpool-it7.fastq | samtools view -Sb | samtools sort -@ 50 > ./Maccli_23b-7d_mitobim_mtgenome_noIUPAC_iteration7IlluminareadsMapped.srt.bam
	samtools index Maccli_23b-7d_mitobim_mtgenome_noIUPAC_iteration7IlluminareadsMapped.srt.bam
	samtools flagstat Maccli_23b-7d_mitobim_mtgenome_noIUPAC_iteration7IlluminareadsMapped.srt.bam
	
	bamCoverage -p 50 --bam Maccli_23b-7d_mitobim_mtgenome_noIUPAC_iteration7IlluminareadsMapped.srt.bam --outFileFormat bedgraph --outFileName Maccli_23b-7d_mitobim_mtgenome_noIUPAC_iteration7IlluminareadsMapped.srt.cov.bedgraph

Map Maccli PacBio reads to the Maccli mt_genome	

	PBreadsDir="${HOME}/data/mac_genomes/genome_Maccli/raw_reads/gDNA/PacBio/r54273_20200330_131313"
	ls -1 ${PBreadsDir}/*/*subreads.fastq > fastq.fofn

	pbmm2 align --sort -J 20 -j 50 Maccli_23b-7d_mitobim_mtgenome.fasta fastq.fofn Maccli_23b-7d_mitobim_mtgenome_pbmm2readsmapped_srt.bam
	pbindex Maccli_23b-7d_mitobim_mtgenome_pbmm2readsmapped_srt.bam

	conda deactivate

Add the mitobim mitochondrial genome to the assembly

	cat Maccli_GV23d_v2.fasta Maccli_23b-7d_mitobim_mtgenome_noIUPAC.fasta > Maccli_GV23d_v2_wMito.fasta


### 3. Genome assembly [part 2] <a name="genome_assembly2"></a>
#### 3.1 Repeatmasking <a name="repeatmasking"></a>

	BuildDatabase -name maccli Maccli_GV23d_v2.fasta
	RepeatModeler -database maccli -pa 30 -LTRStruct

	RepeatClassifier -consensi maccli-families.fa

Softmasking

	RepeatMasker -pa 20 -xsmall -engine rmblast -lib maccli-families.fa.classified Maccli_GV23d_v2.fasta	

Hardmasking

	RepeatMasker -pa 20 -gff -engine rmblast -lib maccli-families.fa.classified Maccli_GV23d_v2.fasta


#### 3.2 Re-mapping reads to final genome assembly <a name="re_mapping"></a>

Use the `nomito` or `nomito_rptsftmsk` genomes for this.

Re-map the Illumina reads

	trimReadsDir=${HOME}/data/mac_genomes/genome_Maccli/raw_reads/gDNA/Illumina
	bwa index Maccli_GV23d_v2.fasta

	bwa mem -t 50 Maccli_GV23d_v2.fasta ${trimReadsDir}/BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R1_001_MM_1_adTr.fq.gz ${trimReadsDir}/BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_R2_001_MM_1_adTr.fq.gz | samtools view -Sb | samtools sort -@ 50 > Maccli_GV23d_v2_D128_129-bwaIlluminaReadsMapped.srt.bam

	samtools index Maccli_GV23d_v2_D128_129-bwaIlluminaReadsMapped.srt.bam

	bwa mem -t 50 Maccli_GV23d_v2.fasta ${trimReadsDir}/BSSE_QGF_152430_HTWN3DRXX_2_D130_131_TCTCTACT_CGCGGTTC_S3_L002_R1_001_MM_1_adTr.fq.gz ${trimReadsDir}/BSSE_QGF_152430_HTWN3DRXX_2_D130_131_TCTCTACT_CGCGGTTC_S3_L002_R2_001_MM_1_adTr.fq.gz | samtools view -Sb | samtools sort -@ 50 > Maccli_GV23d_v2_D130_131-bwaIlluminaReadsMapped.srt.bam

	samtools index Maccli_GV23d_v2_D130_131-bwaIlluminaReadsMapped.srt.bam

Re-map PacBio reads

	PBreadsDir="${HOME}/data/mac_genomes/genome_Maccli/raw_reads/gDNA/PacBio/r54273_20200330_131313"
	ls -1 ${PBreadsDir}/*/*subreads.bam > bams.fofn

	pbmm2 align --sort -J 20 -j 50 Maccli_GV23d_v2.fasta bams.fofn ./Maccli_GV23d_v2_pbmm2readsmapped_srt.bam

	pbindex Maccli_GV23d_v2_pbmm2readsmapped_srt.bam


#### 3.3 QC the final and initial assemblies using merqury <a name="merqury"></a>

First get k-mers from the reads

	sh $MERQURY/best_k.sh 230000000

Here I use the reads from the strain that was sequenced, and also used for polishing.

	meryl k=19 count ${trimReadsDir}/BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002_*_001_MM_1_adTr.fq.gz output ${trimReadsDir}/BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002.meryl
	
	./merqury.sh ${trimReadsDir}/BSSE_QGF_152429_HTWN3DRXX_2_D128_129_GACCTGAA_TTGGTGAG_S2_L002.meryl Maccli_23b-7d.contigs.fasta Maccli_GV23d_v2_wMito.fasta merq_D128_129_v_Maccli_final_wmito_and_raw.out



### 4. Annotation <a name="annotation"></a>
#### 4.1 SL sequence identification <a name="SL_sequence_id"></a>
Here I use a simple approach to identify spliced leader (SL) sequences.
This approach takes very common k-mers from previously published *de novo* transcriptome assemblies. 
These k-mers are then assembled with velvet. 
The resulting sequences look very similar to previously identified *M. lignano* SL sequences.


Collect the available transcriptome assemblies for Maccli from: 
[https://doi.org/10.5281/zenodo.4543289](https://doi.org/10.5281/zenodo.4543289)
I use the `noPrimer` data

	sample=2930
	python3 ${scriptsDir}/100bpFasta.py -f ~/data/assemblies/Brand_et_al_2022_MPE/Maccli_${sample}_20180705.fa -n 100 > Maccli_${sample}_20180705.100bp.fasta
	jellyfish count -s 100M -m 19 -t 10 ./Maccli_${sample}_20180705.100bp.fasta
	jellyfish dump -c mer_counts.jf > Maccli_${sample}_20180705_mer_counts.fas
	awk 'BEGIN{FS=OFS=" "}{if($2 > 1000) print ">"$2"\n"$1}' Maccli_${sample}_20180705_mer_counts.fas > Maccli_${sample}_20180705_common_mers.fasta
	velveth ./ 10 Maccli_${sample}_20180705_common_mers.fasta -short
	velvetg ./ -cov_cutoff 1 -min_contig_lgth 10 -coverage_mask 0 -exp_cov 3
	cat contigs.fa > Maccli_${sample}_20180705_SL.fasta
	cat Maccli_${sample}_20180705_SL.fasta


#### 4.2 RNA mapping to final assembly <a name="rna_mapping"></a>

Trimming and mapping of RNA-seq reads to the genome assembly.

##### 4.2.1 Trimming: adapters <a name="adapter_trimming"></a>

This was done for many individual read files. 
Here I only give example code for a single instance.

	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	--cores 20 --info-file ${outDir}/${newRead}_adTr_cutadapt.info \
	-o ${outDir}/${newRead}_adTr_R1_001_MM_1.fastq.gz -p ${outDir}/${newRead}_adTr_R2_001_MM_1.fastq.gz \
	${read}_R1_001_MM_1.fastq.gz ${read}_R2_001_MM_1.fastq.gz

##### 4.2.2 Trimming: SL sequence <a name="SL_trimming"></a>

In stranded TruSeq libraries, the SL sequence should only be present on R2 reads. 
Nevertheless, to trim, we run cutadapt as for Illumina adapters. 
This ensures that reads which end up too short are thrown out as pairs.

This was done for many individual read files. 
Here I only give example code for a single instance.


	cutadapt -G GCCGTAAAGACGGTCTCTTACTGCGAAGACTCAATTTATTGCATG -O 6 --minimum-length 50 --cores 20 \
	-o ${READ}_adTrSLTr_R1_001_MM_1.fastq.gz -p ${READ}_adTrSLTr_R2_001_MM_1.fastq.gz \
	${READ}_adTr_R1_001_MM_1.fastq.gz ${READ}_adTr_R2_001_MM_1.fastq.gz

However, to identify which reads *contain* the SL sequence, we can run cutadapt for only the R2 reads in single-end mode and create an "info file". 
Parsing the info file allows us to then extract which reads contain the SL sequence and to map only these to the genome.

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

	STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles Maccli_GV23d_v2.rptsftmsk.fasta --runThreadN 10 -sjdbOverhang ReadLength-1 --genomeSAindexNbases 13

Map the trimmed reads

	rnaDir=${HOME}/data/genome_Maccli/canu_assembly/rna_mapping_v2/SL_trimmed_reads/cutadapt
	readFiles_R1=$(echo $(ls ${rnaDir}/*R1*fastq.gz) | grep -v "SLreads" | sed 's; ;,;g')
	readFiles_R2=$(echo ${readFiles_R1} | sed 's;_R1_001_MM_1.fastq.gz;_R2_001_MM_1.fastq.gz;g')

	ulimit -n 2048 # Need to set a higher limit to the maximum number of open files (originally at 1024)

	STAR --runThreadN 50 --genomeDir ./ --readFilesIn ${readFiles_R1} ${readFiles_R2} \
	--readFilesCommand zcat --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0.5 --outFilterMatchNmin 0 --alignEndsType EndToEnd \
	--outFileNamePrefix Maccli_GV23d_v2_rptsftmsk_filt --outSAMstrandField intronMotif --outSAMtype BAM Unsorted

	samtools sort -@ 40 -o Maccli_GV23d_v2_rptsftmsk_filtAligned.out_srt.bam Maccli_GV23d_v2_rptsftmsk_filtAligned.out.bam

	samtools index Maccli_GV23d_v2_rptsftmsk_filtAligned.out_srt.bam


#### 4.3 BRAKER <a name="braker"></a>

	braker.pl --verbosity=4 --cores 20 --gff3 --UTR=on --softmasking --species=Maccli --genome=Maccli_GV23d_v2_rptsftmsk.rptsftmsk.fasta --bam=Maccli_GV23d_v2_rptsftmsk_filtAligned.out_srt.bam

Run busco

	conda activate busco
	busco -f -i augustus.hints_utr.codingseq -c 30 -m tran --out augustus.hints_utr.codingseq-busco -l metazoa_odb10
	conda deactivate


#### 4.4 IsoSeq <a name="isoseq"></a>
These steps run the PacBio Isoseq3 pipeline as well as pos-processing steps for the IsoSeq data to produce full length expressed transcripts.

	conda activate pb-assembly

Data come from two "sources". 
The same library was sequenced twice.

	sample1 # SEQUEL I, D-BSSE
	sample2 # SEQUEL II, Lausanne

**IsoSeq3 Pipeline**

	readsDir=${HOME}/data/genome_Maccli/canu_assembly/isoseq/raw_reads

Merge sample1 and sample2 bams

	ls -1 ${readsDir}/sample*.ccs.bam > ${readsDir}/ccs_bamlist.fofn
	samtools merge -@ 10 -b ${readsDir}/ccs_bamlist.fofn ${readsDir}/sample1_2.ccs.bam
	pbindex ${readsDir}/sample1_2.ccs.bam
	ls -1 ${readsDir}/sample*.subreads.bam > ${readsDir}/subreads_bamlist.fofn
	samtools merge -@ 10 -b ${readsDir}/subreads_bamlist.fofn ${readsDir}/sample1_2.subreads.bam
	pbindex ${readsDir}/sample1_2.subreads.bam

Remove primers from ccs reads

	lima ${readsDir}/sample1_2.ccs.bam ./primers_empirical.fasta ./sample1+2/sample1_2.noprimer_empirical.ccs.bam --isoseq

Refine reads by removing poly A tails and solving chimaeric reads (sample1+2)

	isoseq3 refine -j 20 ./sample1+2/sample1_2.noprimer_empirical.ccs.IsoSeq_Express_Primer_3p--NEB_5p.bam ./primers_empirical.fasta ./sample1+2/sample1_2.noprimer_empirical.ccs.flnc.bam --require-polya

Cluster similar isoforms (sample1+2)

	isoseq3 cluster -j 20 ./sample1+2/sample1_2.noprimer_empirical.ccs.flnc.bam ./sample1+2/sample1_2.noprimer_empirical.ccs.clustered.bam --verbose --use-qvs

Polish isoforms (sample1+2)

	isoseq3 polish -j 30 ./sample1+2/sample1_2.noprimer_empirical.ccs.clustered.bam ${readsDir}/sample1_2.subreads.bam ./sample1+2/sample1_2.noprimer_empirical.ccs.polished.bam

Map the polished transcripts to the current (v2) genome assembly

	genomeDir=${HOME}/data/genome_Maccli/canu_assembly/final_assembly

	minimap2 -ax splice -t 10 -uf --secondary=no -C5 Maccli_GV23d_v2.fasta ./sample1_2.noprimer_empirical.ccs.polished.hq.fastq > Maccli_GV23d_v2_IsoSeq-sample1_2PolishedHQFLTrans.sam

Use this .sam file for cDNA_Cupcake

**cDNA_Cupcake ToFU Pipeline**

	Collapse the transcripts by their mapping locations

	sort -k 3,3 -k 4,4n Maccli_GV23d_v2_IsoSeq-sample1_2PolishedHQFLTrans.sam > Maccli_GV23d_v2_IsoSeq-sample1_2PolishedHQFLTrans.srt.sam

	gunzip sample1_2.noprimer_empirical.ccs.polished.hq.fastq.gz

	collapse_isoforms_by_sam.py -c 0.98 -i 0.95 --input sample1_2.noprimer_empirical.ccs.polished.hq.fastq --fq -s Maccli_GV23d_v2_IsoSeq-sample1_2PolishedHQFLTrans.srt.sam --dun-merge-5-shorter -o sample1_2.noprimer_empirical.ccs.polished.hq

	seqtk seq -A sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fq | sed 's;/;_;g' > sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa

	gffread sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.gff > sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.gff3

These are the final "non-redundant" transcripts/isoforms.


#### 4.5 PASA <a name="pasa"></a>
These steps first aligns IsoSeq collapsed transcripts to genome with PASA and uses them as "known" expressed sequence tags (ESTs). 
These alignments are then "updated" using the BRAKER annotation.

Copy over the sqlite config files and modify the database name

	mkdir sqlite.confs
	cp ${HOME}/packages/PASApipeline/sample_data/sqlite.confs/alignAssembly.config ./sqlite.confs/
	cp ${HOME}/packages/packages/PASApipeline/sample_data/sqlite.confs/annotCompare.config ./sqlite.confs/

In these files, change database name and path: 

	"DATABASE=/tmp/sample_mydb_pasa.sqlite" > "DATABASE=${HOME}/data/genome_Maccli/canu_assembly/pasa/Maccli_GV23d_v2_mydb_pasa.sqlite"

	cp ${HOME}/data/genome_Maccli/canu_assembly/final_assembly/Maccli_GV23d_v2.fasta ./

The PASA step usese the output from BRAKER and isoseq3 files and prepare them

**ISOSEQ3**

	cp ${HOME}/data/genome_Maccli/canu_assembly/isoseq/isoseq3/sample1+2/sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa ./	
	sed -i 's;|.*;;'g sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa
	grep ">" sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa | sed 's/>//g' > sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.names

**BRAKER**

	cp ${HOME}/data/genome_Maccli/canu_assembly/braker_annotation/augustus.hints_utr.gff3 ./

Start the `pasa` container
e.g.

	docker run --rm -it -v ${HOME}/data/genome_Maccli/canu_assembly/pasa:${HOME}/data/genome_Maccli/canu_assembly/pasa pasapipeline/pasapipeline:latest bash

The next few lines are run within the docker container:

	cd ${HOME}/data/genome_Maccli/canu_assembly/pasa  
	PASAHOME="/usr/local/src/PASApipeline"	

Align and assemble the IsoSeq data using both gmap and blat. 
This also creates a transcripts database that we will update later on.

	${PASAHOME}/Launch_PASA_pipeline.pl --replace -c sqlite.confs/alignAssembly.config -C -R --ALIGNER gmap,blat --CPU 7 -g Maccli_GV23d_v2.fasta -t sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa

	${PASAHOME}/scripts/build_comprehensive_transcriptome.dbi -c sqlite.confs/alignAssembly.config -t sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa --min_per_ID 95 --min_per_aligned 30

**PASA Update: Round 1**  
Add the BRAKER transcript file to the database and compare/update the annotation using the IsoSeq data
	
	${PASAHOME}/misc_utilities/pasa_gff3_validator.pl augustus.hints_utr.gff3
	
	${PASAHOME}/scripts/Load_Current_Gene_Annotations.dbi -c sqlite.confs/alignAssembly.config -g Maccli_23b-7d_curated_l5_m55_h190_a60_arrow2-pilon5_nomito.fasta -P augustus.hints_utr.gff3

Compare the annotations

	${PASAHOME}/Launch_PASA_pipeline.pl -c sqlite.confs/annotCompare.config -A -g Maccli_GV23d_v2.fasta -t sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa

**PASA Update: Round 2**  
Add the PASA output from the previous round to the database and compare/update the annotation using the IsoSeq data

	${PASAHOME}/misc_utilities/pasa_gff3_validator.pl Maccli_GV23d_v2_mydb_pasa.sqlite.gene_structures_post_PASA_updates.19113.gff3
	
	${PASAHOME}/scripts/Load_Current_Gene_Annotations.dbi -c sqlite.confs/alignAssembly.config -g Maccli_GV23d_v2.fasta -P Maccli_GV23d_v2_mydb_pasa.sqlite.gene_structures_post_PASA_updates.19113.gff3

Compare the annotations

	${PASAHOME}/Launch_PASA_pipeline.pl -c sqlite.confs/annotCompare.config -A -g Maccli_GV23d_v2.fasta -t sample1_2.noprimer_empirical.ccs.polished.hq.collapsed.rep.fa

Make sure I'm allowed to use the files that have been created 

	chown --recursive 1010:1004 ./*
	
	exit


#### 4.6 Final annotation <a name="final_annotation"></a>

#### 4.6.1 Annotation cleaning and collecting stats  <a name="final_annotation_stats"></a>

	cp ${HOME}/data/mac_genomes/genome_Maccli/canu_assembly/pasa/Maccli_GV23d_v2_mydb_pasa.sqlite.gene_structures_post_PASA_updates.18.gff3 ./

Remove superfluous entries and empty lines in the gff3 file

	grep -v "^#" accli_GV23d_v2_mydb_pasa.sqlite.gene_structures_post_PASA_updates.18.gff3 | sed '/^$/d' > Maccli_GV23d_v2_mydb_pasa.sqlite.gene_structures_post_PASA_updates.18_reduced.gff3
	
Make some changes to the final annotation file.  
Assures that "ID=" attributes field has a unique value for each CDS feature

	python3 ${scriptsDir}/gff_mod.py -gff Maccli_GV23d_v2_mydb_pasa.sqlite.gene_structures_post_PASA_updates.18_reduced.gff3 > Maccli_GV23d_v2.gff3

Remove all "Name=" attributes for gene and mRNA features

	sed -i 's/;Name=.*//g' Maccli_GV23d_v2.gff3

Then re-name all the genes so that they follow the current *M. lignano* convention.
i.e. Mcli[XXXX].g[YYY].t[ZZZ]

	python3 ${scriptsDir}/gff3_rename.py -gff Maccli_GV23d_v2.gff3 -prefix Mcli > Maccli_GV23d_v2_rnm.gff3

Extract cds and aa sequences

	maccli_genome=${HOME}/data/mac_genomes/genome_Maccli/canu_assembly/final_assembly/Maccli_GV23d_v2.rptsftmsk.fasta	

	gffread Maccli_GV23d_v2_rnm.gff3 -x Maccli_GV23d_v2.cds.fa -g Maccli_GV23d_v2.rptsftmsk.fasta
	gffread Maccli_GV23d_v2_rnm.gff3 -y Maccli_GV23d_v2.pep.fa -g Maccli_GV23d_v2.rptsftmsk.fasta

Run BUSCO for the final annotated transcripts

	conda activate busco

	busco -i Maccli_GV23d_v2.cds.fa --cpu 20 --out Maccli_GV23d_v2.cds.busco -m trans -l metazoa_odb10

	conda deactivate

Get some stats for the annotations.

	python3 ${scriptsDir}/gff3_stats.py -gff Maccli_GV23d_v2_rnm.gff3 -mode g > Maccli_GV23d_v2_rnm.gff3.genes.stats
	python3 ${scriptsDir}/gff3_stats.py -gff Maccli_GV23d_v2_rnm.gff3 -mode t > Maccli_GV23d_v2_rnm.gff3.transcripts.stats

Identify the genes with "invalid" coding sequences (length not multiple of 3)

	awk 'BEGIN{FS=OFS="\t"}{if ($6 == 0) print $2}' Maccli_GV23d_v2_rnm.gff3.transcripts.stats | sort -k1,1 | uniq > Maccli_GV23d_v2_rnm_invalid.list

I also produce a filtered annotation file where I have remove genes with coding-sequences whose length is not a multiple of 3.

	python3 ${scriptsDir}/gff3_remove.py -gff Maccli_GV23d_v2_rnm.gff3 -list Maccli_GV23d_v2_rnm_invalid.list > Maccli_GV23d_v2_rnm_rmnc.gff3

Finally, another filtered annotation file is produced where I identify and extract only the longest isoform for each gene

	python3 ${scriptsDir}/gff3_longest.py -gff Maccli_GV23d_v2_rnm_rmnc.gff3 > Maccli_GV23d_v2_rnm_rmnc_longest.gff3


##### 4.6.1 Mapping: SL RNA-seq reads <a name="SL_reads_mapping"></a>
Here I identify and map just the reads that contain SL sequence (SLreads)

	rnaDir=${HOME}/data/mac_genomes/genome_Maccli/canu_assembly/rna_mapping_v2/SL_trimmed_reads/cutadapt
	readFiles=$(echo $(ls ${rnaDir}/*_adTrSLTr_R2_001_MM_1-SLreads.fastq) | sed 's; ;,;g')

	ulimit -n 2048 # Need to set a higher limit to the maximum number of open files (originally at 1024)

	STAR --runThreadN 25 --genomeDir ./ --readFilesIn ${readFiles} --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0.5 --outFilterMatchNmin 0 --alignEndsType EndToEnd --outFileNamePrefix .Maccli_GV23d_v2_rptsftmsk_filt_SLreads --outSAMstrandField intronMotif --outSAMtype BAM Unsorted

	samtools sort -@ 25 -o Maccli_GV23d_v2_rptsftmsk_filt_SLreadsAligned.out_srt.bam Maccli_GV23d_v2_rptsftmsk_filt_SLreadsAligned.out.bam

	samtools index Maccli_GV23d_v2_rptsftmsk_filt_SLreadsAligned.out_srt.bam

	bamCoverage -p 25 --bam Maccli_GV23d_v2_rptsftmsk_filt_SLreadsAligned.out_srt.bam --outFileFormat bedgraph --outFileName Maccli_GV23d_v2_rptsftmsk_filt_SLreadsAligned.out_srt.cov.bedgraph

Extract the regions with coverage > 100x (also did >10x)

	awk 'BEGIN{FS=OFS="\t"}{if ($4 > 100) print $0 }' Maccli_GV23d_v2_rptsftmsk_filt_SLreadsAligned.out_srt.cov.bedgraph > Maccli_GV23d_v2_rptsftmsk_filt_SLreadsAligned.out_srt.cov_min100.bedgraph
	
	bedtools merge -d 10 -i Maccli_GV23d_v2_rptsftmsk_filt_SLreadsAligned.out_srt.cov_min100.bedgraph > Maccli_GV23d_v2_rptsftmsk_filt_SLreadsAligned.out_srt.cov_min100.merged.bedgraph
	
	ulimit -n 1024 # Re-set the limit to the maximum number of open files to the default

Find the overlapping between the bedgraph of SL reads and the annotation .gff

	python3 ${scriptsDir}/gff2coding_region_bed.py -gff Maccli_GV23d_v2_rnm.gff3 > coding_annotation.bed
	bedtools intersect -wo -a coding_annotation.bed -b Maccli_GV23d_v2_rptsftmsk_filt_SLreadsAligned.out_srt.cov_min100.merged.bedgraph > overlaps.tab
