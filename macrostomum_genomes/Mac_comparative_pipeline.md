# Comparative analyses

### Table of contents
0. [Install packages](#0packages)
1. [Collect CDS/AA sequences](#collectcds)
	1. [Prepare Maclig annotation file + genes](#prepare_maclig)
	2. [Get CDS and AA sequences for genes annotated in all species](#get_all_cds)
2. [Proteinortho+PoFF ortholog identification](#2proteinortho)

This document details example code used for the identifying of orthologs and building of orthogroups.
For clarity, I have ommitted, in most cases, full paths to input and output files.

## 0. Install packages, make environments, obtain docker containers <a name="0packages"></a>

**proteinortho**

	docker pull quay.io/biocontainers/proteinortho:6.0.28--hfd40d39_0
	docker pull quay.io/biocontainers/proteinortho:6.0.28--hb0e25da_1 #newest

	conda create -n proteinortho
	conda activate proteinortho
	conda install proteinortho
	conda deactivate
	conda list --export -n proteinortho > proteinortho_env.txt

My own helper scripts that are used throughout can be found in the folder `scripts`
Where these were used I prefix with a variable as a placeholder for the path to these scripts: `${scriptsDir}`

## 1. Collect CDS/AA sequences <a name="collectcds"></a>  
### 1.1. Pepare Maclig annotation file + genes  <a name="prepare_maclig"></a>  
The Maclig annotation file `Mlig_RNA_3_7_DV1.v3.coregenes.bestORF.gff3` contains several genes with no predicted coding sequence (open reading frame).  
These seem to be all the `tbone` entries which, if they contain an ORF, will also have a `transdecoder` entry. Thus I remove all `tbone` entries.

	sed '/.*tbone.*/d' Mlig_RNA_3_7_DV1.v3.coregenes.bestORF.gff3 > Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_notbone.gff3
	python3 ${scriptsDir}/gff3_longest.py -gff Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_notbone.gff3 > Mlig_RNA_3_7_DV1.v3.coregenes_longest.gff3
	python3 ${scriptsDir}/gff_mod.py -gff Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest.gff3 > Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest_mod.gff3
	sed -i '/scaf0/d' Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest_mod.gff3


### 1.2. Get CDS and AA sequences for genes annotated in all species <a name="get_all_cds"></a>  
Subset .gff3 file to longest isoform per gene

	python3 ${scriptsDir}/gff3_longest.py -gff Maccli_GV23d_v2_rnm.gff3 > Maccli_v2_longest.gff3
	
	python3 ${scriptsDir}/gff3_longest.py -gff Machtx_SR1_v2_rnm.gff3 > Machtx_v2_longest.gff3

Collect the CDS sequences for longest isoform of each gene.

	maccli_genome=Maccli_GV23d_v2.rptsftmsk.fasta
	maccli_annotation=Maccli_GV23d_v2_longest.gff3
	
	gffread $maccli_annotation -x Maccli_v2_longest.cds -g $maccli_genome

	machtx_genome=Machtx_SR1_v2.rptsftmsk.fasta
	machtx_annotation=Machtx_SR1_v2_longest.gff3
	
	gffread $machtx_annotation -x Machtx_v2_longest.cds -g $machtx_genome

	maclig_genome=Mlig_3_7_rptsftmsk.fasta
	maclig_annotation=Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest_mod.gff3
	
	gffread $maclig_annotation -x Mlig_RNA_3_7_DV1.v3.coregenes_longest.cds -g $maclig_genome

Collect the aa sequences for each gene

	cd /home/axel/data/mac_genomes/orthologs/proteinortho/aa

	gffread $machtx_annotation -y Machtx_v2_longest.pep -g $machtx_genome
	sed -i 's;\.t[[:digit:]]*$;;g' Machtx_v2_longest.pep
	
	gffread $maccli_annotation -y Maccli_v2_longest.pep -g $maccli_genome
	sed -i 's;\.t[[:digit:]]*$;;g' Maccli_v2_longest.pep
	
	gffread $maclig_annotation -y Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest.pep -g $maclig_genome

Turn gff3 files from pasa into gff files for Proteinortho+PoFF.
This is done with the `gff2gene_order.py` script.
These gene order files need to have the same names as the .pep files, except for the suffix (.gff3).

	cp $maccli_annotation ./
	sort -k1,1 -k4,4n ./Maccli_GV23d_v2_longest.gff3 > ./Maccli_v2_longest_srt.gff3
	awk '{print $1}' ./Maccli_GV23d_v2_longest.gff3 | sort -k1,1 | uniq > ./Maccli_v2_longest_contigs.list
	python3 ${scriptsDir}/gff2gene_order.py -gff ./Maccli_v2_longest_srt.gff3 -chroms ./Maccli_v2_longest_contigs.list > ./Maccli_v2_longest_geneorder.gff3
	cp ./Maccli_v2_longest_geneorder.gff3 Maccli_v2_longest.gff3

	cp $machtx_annotation ./
	sort -k1,1 -k4,4n ./Machtx_SR1__v2_longest.gff3 > ./Machtx_v2_longest_srt.gff3
	awk '{print $1}' ./Machtx_SR1_v2_longest.gff3 | sort -k1,1 | uniq > ./Machtx_v2_longest_contigs.list
	python3 ${scriptsDir}/gff2gene_order.py -gff ./Machtx_v2_longest_srt.gff3 -chroms ./Machtx_v2_longest_contigs.list > ./Machtx_v2_longest_geneorder.gff3
	cp ./Machtx_v2_longest_geneorder.gff3 ./aa/Machtx_v2_longest.gff3

	cp $maclig_annotation ./Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest.gff3
	sort -k1,1 -k4,4n ./Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest.gff3 > ./Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest_srt.gff3
	awk '{print $1}' ./Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest_srt.gff3 | sort -k1,1 | uniq > ./Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest_contigs.list
	python3 ${scriptsDir}/gff2gene_order.py -gff ./Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest_srt.gff3 -chroms ./Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest_contigs.list > ./Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest_geneorder.gff3
	cp ./Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest_geneorder.gff3 ./aa/Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest.gff3


### 2. Proteinortho+PoFF ortholog identification <a name="2proteinortho"></a>

	conda activate proteinortho

	proteinortho -keep -force -cpus=30 -dups=3 -cs=5 -alpha=0.7 -project=mac_orths_dups3-cs5-alpha0.7 -synteny Machtx_v2_longest.pep Maccli_v2_longest.pep Mlig_RNA_3_7_DV1.v3.coregenes.bestORF_longest.pep
	
	proteinortho_history.pl -project=mac_orths 
	proteinortho_history.pl -project=mac_orths [] -step=1
	proteinortho_history.pl -project=mac_orths [] -step=2
	proteinortho_history.pl -project=mac_orths [] -step=3

	conda deactivate