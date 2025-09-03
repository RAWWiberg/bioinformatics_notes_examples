# Genome assemblies of the simultaneously hermaphroditics flatworms *Macrostomum cliftonense* and *Macrostomum hystrix*

**R. Axel W. Wiberg\***, Jeremias N. Brand, Gudrun Viktorin, Jack Mitchell, Christian Beisel, Lukas Schärer

**\***Correspondence: raw.wiberg@protonmail.com

This repository contains pipelines (and all mentioned custom scripts) for the assembly and population genomic analysis of the *M. hystrix* and *M. cliftonense* genomes.  

Notes on the genome assembly, annotation, snp-calling, and comparative analysis pipelines are here:

- [Assembly, and annotation for *Macrostomum cliftonensis*](Maccli_23b-7d_pipeline.md)
- [Assembly, and annotation for *Macrostomum hystrix*](Machtx_SR1_pipeline.md)
- [Comparative analyses](Mac_comparative_pipeline.md)

In addition, there are pipelines, scripts, and results for gene expression analyses also included in this project.

- [Positional RNA-seq data analyses](./positional_rnaseq)  
See also the folder `positional_rnaseq`  
- [Adult vs. hatchling contrasts](./AvH_contrasts)  
See also the folder `AvH_contrasts`  

This repository also contains an R project file (`macrostomum_genomes.Rproj`). 
This project is associated with the R notebooks in the `Rnotebooks` folder that detail the summaries and analyses presented in the paper.
These R notebooks rely on various packages, details of which are given at the top of each notebook.


export PATH=$PATH:/data/programs/FastQC:/data/programs/samtools-1.9:/data/programs/cutadapt-1.8.3/bin:/data/programs/fastsimcoal26:/data/programs/plink1.9:/data/programs/bwa-mem2-2.2.1_x64-linux:/data/programs/gatk-4.6.2.0:/data/programs/admixture_linux-1.3.0
export JAVA_HOME=/opt/jdk17
export PATH=$JAVA_HOME/bin:$PATH
