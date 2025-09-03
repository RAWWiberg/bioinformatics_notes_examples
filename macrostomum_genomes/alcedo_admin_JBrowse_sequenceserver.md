This document describes how sequenceserver (https://sequenceserver.com) and JBrowse (https://jbrowse.org) were set up on evo-alcedo.

R. Axel W. Wiberg, 17.06.2020

####################

# 1. sequenceserver
Followed instructions from: https://sequenceserver.com/doc/
# Install sequenceserver
Sequenceserver was already installed on alcedo

$ which sequenceserver
/usr/local/bin/sequenceserver

Running sequenceserver prompts to install NCBI BLAST+ 2.2.30. I have done this. 

--------------------
## 1.1 Create a user for sequenceserver that will start an sequenceserver instance on boot.

	useradd -s /sbin/nologin seqserveruser

Copy blast binaries to seqserveruser home (as root)

	cp -r /root/.sequenceserver/ncbi-blast-2.2.30+ /home/seqserveruser/

Make sure permissions are good

	chmod go+rc /home/seqserveruser/ncbi-blast-2.2.30+

Copy the directory with blast databases to seqserveruser home (as root)

	cp -r /home/scharer_group/dbs/blast_db /home/seqserveruser/

Make sure permissions are good

	chmod go+rx /home/seqserveruser/blast_db

The database directory will be set as "/home/seqserveruser/blast_db", [!N.B. this should not change!]
[!N.B. This directory also contains the links.rb file that is needed to integrate with JBrowse!]

--------------------
## 1.2 Autostart configuration
Created an autostartup systemd service by adding the file "sequenceserver.service" to "/etc/systemd/system/"

Contents of "sequenceserver.service":

	[Unit]
	Description=SequenceServer server daemon
	Documentation="file://sequenceserver --help" "http://sequenceserver.com/doc"
	After=network.target
	
	[Service]
	Type=simple
	User=seqserveruser
	ExecStart=/usr/local/bin/sequenceserver -bin /home/seqserveruser/ncbi-blast-2.2.30+/bin -c /home/seqserveruser/.sequenceserver.conf
	KillMode=process
	Restart=on-failure
	RestartSec=42s
	RestartPreventExitStatus=255
	
	[Install]
	WantedBy=multi-user.target

To get systemd to start the sequenceserver on boot need to tell systemctl about the new service file 
[!N.B. make sure no instance of sequenceserver is running]

let systemd know about changed files

	systemctl daemon-reload

enable service for automatic start on boot

	systemctl enable sequenceserver.service

start service immediately

	systemctl start sequenceserver.service

If something doesn't work or sequenceserver is down, check the logs with: 

	journalctl -u sequenceserver.service 

the most recent entries are at the bottom of the stream

The file `.sequenceserver.conf` for sequenceserver has been placed in `/home/seqserveruser/`.
This file tells sequenceserver which port to serve data out of and also which directory to use for databases.

Contents of `.sequenceserver.conf`:

	---
	:num_threads: 1
	:port: 4567
	:host: 0.0.0.0
	:database_dir: /home/seqserveruser/blast_db


## 1.3 Make sure that the correct port (4567) is open on alcedo.
I followed this procedure (https://stackoverflow.com/questions/24729024/open-firewall-port-on-centos-7), to open port 4567.

	firewall-cmd --get-active-zones
	firewall-cmd --zone=public --add-port=2888/tcp --permanent
	firewall-cmd --reload


## 1.4 Make a blast database
To make new blast databases use the ncbi-blast-2.2.30+ package that was downloaded by sequenceserver:

For example:

	/home/seqserveruser/ncbi-blast-2.2.30+/bin/makeblastdb -in Mlig_3_7.fasta -title Mlig_3_7_genome -out Mlig_3_7_genome -dbtype nucl -parse_seqids -hash_index



# 2. JBrowse
An Apache webserver is already running on alcedo (set up by Lukas Zimmerman)

## 2.1 Download and unzip the Jbrowse package

## 2.2 Move this directory to a location that is served by Apache (/var/www/html)

	su
	mv ./JBrowse-1.16.9 /var/www/html/jbrowse

Change ownership of the directory to root

	chown `whoami` /var/www/html/jbrowse


Plugins (following instructions on: http://gmod.org/wiki/JBrowse_Configuration_Guide#Plugins)

Apollo:



Run the setup script, this checks the system and sets up the example files. Also downloads a bunch of necessary perl moduls and places them in a local folder (hopefully this doesn't mess up the system). 

	cd /var/www/html/jbrowse
	./setup.sh

## 2.3 Change CentOS permissions
Now it should work, but on Centos things are complicated. There is an extra level of security for permissions of directories (SELinux).Lukas Zimmerman changed the permissions for the directory `/var/www/html/jbrowse`

	chcon -Rv --type=httpd_sys_content_t /var/www/html/jbrowse

## 2.4 Make a new .conf file for apache to properly serve .bam, .bami, and .bai files

	su
	touch /etc/httpd/conf.d/jbrowse.conf
	nano /etc/httpd/conf.d/jbrowse.conf

Contents of "jbrowse.conf":

	AddType application/octet-stream .bam .bami .bai

Test the configurations

	apachectl configtest

N.B. this test warns about not being able to determine the servers domain name, this is probably (?) not an issue.

Reload apache

	systemctl reload httpd.service

Jbrowse is now set up!

# 3. Link sequenceserver to jbrowse

Place the file `links.rb` in the directory of databases (see 1.1)

Modify the contents of `sequenceserver.service` (in `/etc/systemd/system/`) so that sequenceserver loads this `links.rb` file

Contents of `sequenceserver.service`:

	[Unit]
	Description=SequenceServer server daemon
	Documentation="file://sequenceserver --help" "http://sequenceserver.com/doc"
	After=network.target
	[Service]
	Type=simple
	User=seqserveruser
	ExecStart=/usr/local/bin/sequenceserver -bin /home/seqserveruser/ncbi-blast-2.2.30+/bin -r /home/seqserveruser/blast_db/links.rb -c /home/seqserveruser/.sequenceserver.conf
	KillMode=process
	Restart=on-failure
	RestartSec=42s
	RestartPreventExitStatus=255
	
	[Install]
	WantedBy=multi-user.target

Restart sequenceserver

	systemctl restart sequenceserver.service


# 4. Set up a new genome for jbrowse (brief example with Mlig_3_7)
Due to issues with SELinux security protocols I have placed data directories directly in `/var/www/html/jbrowse/`

If this becomes an issue, we will need to think of something else.

	mkdir /var/www/html/jbrowse/Mlig_3_7

Copy the genome fasta file to this directory

	cp /path/to/Mlig_3_7.fasta /var/www/html/jbrowse/Mlig_3_7/
	cd /var/www/html/jbrowse/Mlig_3_7

Prepare the genome .fasta file for loading

	/var/www/html/jbrowse/bin/prepare-refseqs.pl --fasta ./Mlig_3_7.fasta

This creates a sub-directory called `data`, as well as an initial file called `trackList.json` and `tracks.conf`
These can then be modified in various ways to display data in the browse.

Copy the v. 3 annotation to the data folder

	cp /path/to/Mlig_RNA_3_7_DV1.v3.coregenes.bestORF.gff3.gz /var/www/html/jbrowse/Mlig_3_7/
	gunzip ./Mlig_RNA_3_7_DV1.v3.coregenes.bestORF.gff3.gz
	/var/www/html/jbrowse/bin/flatfiles-to-json.pl --gff ./Mlig_RNA_3_7_DV1.v3.coregenes.bestORF.gff3

This will add a track entry to the "trackList.json" file. the .json file or the .conf file can also be edited manually.

For more info on how to set up tracks and all manner of custom configuration, see [http://jbrowse.org/docs/installation.html](http://jbrowse.org/docs/installation.html)


# 5. Set up Maccli draft genome for jbrowse (example below is for Maccli)

	mkdir /var/www/html/jbrowse/Maccli
	cd /var/www/html/jbrowse/Maccli

Copy the genome fasta file to this directory

	cp ~/data/genome_Maccli/canu_assembly/polished/Maccli_genome_polished.fasta ./
	sed -i 's;_pilon;;g' Maccli_genome_polished.fasta

###############################################################################
### OPTIONAL: Set up a sequencesserver blast database for the genome
Copy the genome fasta file to the sequencesserver database directory

	cp Maccli_genome_polished.fasta /home/seqserveruser/blast_db
	cd /home/seqserveruser/blast_db
	sed -i 's;>;>Maccli-;g' Maccli_genome_polished.fasta
	/home/seqserveruser/ncbi-blast-2.2.30+/bin/makeblastdb -in Maccli_genome_polished.fasta -title Maccli_genome -out Maccli_genome -dbtype nucl -parse_seqids -hash_index

Restart sequenceserver

	su
	systemctl restart sequenceserver.service
	exit
	
	cd /var/www/html/jbrowse/Maccli
###############################################################################


## Prepare the genome .fasta file for loading to JBrowse
$ /var/www/html/jbrowse/bin/prepare-refseqs.pl --fasta ./Maccli_genome_polished.fasta


Copy the stringtie annotation to the data folder

	cp ~/data/genome_Maccli/canu_assembly/gg_transcriptome_assembly/stringtie_assembly/stringtie.gff ./
	sed -i 's;_pilon;;g' stringtie.gff
	gffread stringtie.gff > stringtie.gff3
	/var/www/html/jbrowse/bin/flatfile-to-json.pl --nameAttributes "id" --gff ./stringtie.gff3 --trackLabel stringtie_genes

	cp ~/data/genome_Maccli/canu_assembly/gg_transcriptome_assembly/stringtie_assembly/transdecoder/longest_orfs.gff3 ./
	sed -i 's;_pilon;;g' longest_orfs.gff3
	/var/www/html/jbrowse/bin/flatfile-to-json.pl --nameAttributes "id" --gff ./stringtie_transcripts.fasta.transdecoder.genome.gff3 --trackLabel stringtie_transdecoder


For adding a generic .bed file

	/var/www/html/jbrowse/bin/flatfile-to-json.pl --bed ./stringtie.bed --trackLabel [LABEL]

Generate searchable names from the genes in the .gff3 files

	../bin/generate-names.pl

This will add a track entry to the `trackList.json` file. 
The `.json` file or the .conf file can also be edited manually to add tracks, e.g:

Copy the RNA-seq read mappings to the data folder

	cd ./data
	
	cp ~/data/genome_Maccli/canu_assembly/rna_mapping/polished/Maccli_genome_polished_Mac084_2930_adSLTr_filtAligned.sortedByCoord.out.* ./
	
	samtools view -h ./Maccli_genome_polished_Mac084_2930_adSLTr_filtAligned.sortedByCoord.out.bam | sed 's;_pilon;;g' | samtools view -b - > ./Maccli_genome_polished_Mac084_2930_adSLTr_filtAligned.sortedByCoord.out_mod.bam

	rm ./Maccli_genome_polished_Mac084_2930_adSLTr_filtAligned.sortedByCoord.out.bam*x`x`
	
	samtools index Maccli_genome_polished_Mac084_2930_adSLTr_filtAligned.sortedByCoord.out_mod.bam
	
	cp ~/data/genome_Maccli/canu_assembly/rna_mapping/polished/Maccli_genome_polished_Mac084_2911_adSLTr_filtAligned.sortedByCoord.out.* ./
	
	samtools view -h ./Maccli_genome_polished_Mac084_2911_adSLTr_filtAligned.sortedByCoord.out.bam | sed 's;_pilon;;g' | samtools view -b - > ./Maccli_genome_polished_Mac084_2911_adSLTr_filtAligned.sortedByCoord.out_mod.bam
	
	rm ./Maccli_genome_polished_Mac084_2911_adSLTr_filtAligned.sortedByCoord.out.bam*
	
	samtools index Maccli_genome_polished_Mac084_291_adSLTr_filtAligned.sortedByCoord.out_mod.bam
	
	cd ../

Modify the tracks.conf file manually (according to: http://jbrowse.org/docs/faq_setup.html)

	nano ./data/tracks.conf

Contents of `tracks.conf`:

	[general]
	dataset_id = maccli

	[tracks.alignments] 
	storeClass = JBrowse/Store/SeqFeature/BAM 
	urlTemplate = ./Maccli_genome_polished_Mac084_2930_adSLTr_filtAligned.sortedByCoord.out_mod.bam 
	type = Alignments2 
	key = RNA-seq (Sample 2930)

	[tracks.alignments] 
	storeClass = JBrowse/Store/SeqFeature/BAM 
	urlTemplate = ./Maccli_genome_polished_Mac084_2911_adSLTr_filtAligned.sortedByCoord.out_mod.bam 
	type = Alignments2 
	key = RNA-seq (Sample 2911)


	[tracks.mytrack] storeClass=JBrowse/Store/SeqFeature/BAM urlTemplate=Maccli_genome_polished_Mac084_2930_adSLTr_filtAligned.sortedByCoord.out_mod.bam type=Alignments2 key=RNA-seq (Sample 2930)



To delete a track and remove the associated data directory

	../bin/remove-track.pl -trackLabel stringtie_longest_orfs --delete


Add dataset to jbrowse.conf

	[datasets.maccli]
	url = ?data=Maccli%2Fdata
	name = Maccli Genome v1
