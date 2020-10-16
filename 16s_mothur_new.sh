#!/bin/bash

# Changed by CJB to accomodate new python scripts that can work on gzipped sequence files and 
# sequence files nested in sample folders without decompressing or making unnecessary copies of the data 8/14/19 and 1/6/2020
# Script for processing 16s samples and running Mothur and BLAST on them:
# Usage:  /Volumes/GriffenLeysLab/zShared/Scripts_dbs/./16s_mothur_new.sh <1. sample folder name> <2. file name> <3.processor cores>

# Standardized Process for Bacterial 16S Processing:

# make folders (lines commented out by CB, uses unaltered compressed fastq as input)

#mkdir fastq

# Unzip fastq files from sequencing facility recursively:

#find $1/  -name "*.gz" | grep -v Apple | while read file; do gunzip "$file"; done

# Move Fastq files to folder:

# mv `find $1 -type f -name "*.fastq.gz" | grep -v Apple` fastq

#Change: Removed # removed what? Ran (1/4/17)

python /Volumes/GriffenLeysLab/zShared/Scripts_dbs/make_file_new.py $1 $2.txt


# Processing in Mothur

mothur "#make.contigs(file=$2.txt,processors=$3); summary.seqs();screen.seqs(minlength=450, maxlength=570, maxambig=10); trim.seqs(oligos=/Volumes/GriffenLeysLab/zShared/Scripts_dbs/miseq.oligos, pdiffs=2); summary.seqs(); summary.seqs()"

# Remove bdiffs/pdiffs from fasta file:

cat $2.trim.contigs.good.trim.fasta | sed 's/bdiffs.*//g' > $2.trimmed.fasta

# Filter by Mean quality score of 28

python /Volumes/GriffenLeysLab/zShared/Scripts_dbs/mean_qual_screen_new.py $2.txt $2.trimmed.fasta $2.contigs.groups $2.q28.fasta $2.q28.groups


# BLAST with CORE

mkdir blast 

mv $2.q28.fasta blast/

cd blast

/Volumes/GriffenLeysLab/zShared/Scripts_dbs/./blast.bsh $2.q28.fasta /Volumes/GriffenLeysLab/zShared/Scripts_dbs/core_vag_fm_2020_01_06.fasta $3 400

# Uncomment line below to load results into SQL database
# /Volumes/GriffenLeysLab/zShared/Scripts_dbs/./load_blast.sh $2
# Uncomment and modify line below to generate text file that is readable in R:
# python /Volumes/GriffenLeysLab/zShared/Scripts_dbs/./sum_to_taxa4.py $2.q28.fasta.sum /Volumes/GriffenLeysLab/zShared/Scripts_dbs/core_vag_fm_taxonomy_2020_01_06.tab ../RS_odont.q28.groups RS_odont_blast_table.txt
