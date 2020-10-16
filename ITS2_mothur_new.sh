#!/bin/bash

# Changed by CJB to accomodate new python scripts that can work on gzipped sequence files and 
# sequence files nested in sample folders without decompressing or making unnecessary copies of the data 8/14/19
# Script for processing 16s samples and running Mothur on them:
# Usage:  /Volumes/GriffenLeysLab/zShared/Scripts_dbs/./ITS2_mothur_new.sh <1. sample folder name> <2. file name> <3.processor cores>

# Standardized Process for Fungal ITS Processing:

# make folders

#mkdir fastq

# Unzip fastq files from sequencing facility recursively:

#find $1/  -name "*.gz" | grep -v Apple | while read file; do gunzip "$file"; done

# Move Fastq files to folder:

# mv `find $1 -type f -name "*.fastq.gz" | grep -v Apple` fastq

#Change: Removed # removed what? Ran (1/4/17)

python /Volumes/GriffenLeysLab/zShared/Scripts_dbs/make_file_new.py $1 $2.txt


# Processing in Mothur

mothur "#make.contigs(file=$2.txt,processors=$3); summary.seqs();screen.seqs(minlength=150, maxambig=10); trim.seqs(oligos=/Volumes/GriffenLeysLab/zShared/Scripts_dbs/fungal.ITS2.oligos.txt, pdiffs=2); summary.seqs(); summary.seqs()"

# Remove bdiffs/pdiffs from fasta file:

cat $2.trim.contigs.good.trim.fasta | sed 's/bdiffs.*//g' > $2.trimmed.fasta

# Filter by Mean quality score of 28

python /Volumes/GriffenLeysLab/zShared/Scripts_dbs/mean_qual_screen_new.py $2.txt $2.trimmed.fasta $2.contigs.groups $2.q28.fasta $2.q28.groups


# BLAST with ISHAM

mkdir blast 

mv $2.q28.fasta blast/

cd blast

/Volumes/GriffenLeysLab/zShared/Scripts_dbs/./blast.bsh $2.q28.fasta /Volumes/GriffenLeysLab/Daniel/database/Fungal_db/ISHAM_ITS_db/ITSSeq.clean.fas $3 100

cd ..


#/Volumes/GriffenLeysLab/zShared/Scripts_dbs/./load_blast.sh $2
