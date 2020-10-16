#This program generates a fasta and group file with 28% mean quality (phred) score
#Parameters: 1. sequence_info_file 2. processed_fasta_file 3. groups_file 4. output_file 5. output_groups
import sys
from Bio import SeqIO
from mimetypes import guess_type
import gzip
meanquals = {}
with open(sys.argv[1], 'r') as infile:
    fastqfiles = [(line.rstrip().split('\t')[1], line.rstrip().split('\t')[2]) for line in infile]
for (forfile, revfile) in fastqfiles:
    if guess_type(forfile)[1] == 'gzip':
        forhandle = gzip.open(forfile, 'r')
        revhandle = gzip.open(forfile, 'r')
    else:
        forhandle = open(forfile, 'r')
        revhandle = open(revfile, 'r')
    for (forrec, revrec) in zip(SeqIO.parse(forhandle, 'fastq'), SeqIO.parse(revhandle, 'fastq')):
        if forrec.id == revrec.id:
            quals = forrec.letter_annotations['phred_quality'] + revrec.letter_annotations['phred_quality']
            meanqual = float(sum(quals))/float(len(quals))
            meanquals[forrec.id] = meanqual
        else:
            print 'unmatched pairs for.id {1} rev.id {2} file {3}'.format(forrec.id, revrec.id, forfile)
            break
            
with open('meanquals.txt', 'w') as outfile:
    for seqid in meanquals:
        outfile.write(seqid + '\t' + str(meanquals[seqid]) + '\n')
        
total = 0
kept = 0
with open(sys.argv[2], 'r') as infasta, open(sys.argv[3], 'r') as ingroups, open(sys.argv[4], 'w') as outfasta, open(sys.argv[5], 'w') as outgroups:
    groupdict = {}
    for line1 in ingroups:
        linelist1 = line1.rstrip().split('\t')
        groupdict[linelist1[0]] = linelist1[1]
    for line2 in infasta:
        if line2.startswith('>'):
            lastid = line2.split('\t')[0] + '\n'
            total += 1
        elif lastid is not None:
            key = lastid.rstrip()[1:].replace('_', ':')
            if meanquals[key] >= 28.0:
                kept += 1
                outfasta.write(lastid)
                outfasta.write(line2)
                outgroups.write('\t'.join((lastid.rstrip()[1:],
groupdict[lastid.rstrip()[1:]])) + '\n')
print 'Input sequences: ', total
print 'Kept sequences: ', kept
