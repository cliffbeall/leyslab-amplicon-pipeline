"""sum_to_taxa3.py - 5/22/14 Cliff Beall
takes the files from ncbiSummarize script by N. Firestone
and substitutes taxonomic an sample information so data
can be processed in R.
Major mods 8/2019 to use MiSeq data (no barcode) - intended mainly for the fungal ITS2 data at this stage
Usage: python sum_to_taxa.py <sum_folder> <tax_file> <group_file> <outfile>"""
import sys
import os

def parse_file(filename):
    ret_dict = {}
    with open(filename, 'r') as handle:
        for line in handle:
            linelist = line.rstrip().split('\t')
            ret_dict[linelist[0]] = linelist[1:]
    return ret_dict

tax_dict = parse_file(sys.argv[2])
group_dict = parse_file(sys.argv[3])

outfile = open(sys.argv[4], 'w')
outfile.write('read_id\tacc_no\tsample\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n')
total_lines = 0
total_matches = 0

filelist = [f for f in os.listdir(sys.argv[1]) if not f.startswith('.')]
for filename in filelist:
    with open(os.path.join(sys.argv[1], filename), 'r') as sum_file:
        for s_line in sum_file:
            s_list = s_line.rstrip().split('\t')
            total_lines += 1
            if float(s_list[2]) >= 0.98:
                total_matches += 1
                outfile.write('\t'.join(s_list[:2]))
                outfile.write('\t' + group_dict[s_list[0]][0] + '\t')
                outfile.write('\t'.join(tax_dict[s_list[1]]) + '\n')
outfile.close()
print "Total reads: ", total_lines
print "Matches over 98%: ", total_matches
print "Percentage w. matches: %", float(total_matches) * 100 / float(total_lines)

