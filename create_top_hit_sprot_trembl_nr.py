#!/usr/bin/env python3

"""
> create_top_hit_sprot_trembl_nr.py <

Based on the BLASTP/BLASTX results of the proteins/transcripts vs. 
sprot/trembl/nr, create a file containing top hits from those three databases.

Requires five files, four of which has to follow naming conventions:
- xxxx_vs_sprot.tGO.tsv
- xxxx_vs_trembl.tGO.tsv
- xxxx_vs_nr.t1.tsv
- xxxx_go_annots.all.tsv

Hit threshold = 1e-5

Workflow goes:
1. If there's a hit from Swiss-Prot below the threshold, take it.
2. If not, check TrEMBL.
3. If not, then nr.

Program uses the tabular form of the BLAST results (parse_blast_xml.py),
not the raw XML format.

Uncomment the print statements to view how many annotations came from
SwissProt/TrEMBL/nr.
"""
import argparse
import re

import natural_sort
import parse_fasta

parser = argparse.ArgumentParser(description="""
Based on the blastp results of the transcripts vs. sprot/trembl,
create a GO annotation file for the transcripts.""")

parser.add_argument('species', metavar="species_code",
                    help="3/4-letter code for species in question.")
parser.add_argument('prot_file', metavar="fasta_file",
                    type=argparse.FileType('r'),
                    help="protein FASTA file of the gene models.")
parser.add_argument('-n', '--no_nr', action='store_true',
                    help='script works without nr too!')
parser.add_argument('-p', '--no_parents', action='store_true',
                    help='use the annotation file that do not contain parents.')
args = parser.parse_args()

all_transcripts = parse_fasta.get_all_sequences(args.prot_file, 'fasta')

sprot_tsv = open('{}_vs_sprot.tGO.tsv'.format(args.species))
trembl_tsv = open('{}_vs_trembl.tGO.tsv'.format(args.species))
if not args.no_nr:
    nr_tsv = open('{}_vs_nr.t1.tsv'.format(args.species))
go_tsv = open('{}_go_annots.{}all.tsv'.format(args.species, 
        'no_parents.' if args.no_parents else ''))

transcript_go_terms = {}
for line in go_tsv:
    cols = line.strip().split('\t')
    transcript_go_terms[cols[0]] = cols[1]

def get_go_terms(transcript):
    if transcript in transcript_go_terms:
        return '\t' + transcript_go_terms[transcript]
    else:
        return ''

# create dictionary to store the BLASTX details of each transcript
#   blast_details[transcript_name] = 'XXXXXX'
blast_details = {}

for line in sprot_tsv:
    # exclude header line
    if '\tsp|' in line:
        annot = line.split('\t')[0]
        blast_details[annot] = line.strip() + get_go_terms(annot)

for line in trembl_tsv:
    # exclude header line
    if '\ttr|' in line:
        annot = line.split('\t')[0]
        if annot not in blast_details:
            blast_details[annot] = line.strip() + get_go_terms(annot)

if not args.no_nr:
    for line in nr_tsv:
        # exclude header line
        if '\tgi|' in line:
            annot = line.split('\t')[0]
            if annot not in blast_details:
                blast_details[annot] = line.strip()

# assign 'no_hits' to transcripts that does not hits to Swiss-Prot/TrEMBL/nr
for a in all_transcripts:
    if a not in blast_details:
        blast_details[a] = a

# output a table that joins top hits from Swiss-Prot/TrEMBL/nr
print ('\t'.join(['Query', 'Source', 'Hit accession', 'Hit description', 
                  'Query length', 'Hit length', 'Query (start, end)', 
                  'Hit (start, end)', 'Frame', 'Max bit score', 
                  'Total bit score', 'Identity', 'Identity %', 'Coverage %', 
                  'Expect', 'GO terms']))

for t in natural_sort.natural_sort(blast_details):
    # find out the source (i.e. sprot/trembl/nr?) of the line
    if 'sp|' in blast_details[t]:
        src = 'Swiss-Prot'
    elif 'tr|' in blast_details[t]:
        src = 'TrEMBL'
    elif 'gi|' in blast_details[t]:
        src = 'nr'
    else:
        src = 'no_hit'

    cols = blast_details[t].strip().split('\t')

    if len(cols) > 0:
        output_list = [cols[0], src] + cols[1:]
    else:
        output_list = [cols[0], src]

    print ('\t'.join(output_list))
