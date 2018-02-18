#!/usr/bin/env python3

"""
> create_go_annots_sprot_trembl.py <

Based on the blastp results of the transcripts vs. sprot/trembl,
create a GO annotation file for the transcripts.

Assumes input files follow the standard pattern of "xxxx_vs_sprot.tGO.tsv" and
"xxxx_vs_trembl.tGO.tsv".

Hit threshold = 1e-5

Workflow goes:
1. If there's a hit from Swiss-Prot below the threshold, take it.
2. If not, check TrEMBL.
3. nr doesn't have GO annotations... too bad.

Program uses the tabular form of the blastp results (parse_blast_xml.py's
--table), not the raw XML format.
"""
import argparse
import re

import natural_sort
import parse_go_obo

parser = argparse.ArgumentParser(description="""
Based on the blastp results of the transcripts vs. sprot/trembl,
create a GO annotation file for the transcripts.""")

parser.add_argument('species', metavar="species_code",
                    help="3/4-letter code for species in question.")
parser.add_argument('-n', '--no_parents', action='store_true',
                    help='suppresses the search for parent GO terms.')
args = parser.parse_args()

E_VALUE_THRESHOLD = 1e-5

sprot_tsv = open('{}_vs_sprot.tGO.tsv'.format(args.species))
trembl_tsv = open('{}_vs_trembl.tGO.tsv'.format(args.species))

# NOTE: UniProtKB IDs can now be in 6-char or 10-char format. script _can_ deal
#       with both. create dictionary to store UniProtKB ID of the homologue that is
# the most similar to the transcripts
#   transcript_uniprot_id[transcript_name] = 'XXXXXX'
transcript_uniprot_id = {}

for line in sprot_tsv:
    # exclude header line
    if 'sp|' in line:
        cols = line.split('\t')
        transcript_uniprot_id[cols[0]] = re.search('sp\|(\w+)\|', cols[1]).group(1)

for line in trembl_tsv:
    # exclude header line
    if 'tr|' in line:
        cols = line.split('\t')
        if cols[0] not in transcript_uniprot_id:
            transcript_uniprot_id[cols[0]] = re.search('tr\|(\w+)\|', cols[1]).group(1)

all_uniprot_id = set(transcript_uniprot_id.values())

goa_tsv = open('/lithium/data_repo/goa/goa_uniprot_all.parsed.gaf')
# create dictionary to translate the six-character UniProt ID into GO terms
goa_data = {}
for line in goa_tsv:
    if 'GO' in line:
        cols = line.strip().split('\t')
        if cols[0] in all_uniprot_id:
            if cols[0] not in goa_data:
                goa_data[cols[0]] = []

            goa_data[cols[0]].append(cols[1])

# search for parent GO terms for each child term
associated_go_terms = {}
for t in transcript_uniprot_id:
    expanded_go_terms = []
    # not all six-char UniProt IDs have GO terms assigned to them
    if transcript_uniprot_id[t] in goa_data:
        if args.no_parents:
            # suppress search for parent GO terms
            expanded_go_terms = goa_data[transcript_uniprot_id[t]]
        else:
            for u in goa_data[transcript_uniprot_id[t]]:
                # NOTE: if the uncommented line crashes, modify the line to
                # expanded_go_terms += parse_go_obo.get_all_parent_terms(parse_go_obo.translate_alt_id(u))
                expanded_go_terms += parse_go_obo.get_all_parent_terms(u)

    associated_go_terms[t] = sorted(list(set(expanded_go_terms)))

# to create three separate files, only with BP, CC OR MF terms.
binned_go_terms = {}
for t in associated_go_terms:
    binned_go_terms[t] = parse_go_obo.bin_go_terms_by_namespace(associated_go_terms[t])

# print it out in a list
with open('{}_go_annots.{}all.tsv'.format(args.species, 'no_parents.' if args.no_parents else ''), 'w') as f:
    for t in natural_sort.natural_sort(associated_go_terms):
        if associated_go_terms[t]:
            # only print transcripts with non-empty GO term lists
            print (t, ','.join(associated_go_terms[t]), sep='\t', file=f)

bins = ['biological_process', 'cellular_component', 'molecular_function']
for b in bins:
    abbrev_fn = ''.join([x[0] for x in b.split('_')])
    with open('{}_go_annots.{}{}.tsv'.format(args.species, 'no_parents.' if args.no_parents else '', abbrev_fn), 'w') as f:
        for t in natural_sort.natural_sort(binned_go_terms):
            if binned_go_terms[t][b]:
                # only print transcripts with non-empty GO term lists
                print (t, ','.join(binned_go_terms[t][b]), sep='\t', file=f)
