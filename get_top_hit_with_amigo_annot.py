#!/usr/bin/env python3

"""
> get_top_hit_with_amigo_annot.py <

Solves the annoying problem where the top hit to Swiss-Prot / TrEMBL might not
have an annotation, but the subsequent few hits do. This script takes in a 
parse_blast_xml.py --table file, and for each gene/transcript, retain the 
best hit with an AmiGO annotation.

This file can then be fed into create_go_annots and create_top_hit to
produce files that contain "best hit that has an annotation".
"""

import argparse
import csv
import re
import sys
import time

parser = argparse.ArgumentParser(description="""
Solves the annoying problem where the top hit to Swiss-Prot / TrEMBL might not
have an annotation, but the subsequent few hits do. This script takes in a 
parse_blast_xml.py --table file, and for each gene/transcript, retain the 
best hit with an AmiGO annotation.""")

parser.add_argument('blast_tsv', metavar="tsv_filename",
                    type=argparse.FileType('r'),
                    help="parsed BLAST XML in .tsv format.")
args = parser.parse_args()

start_time = time.time()
print ('start!', file=sys.stderr)

uniprot_ids = {}
# read data
tsv_reader = csv.reader(args.blast_tsv, delimiter='\t')

# skip header row
next(tsv_reader)

for row in tsv_reader:
    if row:     # skips empty rows
        if row[0] not in uniprot_ids:
            uniprot_ids[row[0]] = []
        
        id = re.search('[sp|tr]\|(\w+)\|', row[1]).group(1)
        uniprot_ids[row[0]].append(id)

# create list of UniProtKB IDs with GO terms
uniprot_ids_with_go = []

goa_tsv = '/lithium/data_repo/goa/goa_uniprot_all.unique_ids.txt'
with open(goa_tsv) as f:
    for line in f:
        line = line.strip()
        
        if line:     # skips empty rows
            uniprot_ids_with_go.append(line)

print ('finish go reading!', time.time() - start_time, file=sys.stderr)

# create dictionary to store best hit with a GO term
best_hit_with_go = {}
temp = 0
for x in uniprot_ids:
    temp += 1
    for y in uniprot_ids[x]:
        if y in uniprot_ids_with_go:
            best_hit_with_go[x] = y
            print ('found best hit for {}! [{}/{}]'.format(x, len(best_hit_with_go), temp), time.time() - start_time, file=sys.stderr)
            break

print ('finish best hit dictionary!', time.time() - start_time, file=sys.stderr)

# re-read data, and print lines that correspond to best hit that has at least
# a GO term
args.blast_tsv.seek(0)
tsv_reader = csv.reader(args.blast_tsv, delimiter='\t')

# preserve header row
print ('\t'.join(next(tsv_reader)))

for row in tsv_reader:
    if row:     # skips empty rows
        if row[0] in best_hit_with_go:
            if re.search('[sp|tr]\|{}\|'.format(best_hit_with_go[row[0]]), row[1]):
                print ('\t'.join(row))