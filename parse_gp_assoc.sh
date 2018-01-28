#!/bin/bash

# deflate file temporarily
gzip -d goa_uniprot_all.gaf.gz -k

# NOTE: it has come to my attention that the newer goa_uniprot_all.gaf.gz files 
# have more lines starting with a !. my original code assumed 8 lines; it's 
# 11 now. in order to future-proof the script, i'll write something that 
# will work, but will take far longer to run
# Regardless, the idea is to remove lines with !, and take columns 2 and 5
grep -v '^!' goa_uniprot_all.gaf | cut -f 2,5 > goa_uniprot_all.parsed.gaf

# remove the non-compressed file
rm -f goa_uniprot_all.gaf

# parse out unique IDs into another file. why is this done? while my python
# script could technically get unique IDs by reading goa_uniprot_all.parsed.gaf,
# in practice, it took freaking ages
# This command aims to produce a file containing ONLY one UniProt ID per line, 
# i.e.
#   A0A000
#   A0A001
#   ...
cut -f 1 goa_uniprot_all.parsed.gaf | uniq > goa_uniprot_all.unique_ids.txt
