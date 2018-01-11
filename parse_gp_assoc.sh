#!/bin/bash

# deflate file temporarily
gzip -d goa_uniprot_all.gaf.gz -k

# parse file
sed 1,8d goa_uniprot_all.gaf | cut -f 2,5 > goa_uniprot_all.parsed.gaf

# remove non-compresed file
rm -f goa_uniprot_all.gaf

# parse out unique IDs
uniq goa_uniprot_all.parsed.gaf > goa_uniprot_all.unique_ids.txt
