#!/bin/bash

for a in {00..79}
do
    echo "Grabbing nr.${a}.tar.gz"
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.${a}.tar.gz
done
