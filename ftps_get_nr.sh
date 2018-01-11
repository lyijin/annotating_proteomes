#!/bin/bash

for a in {00..78}
do
    echo "Grabbing nr.${a}.tar.gz"
    ssh liewy@ftps.kaust.edu.sa "wget -q -O - ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.${a}.tar.gz" >> nr.${a}.tar.gz
done
