#! /bin/sh

make modules
make tools
make -j

### Lambda:
gfa=test-data/lambda/layout-miniasm.gfa
reads=test-data/lambda/reads.fastq
consensus=temp/gfacons-lambda.fasta
reference=test-data/lambda/NC_001416.fa
threads=4

scripts/gfacons.py $gfa $reads $consensus $threads $reference
