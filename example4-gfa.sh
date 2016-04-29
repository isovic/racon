#! /bin/sh

# make modules
# make tools
# make -j

### Lambda:
gfa=test-data/lambda/layout-miniasm.gfa
# gfa=/home/isovic/work/eclipse-workspace/git/consise/codebase/seqlib/sample-data/test.gfa
reads=test-data/lambda/reads.fastq
consensus=temp/gfacons.fasta

scripts/gfacons.py $gfa $reads $consensus
