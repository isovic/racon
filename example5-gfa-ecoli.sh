#! /bin/sh

make modules
make tools
make -j

### Lambda:
gfa=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/layout.gfa
reads=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/reads.fastq
consensus=temp/gfacons-ecoli.fasta
reference=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/ecoli_K12_MG1655_U00096.3.fasta
threads=4

scripts/gfacons.py $gfa $reads $consensus $threads $reference
