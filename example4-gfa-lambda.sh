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

memtime=temp/gfacons-lambda.memtime
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o $memtime \
	scripts/gfacons.py $gfa $reads $consensus $threads $reference
mkdir -p temp/dnadiff-${dataset}
rm temp/dnadiff-${dataset}/consise-gfacons-${dataset}-${msa}.report
dnadiff -p temp/dnadiff-${dataset}/consise-gfacons-${dataset}-${msa} ${reference} ${consensus}
grep "TotalBases" temp/dnadiff-${dataset}/consise-gfacons-${dataset}-${msa}.report
grep "AlignedBases" temp/dnadiff-${dataset}/consise-gfacons-${dataset}-${msa}.report
grep "AvgIdentity" temp/dnadiff-${dataset}/consise-gfacons-${dataset}-${msa}.report
cat $memtime

tools/edlib/src/aligner ${consensus} ${reference} -p -f NICE > ${consensus}.refalign.txt
head -n 10 ${consensus}.refalign.txt | tail -n 1
