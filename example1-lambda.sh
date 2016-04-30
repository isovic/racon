#! /bin/sh

make modules
make tools
make -j

### Lambda:
# awk '$1 ~/S/ {print ">"$2"\n"$3}' test-data/lambda/layout-miniasm.gfa > test-data/lambda/layout-miniasm.fasta
# contigs=test-data/lambda/layout-miniasm.fasta
contigs=test-data/lambda/layout-miniasm.gfa
reads=test-data/lambda/reads.fastq
sam=test-data/lambda/alignments.sam
dataset=lambda
# msa=mafft
# msa=poav2
msa=poa
consensus=temp/consensus-${dataset}-${msa}.fasta
memtime=temp/consensus-${dataset}-${msa}.memtime
reference=test-data/lambda/NC_001416.fa
tools/graphmap/bin/Linux-x64/graphmap align -a anchor -z 0 -c 40 -B 0 -r ${contigs} -d ${reads} -o ${sam} --extcigar
mkdir -p temp
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o $memtime \
	bin/consise --align 1 -M 1 -X -1 -G -1 -E -1 --bq 10 -w 500 --ovl-margin 0.0 --msa ${msa} -b 200 -t 4 --winpath temp/window.fasta ${contigs} ${sam} ${consensus}
mkdir -p temp/dnadiff-${dataset}
rm temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
dnadiff -p temp/dnadiff-${dataset}/consise-${dataset}-${msa} ${reference} ${consensus}
grep "TotalBases" temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
grep "AlignedBases" temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
grep "AvgIdentity" temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
cat $memtime

tools/edlib/src/aligner ${consensus} ${reference} -p -f NICE > ${consensus}.refalign.txt
head -n 10 ${consensus}.refalign.txt | tail -n 1
