#! /bin/sh

# make modules
# make tools
make -j

threads=4

mkdir -p results/temp
mkdir -p temp

### Lambda:
awk '$1 ~/S/ {print ">"$2"\n"$3}' test-data/lambda/layout-miniasm.gfa > test-data/lambda/layout-miniasm.fasta
reads=test-data/lambda/reads.fastq

dataset=lambda_30x_ont
suffix=poa

reference=test-data/lambda/NC_001416.fa
memtime=results/consensus-${dataset}-${suffix}.memtime

### Run the first iteration ###
contigs=test-data/lambda/layout-miniasm.fasta
sam=results/temp/consensus-${dataset}-${suffix}-iter1.sam
consensus=results/consensus-${dataset}-${suffix}-iter1.fasta

tools/graphmap/bin/Linux-x64/graphmap align -a anchor --rebuild-index -B 0 -r ${contigs} -d ${reads} -o ${sam} --extcigar -t ${threads}
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o $memtime \
	bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 -t ${threads} ${contigs} ${sam} ${consensus}
############################################

#### Run the second iteration ###
contigs=results/consensus-${dataset}-${suffix}-iter1.fasta
sam=results/temp/consensus-${dataset}-${suffix}-iter2.sam
consensus=results/consensus-${dataset}-${suffix}-iter2.fasta

tools/graphmap/bin/Linux-x64/graphmap align -a anchor --rebuild-index -B 0 -r ${contigs} -d ${reads} -o ${sam} --extcigar -t ${threads}
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o $memtime \
	bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 -t ${threads} ${contigs} ${sam} ${consensus}
############################################

### Run dnadiff to get the Avg. Identity ###
mkdir -p results/temp/dnadiff-${dataset}
rm results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix}.report
dnadiff -p results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix} ${reference} ${consensus}
grep "TotalBases" results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix}.report
grep "AlignedBases" results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix}.report
grep "AvgIdentity" results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix}.report
cat $memtime
############################################

### Edit distance calculation - Avg. Identity doesn't take deletions into account ###
echo ""
echo "Evaluating the results."
scripts/edcontigs.py ${reference} ${consensus}
############################################
