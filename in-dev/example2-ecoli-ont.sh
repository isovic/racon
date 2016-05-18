#! /bin/sh

# make modules
# make tools
make -j

threads=12

mkdir -p results/temp
mkdir -p temp

### E. Coli MAP006:
layout_gfa=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/layout.gfa
reads=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/reads.fastq
reference=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/ecoli_K12_MG1655_U00096.3.fasta
dataset=ecoli_map006
suffix=poa

layout_fasta=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/layout.fasta
awk '$1 ~/S/ {print ">"$2"\n"$3}' ${layout_gfa} > ${layout_fasta}

### Run the first iteration ###
contigs=${layout_fasta}
sam=results/temp/consensus-${dataset}-${suffix}-iter1.sam
consensus=results/consensus-${dataset}-${suffix}-iter1.fasta
memtime=results/consensus-${dataset}-${suffix}-iter1.memtime

tools/graphmap/bin/Linux-x64/graphmap align -a anchor --rebuild-index -B 0 -r ${contigs} -d ${reads} -o ${sam} --extcigar -t ${threads}
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o $memtime \
	bin/racon -M 4 -X -5 -G -8 -E -6 --bq 10 -t ${threads} ${contigs} ${sam} ${consensus}
############################################

#### Run the second iteration ###
contigs=results/consensus-${dataset}-${suffix}-iter1.fasta
sam=results/temp/consensus-${dataset}-${suffix}-iter2.sam
consensus=results/consensus-${dataset}-${suffix}-iter2.fasta
memtime=results/consensus-${dataset}-${suffix}-iter2.memtime

tools/graphmap/bin/Linux-x64/graphmap align -a anchor --rebuild-index -B 0 -r ${contigs} -d ${reads} -o ${sam} --extcigar -t ${threads}
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o $memtime \
	bin/racon -M 4 -X -5 -G -8 -E -6 --bq 10 -t ${threads} ${contigs} ${sam} ${consensus}
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
