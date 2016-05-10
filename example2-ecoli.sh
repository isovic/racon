#! /bin/sh

make modules
make subgraph
make tools
make -j

threads=4

### E. Coli:
awk '$1 ~/S/ {print ">"$2"\n"$3}' test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/layout.gfa > test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/layout.fasta
# contigs=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/layout.fasta
reads=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/reads.fastq
# sam=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/alignments.sam
dataset=ecoli_map006
# msa=mafft
# msa=poav2
msa=poa
contigs=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/layout.gfa
sam=temp/consensus-${dataset}-${msa}-v1.sam
consensus=temp/consensus-${dataset}-${msa}-v1.fasta
#
# contigs=temp/consensus-${dataset}-${msa}-v1.fasta
# sam=temp/consensus-${dataset}-${msa}-v2.sam
# consensus=temp/consensus-${dataset}-${msa}-v2.fasta
#
# contigs=temp/consensus-${dataset}-${msa}-v2.fasta
# sam=temp/consensus-${dataset}-${msa}-v3.sam
# consensus=temp/consensus-${dataset}-${msa}-v3.fasta

mkdir -p temp
reference=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/ecoli_K12_MG1655_U00096.3.fasta
memtime=temp/consensus-${dataset}-${msa}-v1.memtime
tools/graphmap/bin/Linux-x64/graphmap align -a anchor --rebuild-index -z 0 -c 40 -B 0 -r ${contigs} -d ${reads} -o ${sam} --extcigar -t ${threads}
# tools/graphmap/bin/graphmap-not_release align -a anchor --rebuild-index -z 0 -c 40 -B 0 -r ${contigs} -d ${reads} -o ${sam} --extcigar -t ${threads} -b 3
# tools/graphmap/bin/graphmap-not_release align -a anchor --rebuild-index --mapq 3 -B 0 -r ${contigs} -d ${reads} -o ${sam} --extcigar -t ${threads}
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o $memtime \
	bin/consise --align 1 -M 5 -X -4 -G -8 -E -6 --bq 10 -w 500 --ovl-margin 0.00 --msa ${msa} -b 20000 -t ${threads} --winpath temp/window.fasta ${contigs} ${sam} ${consensus}
	# bin/consise --align 1 -M 1 -X -1 -G -2 -E -1 --bq 10 -w 100 --ovl-margin 0.00 --msa ${msa} -b 200 -t ${threads} --winpath temp/window.fasta ${contigs} ${sam} ${consensus}
mkdir -p temp/dnadiff-${dataset}
rm temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
dnadiff -p temp/dnadiff-${dataset}/consise-${dataset}-${msa} ${reference} ${consensus}
grep "TotalBases" temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
grep "AlignedBases" temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
grep "AvgIdentity" temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
cat $memtime

# tools/edlib/src/aligner ${consensus} ${reference} -p -f NICE > ${consensus}.refalign.txt
# head -n 10 ${consensus}.refalign.txt | tail -n 1
