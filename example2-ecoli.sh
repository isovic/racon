#! /bin/sh

make modules
make tools
make -j

### E. Coli:
awk '$1 ~/S/ {print ">"$2"\n"$3}' test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/layout.gfa > test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/layout.fasta
contigs=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/layout.fasta
reads=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/reads.fastq
sam=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/alignments.sam
dataset=ecoli
# msa=mafft
# msa=poav2
msa=poa
consensus=temp/consensus-${dataset}-${msa}.fasta
memtime=temp/consensus-${dataset}-${msa}.memtime
reference=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/ecoli_K12_MG1655_U00096.3.fasta
tools/graphmap/bin/Linux-x64/graphmap align -a anchor -z 0 -c 40 -B 0 -r ${contigs} -d ${reads} -o ${sam} --extcigar
mkdir -p temp
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o $memtime \
	bin/consise --align 1 -M 1 -X -1 -G -1 -E -1 -w 500 --ovl-margin 0.15 --msa ${msa} -b 200 -t 4 --winpath temp/window.fasta ${contigs} ${sam} ${consensus}
	# bin/consise --align 1 -M 1 -X -1 -G -1 -E -1 -w 500 --msa ${msa} -b 200 -t 4 --winpath temp/window.fasta ${contigs} ${sam} ${consensus}
mkdir -p temp/dnadiff-${dataset}
rm temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
dnadiff -p temp/dnadiff-${dataset}/consise-${dataset}-${msa} ${reference} ${consensus}
grep "TotalBases" temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
grep "AlignedBases" temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
grep "AvgIdentity" temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
cat $memtime
