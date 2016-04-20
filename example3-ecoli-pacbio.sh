#! /bin/sh

make modules
make tools
make -j

### E. Coli:
folder=test-data/DATASETS_FOR_CONSENSUS/ecoli_pacbio
awk '$1 ~/S/ {print ">"$2"\n"$3}' $folder/layout.gfa > $folder/layout.fasta
contigs=$folder/layout.fasta
reads=$folder/reads.fastq
sam=$folder/alignments.sam
reference=$folder/ecoli_K12_MG1655_U00096.3.fasta
# sam=temp/alignments.sam
dataset=ecoli-pacbio
msa=mafft
# msa=poav2
consensus=temp/consensus-${dataset}-${msa}.fasta
memtime=temp/consensus-${dataset}-${msa}.memtime
tools/graphmap/bin/Linux-x64/graphmap align -a anchor -z 0 -c 40 -B 0 -r ${contigs} -d ${reads} -o ${sam} --extcigar
mkdir -p temp
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o $memtime \
	bin/consise -w 500 --msa ${msa} -b 200 -t 16 --minnewseq 0.80 --maxovl 0.01 --winpath temp/window.fasta ${contigs} ${sam} ${consensus}
mkdir -p temp/dnadiff-${dataset}
rm temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
dnadiff -p temp/dnadiff-${dataset}/consise-${dataset}-${msa} ${reference} ${consensus}
grep "AlignedBases" temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
grep "AvgIdentity" temp/dnadiff-${dataset}/consise-${dataset}-${msa}.report
cat $memtime
