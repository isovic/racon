#! /bin/sh

# make modules
# make tools
make -j

threads=4

mkdir -p results/temp
mkdir -p temp

### Lambda:
layout_gfa=test-data/lambda/layout-miniasm.gfa
reads=test-data/lambda/reads.fastq
reference=test-data/lambda/NC_001416.fa
dataset=lambda_30x_ont
suffix=erc

layout_fasta=${layout_gfa}.fasta
awk '$1 ~/S/ {print ">"$2"\n"$3}' ${layout_gfa} > ${layout_fasta}
readsfasta=${reads}.fasta
codebase/samscripts/src/fastqfilter.py fastq2fasta $reads > $readsfasta

# tools/graphmap/bin/graphmap-not_release align -a anchor --rebuild-index -B 0 -b 3 -r ${contigs} -d ${reads} -o ${sam} --extcigar -t ${threads}

previter=iter0
cp ${layout_fasta} results/consensus-${dataset}-${suffix}-$previter.fasta

### Run the first iteration ###
previter=iter0
curriter=iter1
contigs=results/consensus-${dataset}-${suffix}-$previter.fasta
sam=results/temp/consensus-${dataset}-${suffix}-$curriter.sam
paf=results/temp/consensus-${dataset}-${suffix}-$curriter.paf
mhap=results/temp/consensus-${dataset}-${suffix}-$curriter.mhap
consensus=results/consensus-${dataset}-${suffix}-$curriter.fasta
memtime_minimap=results/consensus-${dataset}-${suffix}-$curriter.minimap.memtime
memtime_racon=results/consensus-${dataset}-${suffix}-$curriter.racon.memtime

echo "tools/minimap/minimap -Sw5 -L100 -m0 $readsfasta $readsfasta > $paf"
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${memtime_minimap} \
	tools/minimap/minimap -Sw5 -L100 -m0 $readsfasta $readsfasta > $paf
echo "tools/miniasm/misc/paf2mhap.pl $contigs $readsfasta $paf > $mhap"
# scripts/paf2mhap.pl $readsfasta $readsfasta $paf > $mhap
echo $reads
echo $paf
echo $mhap

echo "Running Racon:"
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${memtime_racon} \
	bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 --erc --ovl-margin 0.1 -t ${threads} ${reads} ${paf} ${reads} ${consensus}
echo "Racon exited."

rawreadssam=results/temp/consensus-${dataset}-${suffix}-$curriter.rawreads.ref.sam
tools/graphmap/bin/graphmap-not_release align -a anchor --rebuild-index -B 0 -b 3 -r ${reference} -d ${reads} -o ${rawreadssam} --extcigar -t ${threads}
echo "Raw reads statistics:"
codebase/samscripts/src/errorrates.py base ${reference} ${rawreadssam}

ercreadssam=results/temp/consensus-${dataset}-${suffix}-$curriter.rawreads.erc.sam
tools/graphmap/bin/graphmap-not_release align -a anchor --rebuild-index -B 0 -b 3 -r ${reference} -d ${consensus} -o ${ercreadssam} --extcigar -t ${threads}
echo "Error-corrected reads statistics:"
codebase/samscripts/src/errorrates.py base ${reference} ${ercreadssam}

cat ${memtime_racon}

# ### Run the seconditeration ###
# previter=iter1
# curriter=iter2
# reads=results/consensus-${dataset}-${suffix}-$previter.fasta
# sam=results/temp/consensus-${dataset}-${suffix}-$curriter.sam
# paf=results/temp/consensus-${dataset}-${suffix}-$curriter.paf
# mhap=results/temp/consensus-${dataset}-${suffix}-$curriter.mhap
# consensus=results/consensus-${dataset}-${suffix}-$curriter.fasta
# memtime_minimap=results/consensus-${dataset}-${suffix}-$curriter.minimap.memtime
# memtime_racon=results/consensus-${dataset}-${suffix}-$curriter.racon.memtime

# readsfasta=${reads}.fasta
# codebase/samscripts/src/fastqfilter.py fastq2fasta $reads > $readsfasta

# echo "tools/minimap/minimap -Sw5 -L100 -m0 $readsfasta $readsfasta > $paf"
# /usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${memtime_minimap} \
# 	tools/minimap/minimap -Sw5 -L100 -m0 $readsfasta $readsfasta > $paf
# echo "tools/miniasm/misc/paf2mhap.pl $contigs $readsfasta $paf > $mhap"
# scripts/paf2mhap.pl $readsfasta $readsfasta $paf > $mhap
# echo $reads
# echo $mhap

# echo "Running Racon:"
# /usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${memtime_racon} \
# 	bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 --erc -t ${threads} ${reads} ${mhap} ${reads} ${consensus}
# echo "Racon exited."

# ercreadssam=results/temp/consensus-${dataset}-${suffix}-$curriter.rawreads.erc.sam
# tools/graphmap/bin/graphmap-not_release align -a anchor --rebuild-index -B 0 -b 3 -r ${reference} -d ${consensus} -o ${ercreadssam} --extcigar -t ${threads}
# echo "Error-corrected reads statistics:"
# codebase/samscripts/src/errorrates.py base ${reference} ${ercreadssam}
