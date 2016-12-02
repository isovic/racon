#! /bin/sh

# make modules
# make tools
make -j

threads=4
threads=1

mkdir -p results/temp
mkdir -p temp

### Lambda:
layout_gfa=test-data/lambda/layout-miniasm.gfa
reads=test-data/lambda/reads.fastq
reference=test-data/lambda/NC_001416.fa
dataset=lambda_30x_ont
suffix=paf-overlap

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
# mhap=results/temp/consensus-${dataset}-${suffix}-$curriter.mhap
consensus=results/consensus-${dataset}-${suffix}-$curriter.fasta
memtime_minimap=results/consensus-${dataset}-${suffix}-$curriter.minimap.memtime
memtime_racon=results/consensus-${dataset}-${suffix}-$curriter.racon.memtime

# Convert PAF to MHAP:
# echo "tools/minimap/minimap -Sw5 -L100 -m0 $contigs $readsfasta > $paf"
# tools/minimap/minimap -Sw5 -L100 -m0 $contigs $readsfasta > $paf
echo "tools/minimap/minimap $contigs $readsfasta > $paf"
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${memtime_minimap} \
	tools/minimap/minimap $contigs $readsfasta > $paf
# echo "tools/miniasm/misc/paf2mhap.pl $contigs $readsfasta $paf > $mhap"
# scripts/paf2mhap.pl $contigs $readsfasta $paf > $mhap
echo $reads
echo $mhap

echo "Running Racon:"
echo "    bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 -t ${threads} ${reads} ${paf} ${contigs} ${consensus}"
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${memtime_racon} \
 	bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 --ovl-margin 0.1 -t ${threads} ${reads} ${paf} ${contigs} ${consensus}
	# bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 -t ${threads} --mhap --reads $reads ${contigs} ${sam} ${consensus}
	# bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 -t 1 --num-batches 1 --start-window 0 --winbatch 1 ${contigs} ${sam} ${consensus}
echo "Racon exited."

############################################
### Run dnadiff to get the Avg. Identity ###
mkdir -p results/temp/dnadiff-${dataset}
rm results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix}.report
dnadiff -p results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix} ${reference} ${consensus}
grep "TotalBases" results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix}.report
grep "AlignedBases" results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix}.report
grep "AvgIdentity" results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix}.report
cat $memtime_minimap
cat $memtime_racon
############################################

# exit

# Edit distance calculation - Avg. Identity doesn't take deletions into account ###
echo ""
echo "Evaluating the results."
scripts/edcontigs.py ${reference} ${consensus}
##########################################



# #### Run the second iteration ###
previter=iter1
curriter=iter2
contigs=results/consensus-${dataset}-${suffix}-$previter.fasta
sam=results/temp/consensus-${dataset}-${suffix}-$curriter.sam
paf=results/temp/consensus-${dataset}-${suffix}-$curriter.paf
# mhap=results/temp/consensus-${dataset}-${suffix}-$curriter.mhap
consensus=results/consensus-${dataset}-${suffix}-$curriter.fasta
memtime_minimap=results/consensus-${dataset}-${suffix}-$curriter.minimap.memtime
memtime_racon=results/consensus-${dataset}-${suffix}-$curriter.racon.memtime

# # Convert PAF to MHAP:
# # echo "tools/minimap/minimap -Sw5 -L100 -m0 $contigs $readsfasta > $paf"
# # tools/minimap/minimap -Sw5 -L100 -m0 $contigs $readsfasta > $paf
echo "tools/minimap/minimap $contigs $readsfasta > $paf"
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${memtime_minimap} \
	tools/minimap/minimap $contigs $readsfasta > $paf
# echo "tools/miniasm/misc/paf2mhap.pl $contigs $readsfasta $paf > $mhap"
# scripts/paf2mhap.pl $contigs $readsfasta $paf > $mhap
echo $reads
echo $mhap

echo "Running Racon:"
echo "    bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 -t ${threads} ${reads} ${paf} ${contigs} ${consensus}"
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${memtime_racon} \
	bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 --ovl-margin 0.1 -t ${threads} ${reads} ${paf} ${contigs} ${consensus}
	# bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 -t 1 --num-batches 1 --start-window 0 --winbatch 1 ${contigs} ${sam} ${consensus}
echo "Racon exited."
############################################
### Run dnadiff to get the Avg. Identity ###
mkdir -p results/temp/dnadiff-${dataset}
rm results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix}.report
dnadiff -p results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix} ${reference} ${consensus}
grep "TotalBases" results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix}.report
grep "AlignedBases" results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix}.report
grep "AvgIdentity" results/temp/dnadiff-${dataset}/racon-${dataset}-${suffix}.report
cat $memtime_minimap
cat $memtime_racon
############################################

## Edit distance calculation - Avg. Identity doesn't take deletions into account ###
echo ""
echo "Evaluating the results."
scripts/edcontigs.py ${reference} ${consensus}
###########################################
