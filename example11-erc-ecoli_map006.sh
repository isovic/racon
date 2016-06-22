#! /bin/sh

# make modules
# make tools
make -j

threads=4

mkdir -p results/temp
mkdir -p temp

### Lambda:
layout_gfa=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/layout.gfa
reads=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/reads.fastq
reference=test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/ecoli_K12_MG1655_U00096.3.fasta
dataset=ecoli_map006
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
	tools/minimap/minimap -t ${threads} -Sw5 -L100 -m0 $readsfasta $readsfasta > $paf
echo "tools/miniasm/misc/paf2mhap.pl $contigs $readsfasta $paf > $mhap"
scripts/paf2mhap.pl $readsfasta $readsfasta $paf > $mhap
echo $reads
echo $mhap

echo "Running Racon:"
/usr/bin/time --format "Command line: %C\nReal time: %e s\nCPU time: -1.0 s\nUser time: %U s\nSystem time: %S s\nMaximum RSS: %M kB\nExit status: %x" --quiet -o ${memtime_racon} \
	bin/racon -M 5 -X -4 -G -8 -E -6 --bq 10 --mhap --erc -t ${threads} ${reads} ${mhap} ${reads} ${consensus}
echo "Racon exited."

rawreadssam=results/temp/consensus-${dataset}-${suffix}-$curriter.rawreads.ref.sam
# tools/graphmap/bin/graphmap-not_release align -a anchor --rebuild-index -B 0 -b 3 -r ${reference} -d ${reads} -o ${rawreadssam} --extcigar -t ${threads}
echo "Raw reads statistics:"
# codebase/samscripts/src/errorrates.py base ${reference} ${rawreadssam}

ercreadssam=results/temp/consensus-${dataset}-${suffix}-$curriter.rawreads.erc.sam
tools/graphmap/bin/graphmap-not_release align -a anchor --rebuild-index -B 0 -b 3 -r ${reference} -d ${consensus} -o ${ercreadssam} --extcigar -t ${threads}
echo "Error-corrected reads statistics:"
codebase/samscripts/src/errorrates.py base ${reference} ${ercreadssam}

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













# [02:30:00 Run] Automatically setting the maximum allowed number of regions: max. 500, attempt to reduce after 0
# [02:30:00 Run] No limit to the maximum number of seed hits will be set in region selection.
# [02:30:00 Run] Reference genome is assumed to be linear.
# [02:30:00 Run] Only one alignment will be reported per mapped read.
# [02:30:00 ProcessReads] All reads will be loaded in memory.
# [02:30:03 ProcessReads] All reads loaded in 3.02 sec (size around 480 MB). (248550679 bases)
# [02:30:03 ProcessReads] Memory consumption: [currentRSS = 825 MB, peakRSS = 825 MB]
# [02:30:03 ProcessReads] Using 4 threads.
# [02:35:21 ProcessReads] [CPU time: 1221.32 sec, RSS: 846 MB] Read: 25483/25483 (100.00%) [m: 25029, u: 454]                                                        
# [02:35:21 ProcessReads] Memory consumption: [currentRSS = 846 MB, peakRSS = 872 MB]

# [02:35:21 ProcessReads] All reads processed in 1221.38 sec (or 20.36 CPU min).
# Raw reads statistics:
# Processing alignment: 25000...
#                       	mean	std	median	min	max
# Error rate:     	14.70%	5.70%	13.15%	6.13%	51.23%
# Insertion rate: 	6.52%	2.53%	6.00%	1.79%	30.70%
# Deletion rate:  	3.97%	1.89%	3.49%	1.15%	42.00%
# Mismatch rate:  	4.21%	2.47%	3.52%	1.02%	22.23%
# Match rate:     	89.27%	4.55%	90.45%	57.47%	96.67%
# Read length:    	9759.33	4249.78	9346.00	227.00	58229.00
# Difference ratio:	29:44:27 (mismatch:insertion:deletion)
# [02:37:53 Index] Running in normal (parsimonious) mode. Only one index will be used.
# [02:37:53 Index] Index already exists. Loading from file.
# [02:37:59 Index] Index loaded in 2.38 sec.
# [02:37:59 Index] Memory consumption: [currentRSS = 346 MB, peakRSS = 478 MB]

# [02:37:59 Run] Automatically setting the maximum allowed number of regions: max. 500, attempt to reduce after 0
# [02:37:59 Run] No limit to the maximum number of seed hits will be set in region selection.
# [02:37:59 Run] Reference genome is assumed to be linear.
# [02:37:59 Run] Only one alignment will be reported per mapped read.
# [02:37:59 ProcessReads] All reads will be loaded in memory.
# [02:37:59 ProcessReads] All reads loaded in 0.43 sec (size around 63 MB). (65439783 bases)
# [02:37:59 ProcessReads] Memory consumption: [currentRSS = 402 MB, peakRSS = 478 MB]
# [02:37:59 ProcessReads] Using 4 threads.
# [02:39:22 ProcessReads] [CPU time: 330.78 sec, RSS: 417 MB] Read: 5311/5311 (100.00%) [m: 5310, u: 1]                                                              
# [02:39:22 ProcessReads] Memory consumption: [currentRSS = 417 MB, peakRSS = 478 MB]

# [02:39:22 ProcessReads] All reads processed in 330.79 sec (or 5.51 CPU min).
# Error-corrected reads statistics:
# Processing alignment: 5000...
#                       	mean	std	median	min	max
# Error rate:     	9.51%	4.20%	10.06%	0.74%	39.06%
# Insertion rate: 	4.47%	2.11%	4.51%	0.31%	29.52%
# Deletion rate:  	2.64%	1.27%	2.60%	0.32%	36.00%
# Mismatch rate:  	2.41%	1.38%	2.42%	0.02%	15.73%
# Match rate:     	93.13%	3.29%	92.93%	68.60%	99.63%
# Read length:    	12321.73	4205.46	11700.00	1127.00	45107.00
# Difference ratio:	27:46:27 (mismatch:insertion:deletion)
