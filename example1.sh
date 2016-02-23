#! /bin/sh

# mkdir temp
# mkdir tools; cd tools; git clone https://github.com/isovic/graphmap.git && cd graphmap && make -j extcigar && cd ..

### Lambda:
# awk '$1 ~/S/ {print ">"$2"\n"$3}' test-data/lambda/layout-miniasm.gfa > test-data/lambda/layout-miniasm.fasta
# tools/graphmap/bin/Linux-x64/graphmap -a anchor -z 0 -c 40 -r test-data/lambda/layout-miniasm.fasta -d test-data/lambda/reads.fastq -o test-data/lambda/alignments.sam
bin/consise test-data/lambda/layout-miniasm.gfa test-data/lambda/alignments.sam temp/alt_contigs_lambda.sam
../../aligner-comparison/scripts/convert_to_bam.sh temp/alt_contigs_lambda

### E. Coli:
# tools/graphmap/bin/Linux-x64/graphmap -a anchor -z 0 -c 40 -r test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/layout.fasta -d test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/reads.fastq -o test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/alignments.sam
# bin/consise test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/layout.gfa test-data/DATASETS_FOR_CONSENSUS/ecoli_map006_ont/alignments.sam temp/alt_contigs_ecoli.fasta




### Nevermind this:

# grep ">" reads-lambda-R73-without_problematic.fasta > headers.txt
# sed 's/ /:/g' headers.txt > headers1.txt
# sed 's/>//g' headers1.txt > headers.txt
# ../../../samscripts/src/fastqfilter.py header headers.txt reads-1.fastq reads.fastq

# nucmer --maxmatch -c 100 -p dotplot/ecoli ../test-data/ecoli/ecoli_K12_MG1655_U00096.3.fasta alt_contigs.fasta
# show-coords -c dotplot/ecoli.delta > dotplot/ecoli.coords
# mummerplot --png -p dotplot/ecoli dotplot/ecoli.delta -R ../test-data/ecoli/ecoli_K12_MG1655_U00096.3.fasta -Q alt_contigs.fasta
