#! /bin/sh

# fastqfilter.py enumerate ../tests/ecoli-nmeth/reads-nmeth-all_2d.fasta ../tests/ecoli-nmeth/reads-nmeth-all_2d-enumerated.fasta

### Run Pbdagcon to produce th econsensus, and evaluate it with dnadiff.

contigs="../tests/ecoli-nmeth/layout_20160124_020347-ra/contigs_fast.fasta"
afg="../tests/ecoli-nmeth/layout_20160124_020347-ra/contigs.afg"
reads="../tests/ecoli-nmeth/reads-nmeth-all_2d-enumerated.fasta"
reference="../tests/ecoli-nmeth/ecoli_K12_MG1655_U00096.3.fasta"

results="results/minimus/consensus-quiver.fasta"
mappings="results/minimus/mapped"
readsafg="results/minimus/reads.afg"
bank="results/minimus/reads.bnk"
contigafg="results/minimus/contigs.afg"
contigbank="results/minimus/contigs.bnk"

mkdir -p results/minimus

export PATH=$PATH:/home/isovic/assembler/git/ra-consensus/tools/amos-3.1.0/bin/
# export PERLLIB=$PERLLIB:/home/isovic/assembler/git/ra-consensus/tools/amos-3.1.0/src/PerlModules/

# cp $afg results/minimus/
# toAmos -s $reads -o $readsafg
# mkdir -p $bank
# bank-transact -m $readsafg -b $bank
# bank-transact -m $contigafg -b $bank

bank-unlock $bank
make-consensus -v 1 -o 100 -b $bank


# make-consensus -e 0.30 -o 100 $afg $bank
# make-consensus -e 0.30 -o 100 -b $contigbank $bank


# mkdir -p $contigbank
# bank-transact -m $contigafg -b $contigbank