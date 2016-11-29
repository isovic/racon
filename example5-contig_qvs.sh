#! /bin/sh

# contig="test-data/fake_lambda_contig_with_qv/fake_contig.fastq"
contig="test-data/fake_lambda_contig_with_qv/fake_contig2.fastq"
contig="test-data/fake_lambda_contig_with_qv/fake_contig3.fastq"

tools/minimap/minimap $contig test-data/fake_lambda_contig_with_qv/reads_without_fake_contig.fastq > test-data/fake_lambda_contig_with_qv/reads_without_fake_contig.paf

bin/racon -t 4 test-data/fake_lambda_contig_with_qv/reads_without_fake_contig.fastq test-data/fake_lambda_contig_with_qv/reads_without_fake_contig.paf $contig test-data/fake_lambda_contig_with_qv/consensus-no_contig_qv.fasta
mkdir -p test-data/fake_lambda_contig_with_qv/dnadiff
rm test-data/fake_lambda_contig_with_qv/dnadiff/1.report
dnadiff -p test-data/fake_lambda_contig_with_qv/dnadiff/1 test-data/lambda/NC_001416.fa test-data/fake_lambda_contig_with_qv/consensus-no_contig_qv.fasta
grep "TotalBases" test-data/fake_lambda_contig_with_qv/dnadiff/1.report
grep "AlignedBases" test-data/fake_lambda_contig_with_qv/dnadiff/1.report
grep "AvgIdentity" test-data/fake_lambda_contig_with_qv/dnadiff/1.report

bin/racon -t 4 --use-contig-qv test-data/fake_lambda_contig_with_qv/reads_without_fake_contig.fastq test-data/fake_lambda_contig_with_qv/reads_without_fake_contig.paf $contig test-data/fake_lambda_contig_with_qv/consensus-contig_qv.fasta
mkdir -p test-data/fake_lambda_contig_with_qv/dnadiff
rm test-data/fake_lambda_contig_with_qv/dnadiff/2.report
dnadiff -p test-data/fake_lambda_contig_with_qv/dnadiff/2 test-data/lambda/NC_001416.fa test-data/fake_lambda_contig_with_qv/consensus-contig_qv.fasta
grep "TotalBases" test-data/fake_lambda_contig_with_qv/dnadiff/2.report
grep "AlignedBases" test-data/fake_lambda_contig_with_qv/dnadiff/2.report
grep "AvgIdentity" test-data/fake_lambda_contig_with_qv/dnadiff/2.report
